/*************************************************************************
 * Group-DIA
 * Copyright (C) 2015 .
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 *************************************************************************/

#include "MergeTransitions.h"
const double MAX_FDR = 0.01;

template<class T1, class T2> bool large_first(const pair<T1, T2>& i, const pair<T1, T2>& j) {
	return (i.first > j.first);
}

namespace GroupDIA {

} /* namespace GroupDIA */

void GroupDIA::MergeTransitions::add_identification(String pepxml_filename) {

	vector<ProteinIdentification> proteins;
	vector<PeptideIdentification> peptides;

	GroupDIA::PepXMLFile pepxml_file;
	pepxml_file.load(pepxml_filename, proteins, peptides);

	//Decide FDR in PSM level
	vector<pair<double, bool> > score_decoy;
	for (auto it = peptides.begin(); it != peptides.end(); ++it) {
		auto pephit = it->getHits();
		auto jt = &(pephit[0]);
		int cur_id = it->getMetaValue("scan_num").toString().toInt();
		if (cur_id % 2 == 1)
			continue;
		bool is_decoy = true;
		if (ToolsSpectrum::is_target_identification_result(jt->getProteinAccessions())) {
			is_decoy = false;
		}
		score_decoy.push_back(make_pair(jt->getScore(), is_decoy));
	}
	sort(score_decoy.begin(), score_decoy.end(), large_first<double, bool>);
	vector<double> fp(score_decoy.size(), 0);
	if (!score_decoy[0].second)
		fp[0] = 1;
	for (int i = 1; i < fp.size(); i++) {
		fp[i] = fp[i - 1];
		if (score_decoy[i].second)
			fp[i]++;
	}

	double cutoff = score_decoy[0].first;
	for (int i = fp.size() - 1; i >= 0; i--) {
		double fdr = fp[i] / ((double) i);
		if (fdr < MAX_FDR) {
			cutoff = score_decoy[i].first;
			break;
		}
	}

	cout << "FDR = " << MAX_FDR * 100 << "% in PSM level, iProphet score:" << cutoff << endl;
	// load identification result based on cutoff
	for (auto it = peptides.begin(); it != peptides.end(); ++it) {
		auto pephit = it->getHits();
		auto jt = &(pephit[0]);
		int cur_id = it->getMetaValue("scan_num").toString().toInt();
		jt->setMetaValue("scan_num", cur_id);
		int cur_scan_num = cur_id / 2 - cur_windows_num;
		if (cur_scan_num % windows_size == 0) {
			if (jt->getScore() >= cutoff) {
				if (ToolsSpectrum::is_target_identification_result(jt->getProteinAccessions())) {
					map_scan_pep.insert(pair<int, PeptideHit>(cur_id, *jt));
				}
			}
		}
	}
	cout << "Add identification result: " << map_scan_pep.size() << endl;
}

void GroupDIA::MergeTransitions::choose_best_id() {
	//Record max_score
	for (auto it = map_scan_pep.begin(); it != map_scan_pep.end(); ++it) {
		const int& cur_scan_no = it->first;
		int cur_exp = get_exp_num_from_scan_num(cur_scan_no);

		string seq = it->second.getSequence().toString();
		seq += "_" + String(it->second.getCharge());
		if (cur_scan_no % 2 == 1)
			seq += "d";

		double score = it->second.getScore();
		auto i_in_max = max_score.find(seq);
		if (i_in_max == max_score.end()) {
			max_score.insert(
					make_pair(seq, vector<pair<double, int> >(exp_size, make_pair(-1.0, -1))));
			i_in_max = max_score.find(seq);
		}

		auto& cur_pair = i_in_max->second[cur_exp];
		double& max_score = cur_pair.first;
		if (score > max_score) {
			max_score = score;
			cur_pair.second = cur_scan_no;
		}
	}

	//Get max_score;
	for (auto it = max_score.begin(); it != max_score.end(); ++it) {
		double max = -1;
		int max_id = -1;
		for (auto jt = it->second.begin(); jt != it->second.end(); ++jt) {
			if (jt->first > max) {
				max = jt->first;
				max_id = jt->second;
			}
		}
		if (max_id != -1) {
			best_scan_num.insert(max_id);
		}
	}
}

void GroupDIA::MergeTransitions::load_transition(string input_filename, string output_filename) {
	if (!File::exists(input_filename))
		throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, input_filename);
	fstream fi(input_filename.c_str(), fstream::in | fstream::binary);
	while (!fi.eof()) {
		vector<pair<TransitionIdentified, int> > trs;
		if (get_next_identification(fi, trs)) {
			for (auto it = trs.begin(); it != trs.end(); ++it) {
				auto& tr = it->first;
				int scan_num = it->second;
				int exp_num = get_exp_num_from_scan_num(scan_num);

				string seq = tr.get_pep_seq();
				seq += "_" + String(tr.get_precursor_charge());
				if (scan_num % 2 == 1)
					seq += "d";

				if (max_score.find(seq)->second[exp_num].second != scan_num) {
					continue;
				}

				auto i_info = result_identification.find(seq);
				if (i_info == result_identification.end()) {
					result_identification.insert(make_pair(seq, PeptideInformation()));
					i_info = result_identification.find(seq);
					i_info->second.rt_no.resize(exp_size, pair<int, int>(-1, -1));
				}

				if (best_scan_num.find(scan_num) != best_scan_num.end()) {
					i_info->second.tr = tr;
				} else {
					auto& iadd = i_info->second;
					vector<double> product_mz;
					int peak_rt_start_no, peak_rt_length;

					tr.get_product_mz(product_mz);
					tr.get_precursor_rt_no(exp_num, peak_rt_start_no, peak_rt_length);

					iadd.product_mz.insert(iadd.product_mz.begin(), product_mz.begin(),
							product_mz.end());
					iadd.rt_no[exp_num] = make_pair(peak_rt_start_no, peak_rt_length);
				}
			}
		}
	}
	fi.close();

	fstream fo(output_filename.c_str(), fstream::out | fstream::binary);
	//Start merge transition
	merge_transition(fo);
	fo.close();
}

void GroupDIA::MergeTransitions::merge_transition(fstream& fo) {
	for (auto ip = result_identification.begin(); ip != result_identification.end(); ++ip) {
		auto& info = ip->second;
		auto& tr = ip->second.tr;
		tr.add_product_mz(info.product_mz);
		vector<double>().swap(info.product_mz);
//		cout << tr.get_uid() << "\t" << tr.get_exist_exp_num() << endl;
		for (int cur_exp = 0; cur_exp < exp_size; cur_exp++) {
			if (info.rt_no[cur_exp].first > 0) {
				tr.set_precursor_rt_no(cur_exp, info.rt_no[cur_exp].first,
						info.rt_no[cur_exp].second);
			}
		}

		String method = para.getValue("method");
		double identity_threshold = para.getValue("identity_threshold");
		int max_attempts = para.getValue("max_attempts");
		double mz_threshold = para.getValue("mz_threshold");
		double mz_shift = para.getValue("mz_shift");
		double similarity_threshold = para.getValue("similarity_threshold");

		TransitionIdentified tr_decoy(para_ide);
		if (method == "random") {
			tr.generate_decoy_by_simple_copy(tr_decoy);
			tr_decoy.set_uid((++cur_id) * windows_size + cur_windows_num);
			tr.set_product_by_seq(true);
			tr_decoy.set_product_by_seq(false);
			tr_decoy.load_store(fo, LoadSave::IOType::Write, LoadSave::Type::Basic_Precursor_Iden);
		} else {
			bool is_right;
			tr.set_product_by_seq(true);
			tr.generate_decoy_by_openswath_method(tr_decoy, is_right, method, identity_threshold,
					max_attempts, mz_threshold, mz_shift, similarity_threshold);
			if (is_right) {
				tr_decoy.set_product_by_seq(true);
				tr_decoy.set_uid((++cur_id) * windows_size + cur_windows_num);
				tr_decoy.load_store(fo, LoadSave::IOType::Write,
						LoadSave::Type::Basic_Precursor_Iden);
			}
		}

		tr.load_store(fo, LoadSave::IOType::Write, LoadSave::Type::Basic_Precursor_Iden);
	}
}

bool GroupDIA::MergeTransitions::get_next_identification(fstream& fi,
		vector<pair<TransitionIdentified, int> >& result_trans) {
	//Load transition from .chro file
	TransitionIdentified tran(para_ide);
	if (!tran.load_store(fi, LoadSave::IOType::Read, LoadSave::Type::Basic_Precursor)) {
		return false;
	}

	int scan_num = tran.get_scan_num();
	result_trans.clear();

//Set target
	int cur_scan_num = scan_num * 2;
	auto ipep = map_scan_pep.find(cur_scan_num);
	if (ipep != map_scan_pep.end()
			and ToolsSpectrum::is_target_identification_result(
					ipep->second.getProteinAccessions())) {
//Add target transition
		TransitionIdentified tran_id = tran;
		tran_id.set_uid((++cur_id) * windows_size + cur_windows_num);
		tran_id.set_scan_num(ipep->second.getMetaValue("scan_num").toString().toInt());
		tran_id.set_decoy(tran_id.Target);
		tran_id.set_identification_result(ipep->second, false);
		result_trans.push_back(make_pair(tran_id, cur_scan_num));
	}

//Set decoy
	cur_scan_num = scan_num * 2 + 1;
	ipep = map_scan_pep.find(cur_scan_num);

	if (ipep != map_scan_pep.end()
			and ToolsSpectrum::is_target_identification_result(
					ipep->second.getProteinAccessions())) {
//Add decoy transition
		TransitionIdentified tran_id = tran;
		tran_id.set_uid((++cur_id) * windows_size + cur_windows_num);
		tran_id.set_scan_num(ipep->second.getMetaValue("scan_num").toString().toInt());
		tran_id.set_decoy(tran_id.Decoy_Spec);
		tran_id.set_identification_result(ipep->second, false);
		result_trans.push_back(make_pair(tran_id, cur_scan_num));
	}

	if (result_trans.empty())
		return false;
	else
		return true;
}

GroupDIA::MergeTransitions::MergeTransitions(Param para, Param para_ide, int windows_size,
		int cur_windows_num, int exp_size) {
	cur_id = 0;
	this->para = para;
	this->para_ide = para_ide;
	this->windows_size = windows_size;
	this->cur_windows_num = cur_windows_num;
	this->exp_size = exp_size;
}

int GroupDIA::MergeTransitions::get_exp_num_from_scan_num(int scan_num) const {
	return ((scan_num / 2 - cur_windows_num) / windows_size) % exp_size;
}
