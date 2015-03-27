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


#include <AnalysisValidate.h>

namespace GroupDIA {

} /* namespace GroupDIA */

template<class T1, class T2> bool large_first(const pair<T1, T2>& i, const pair<T1, T2>& j) {
	return (i.first > j.first);
}

void GroupDIA::AnalysisValidate::load(vector<string> filename) {
	for (auto it : filename)
		load(it);
}

void GroupDIA::AnalysisValidate::output_mprohet_result(string input_filename,
		string output_filename) {
	struct Info {
		double gscore;
		double p;
		double ide_score;
	};
	//Get the selected row_no for every id
	map<int, Info> row_no_name;
	vector<int> all_ids = mprophet.get_all_id();
	for (int id : all_ids) {
		Info info;
		int no;
		mprophet.get_final_result(id, no, info.gscore, info.p, info.ide_score);

		row_no_name.insert(make_pair(no, info));
	}

	int cur_no = 0;
	bool out_head = true;
	fstream fo(output_filename.c_str(), fstream::out);
	fstream fi(input_filename.c_str(), fstream::in);
	String line;
	getline(fi, line);
	line.trim();
	vector<String> col_names;
	int id_col_num = -1;
	line.split("\t", col_names);
	for (int i = 0; i < col_names.size(); i++) {
		if (col_names[i] == "transition_group_id")
			id_col_num = i;
	}
	if (out_head) {
		fo << line << "\t" << "g_score\tp\tnot_p\tprobability" << endl;
		out_head = false;
	}

	while (true) {
		getline(fi, line);
		line.trim();
		if (line.empty())
			break;

		auto cid = row_no_name.find(cur_no);
		if (cid != row_no_name.end()) {
			fo << line << "\t" << cid->second.gscore << "\t" << cid->second.p << "\t"
					<< 1 - cid->second.p << "\t" << cid->second.p << endl;
		}
		cur_no++;
	}
	fi.close();
	fo.close();
	mprophet = GroupProphet();
}

void GroupDIA::AnalysisValidate::load_mprohet_result(int project_size, string filename) {
	if (!File::exists(filename))
		return;

	double min_tran_p = get_gprophet_min_tran_p(filename);
	double min_spec_p = get_gprophet_min_spec_p(filename, min_tran_p);
	cout << min_tran_p << "\t" << min_spec_p << endl;
	map<string, int> item_name_place;
	fstream fi(filename, fstream::in);
	String line;
	getline(fi, line);
	vector<String> items;
	line.split("\t", items);
	for (int i = 0; i < items.size(); i++) {
		item_name_place[items[i]] = i;
	}
	getline(fi, line);
	while (!line.empty()) {
		getline(fi, line);
		if (line.empty() or fi.eof())
			break;
		vector<String> content;
		line.split("\t", content);

		double tran_p = stod(content[item_name_place["probability"]]);
		double ide_p = stod(content[item_name_place["ide_score"]]);
		int decoy = stoi(content[item_name_place["decoy"]]);
		int decoy_type = stoi(content[item_name_place["decoy_type"]]);

		if (decoy == 1)
			continue;
		if (decoy_type != TransitionIdentified::Target)
			continue;
		if (decoy == 0 and (tran_p < min_tran_p or ide_p * tran_p < min_spec_p))
			continue;

//Get info
		int id = stoi(content[item_name_place["id"]]);
		int exp_num = stoi(content[item_name_place["exp_num"]]);

		MProhetResult cur_result;
		cur_result.left_peak = stod(content[item_name_place["leftWidth"]]);
		cur_result.right_peak = stod(content[item_name_place["rightWidth"]]);
		cur_result.ide_score = stod(content[item_name_place["ide_score"]]);
		cur_result.add_intensity(content[item_name_place["aggr_Peak_Area"]]);
		cur_result.p = tran_p;

		TransitionInfo cur_info;
		cur_info.decoy_type = decoy_type;
		cur_info.accessions = content[item_name_place["accessions"]];
		cur_info.sequence = content[item_name_place["full_sequence"]];
		cur_info.charge = content[item_name_place["charge"]].toInt();
		cur_info.scan_num = content[item_name_place["scan_num"]].toInt();
		cur_info.mz = content[item_name_place["mz"]].toDouble();

		if (cur_info.charge <= 1)
			continue;

		auto it = id_to_no.find(id);
		if (it == id_to_no.end()) {
			TransitionMProhetResult t;
			t.add_result(exp_num, cur_result);
			if (t.size() > 0) {
				id_to_no.insert(pair<int, int>(id, trans_score.size()));
				trans_info.push_back(cur_info);
				trans_score.push_back(t);
			}
		} else {
			trans_score[it->second].add_result(exp_num, cur_result);
		}
	}
}

void GroupDIA::AnalysisValidate::score(int exp_size, string score_filename) {
	vector<int> all_ids;
	for (auto it = id_to_no.begin(); it != id_to_no.end(); ++it)
		all_ids.push_back(it->first);

#pragma omp parallel for
	for (int i = 0; i < all_ids.size(); i++) {
		auto it = id_to_no.find(all_ids[i]);
		if (it != id_to_no.end()) {
			TransitionMProhetResult& t = trans_score[it->second];
			int decoy_type = trans_info[it->second].decoy_type;
			auto& i = trans_info[it->second];
			t.remove_interference_ions();

#pragma omp critical
			{
				result_intensity.add_new_item();
				//Add intensity
				result_intensity.add_item("id", it->first);
				for (int i = 0; i < exp_size; i++) {
					result_intensity.add_item("exp_" + String(i), t.get_intensity(i));
				}
				result_intensity.add_item("swath_score", t.get_max_p());
				result_intensity.add_item("ide_score", t.get_ide_score());
				result_intensity.add_item("exp_size", t.size());
				result_intensity.add_item("decoy",
						decoy_type == TransitionIdentified::Decoy_Rand ? 1 : 0);
				result_intensity.add_item("target",
						decoy_type == TransitionIdentified::Target ? 1 : 0);
				result_intensity.add_item("accession", i.accessions);
				result_intensity.add_item("charge", i.charge);
				result_intensity.add_item("mz", i.mz);
				result_intensity.add_item("scan_num", i.scan_num);
				result_intensity.add_item("full_sequence", i.sequence);
				AASequence aas;
				aas.setStringSequence(i.sequence);
				result_intensity.add_item("sequence", aas.toUnmodifiedString());
				result_intensity.add_item("score", t.get_max_p() * t.get_ide_score());
			}
		}
	}
}

void GroupDIA::AnalysisValidate::calc_fdr() {
	typedef pair<double, bool> ScoreDecoy;
	vector<ScoreDecoy> exp1;

//Record score
	int size = result_intensity.get_size();
	for (int i = 0; i < size; i++) {
		int target = stoi(result_intensity.get_item(i, "target"));
		int exp_size = stoi(result_intensity.get_item(i, "exp_size"));
		double score = stod(result_intensity.get_item(i, "score"));

		bool is_decoy = (target == 1) ? false : true;
		ScoreDecoy c;
		c.first = score;
		c.second = is_decoy;
		exp1.push_back(c);
	}

//Calc FDR
	sort(exp1.begin(), exp1.end(), large_first<double, bool>);
	vector<pair<double, double> > exp1_fdr;
	if (!exp1.empty()) {
		auto* exp = &exp1;
		auto* exp_fdr = &exp1_fdr;
		int tp = 0, fp = 0;
		for (auto it = exp->begin(); it != exp->end(); ++it) {
			if (it->second)
				fp++;
			else
				tp++;
			exp_fdr->push_back(pair<double, double>(it->first, fp * 1.0 / tp));
		}
		double max = exp_fdr->rbegin()->second;
		for (auto it = exp_fdr->rbegin(); it != exp_fdr->rend(); ++it) {
			if (it->second > max)
				it->second = max;
			else
				max = it->second;
		}
	}
//Assign FDR
	map<double, double> exp1_score_fdr;
	for (auto it = exp1_fdr.begin(); it != exp1_fdr.end(); ++it)
		exp1_score_fdr.insert(*it);

	for (int i = 0; i < size; i++) {
		int exp_size = stoi(result_intensity.get_item(i, "exp_size"));
		double score = stod(result_intensity.get_item(i, "score"));

		auto* exp = &exp1_score_fdr;
		if (exp_size == 1)
			exp = &exp1_score_fdr;

		auto it = exp->find(score);
		result_intensity.set_item(i, "fdr", it->second);
	}
}

double GroupDIA::AnalysisValidate::get_gprophet_min_spec_p(string filename,
		double tran_p_cutoff) const {
	map<string, int> item_name_place;
	fstream fi(filename, fstream::in);
	String line;
	getline(fi, line);
	vector<String> items;
	line.split("\t", items);
	for (int i = 0; i < items.size(); i++) {
		item_name_place[items[i]] = i;
	}
	getline(fi, line);

	map<string, pair<double, bool> > seq_p;
	while (!line.empty()) {
		getline(fi, line);
		if (line.empty() or fi.eof())
			break;
		vector<String> content;
		line.split("\t", content);

		if (stod(content[item_name_place["probability"]]) < tran_p_cutoff)
			continue;

		bool decoy_tran = (stoi(content[item_name_place["decoy"]]) == 1);
		if (decoy_tran)
			continue;

		string seq = content[item_name_place["full_sequence"]];
		string charge = content[item_name_place["charge"]];
		bool decoy = (stoi(content[item_name_place["decoy_type"]]) != TransitionIdentified::Target);
		double p = stod(content[item_name_place["probability"]])
				* stod(content[item_name_place["ide_score"]]);

		seq += charge;
		auto it = seq_p.find(seq);
		if (it == seq_p.end()) {
			seq_p.insert(make_pair(seq, make_pair(p, decoy)));
		} else {
			if (p > it->second.first) {
				it->second.first = p;
				it->second.second = decoy;
			}
		}
	}
	fi.close();

	vector<pair<double, bool> > p_decoy;
	for (auto it = seq_p.begin(); it != seq_p.end(); ++it) {
		p_decoy.push_back(it->second);
	}

	sort(p_decoy.begin(), p_decoy.end(), large_first<double, bool>);
	vector<double> fdr(p_decoy.size());
	int fp = 0;
	for (int i = 0; i < p_decoy.size(); i++) {
		if (p_decoy[i].second)
			fp++;
		fdr[i] = ((double) fp) / ((double) p_decoy.size());
	}

	for (int i = fdr.size() - 1; i >= 0; i--) {
		if (fdr[i] <= ACCEPT_FDR)
			return p_decoy[i].first;
	}
	return 0;
}

double GroupDIA::AnalysisValidate::get_gprophet_min_tran_p(string filename) const {
	map<string, int> item_name_place;
	fstream fi(filename, fstream::in);
	String line;
	getline(fi, line);
	vector<String> items;
	line.split("\t", items);
	for (int i = 0; i < items.size(); i++) {
		item_name_place[items[i]] = i;
	}
	getline(fi, line);
	vector<pair<double, bool> > p_decoy;
	double decoy_size = 0;
	while (!line.empty()) {
		getline(fi, line);
		if (line.empty() or fi.eof())
			break;
		vector<String> content;
		line.split("\t", content);

		string seq = content[item_name_place["full_sequence"]];
		string charge = content[item_name_place["charge"]];
		bool decoy = (stoi(content[item_name_place["decoy"]]) == 1);
		double p = stod(content[item_name_place["probability"]]);

		p_decoy.push_back(make_pair(p, decoy));
		if (decoy)
			decoy_size++;
	}
	fi.close();

	sort(p_decoy.begin(), p_decoy.end(), large_first<double, bool>);
	double target_size = p_decoy.size() - target_size;
	vector<double> fdr(p_decoy.size());
	int fp = 0;
	for (int i = 0; i < p_decoy.size(); i++) {
		if (p_decoy[i].second)
			fp++;
		double real_fp = ((double) fp) / decoy_size * (decoy_size + target_size);
		double real_size = 2 * (i - fp);
		fdr[i] = real_fp / real_size;
	}

	for (int i = fdr.size() - 1; i >= 0; i--) {
		if (fdr[i] <= ACCEPT_FDR)
			return p_decoy[i].first;
	}
	return 0;
}

void GroupDIA::AnalysisValidate::merge_gprophet_input_file(const vector<string>& input_filename,
		string output_filename) const {
	fstream fo(output_filename.c_str(), fstream::out);
	bool out_head = true;
	for (auto it = input_filename.begin(); it != input_filename.end(); ++it) {
		fstream fi(it->c_str(), fstream::in);
		string line;
		getline(fi, line);
		if (out_head) {
			fo << line << "\n";
			out_head = false;
		}
		getline(fi, line);
		while (!line.empty()) {
			fo << line << "\n";
			getline(fi, line);
		}
		fi.close();
	}
	fo.close();
}

void GroupDIA::AnalysisValidate::output_intensity(string filename) {
	result_intensity.store(filename);
}

void GroupDIA::AnalysisValidate::load(string filename) {
	int max_id = 0;
	fstream fi(filename.c_str(), fstream::in);
	String line;
	getline(fi, line);
	line.trim();
	if (line.empty())
		return;
	vector<String> col_names;
	vector<int> selected_value_col_num;
	int selected_col_num = -1;
	int ide_score_col_num = -1;
	int decoy_col_num = -1;
	int id_col_num = -1;
	int gid_col_num = -1;
	int group_col_num = -1;
	int int_col_num = -1;
	line.split("\t", col_names);
	for (int i = 0; i < col_names.size(); i++) {
		if (col_names[i].find("var_") == 0) {
			selected_value_col_num.push_back(i);
		}
		if (col_names[i] == "decoy")
			decoy_col_num = i;
		if (col_names[i] == "transition_group_id")
			id_col_num = i;
		if (col_names[i] == "ide_score")
			ide_score_col_num = i;
		if (col_names[i] == "selected")
			selected_col_num = i;
		if (col_names[i] == "id")
			gid_col_num = i;
		if (col_names[i] == "ingroup")
			group_col_num = i;
		if (col_names[i] == "aggr_Peak_Area")
			int_col_num = i;
	}

	if (decoy_col_num == -1 or id_col_num == -1 or ide_score_col_num == -1 or selected_col_num == -1
			or col_names.size() <= 1) {
		cout << "Column names error!" << endl;
		throw;
	}

	int var_num = selected_value_col_num.size();

	map<string, int> map_id;
	while (true) {
		getline(fi, line);
		line.trim();
		if (line.empty())
			break;

		vector<String> cont;
		line.split("\t", cont);

		//Add values
		vector<double> cur_value;
		for (auto it : selected_value_col_num) {
			cur_value.push_back(cont[it].toDouble());
		}

		vector<String> cur_intensity_s;
		cont[int_col_num].split(",", cur_intensity_s);
		vector<double> cur_intensity(cur_intensity_s.size() - 1);
		for (int i = 0; i < cur_intensity_s.size() - 1; i++) {
			cur_intensity[i] = cur_intensity_s[i].toDouble();
		}

		//Decide type
		int cur_type;
		if (cont[decoy_col_num].toInt() == 0)
			if (cont[selected_col_num].toInt() == 1)
				cur_type = GroupProphet::Target;
			else
				cur_type = GroupProphet::Unknown;
		else
			cur_type = GroupProphet::Decoy;

		string& tr_gp_id = cont[id_col_num];
		int cur_id = max_id;
		auto cid = map_id.find(tr_gp_id);
		if (cid == map_id.end()) {
			map_id.insert(make_pair(tr_gp_id, cur_id));
			max_id++;
		} else {
			cur_id = cid->second;
		}
		mprophet.add_values(cont[gid_col_num].toInt(), cur_id, cur_value, cur_intensity, cur_type,
				cont[ide_score_col_num].toDouble(), cont[group_col_num].toInt());
	}
	fi.close();
}

GroupDIA::AnalysisValidate::AnalysisValidate(const Param& para) {
	mprophet.setParameters(para.copy("mProhet:", true));
	ACCEPT_FDR = (DoubleReal) para.getValue("max_accepted_fdr");
}

bool GroupDIA::AnalysisValidate::run_mprohet() {
	return mprophet.learn();
}

void GroupDIA::AnalysisValidate::prepare_for_pyprophet(const vector<string>& input_filename,
		const vector<string>& output_filename) const {
	vector<fstream> fo(output_filename.size());
	for (int i = 0; i < output_filename.size(); i++) {
		fo[i].open(output_filename[i].c_str(), fstream::out);
	}

	bool head_out = true;

	for (int i = 0; i < input_filename.size(); i++) {
		if (!File::exists(input_filename[i])) {
			continue;
		}

		fstream fi(input_filename[i].c_str(), fstream::in);

		String line;
		int exp_num_col_no = -1;
		getline(fi, line);
		vector<String> col;
		line.split("\t", col);
		for (int j = 0; j < col.size(); j++)
			if (col[j] == "exp_num")
				exp_num_col_no = j;
		if (head_out) {
			for (int j = 0; j < fo.size(); j++)
				fo[j] << line << "\n";
			head_out = false;
		}

		if (exp_num_col_no < 0)
			cout << "Error: not find exp_num" << endl;

		while (!line.empty()) {
			getline(fi, line);
			if (line.empty() or fi.eof())
				break;
			vector<String> cur;
			line.split("\t", cur);
			int exp = stoi(cur[exp_num_col_no]);
			fo[exp] << line << "\n";
		}
		fi.close();
	}

	for (int j = 0; j < fo.size(); j++)
		fo[j].close();
}
