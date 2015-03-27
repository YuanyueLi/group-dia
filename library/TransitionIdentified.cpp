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

#include "TransitionIdentified.h"
//#define SINGLE
template<class T1, class T2> bool large_first(const pair<T1, T2>& i, const pair<T1, T2>& j) {
	return (i.first > j.first);
}

namespace GroupDIA {

void TransitionIdentified::get_product_max_intensity(vector<double>& intensity) const {
	intensity.clear();
	for (int i = 0; i < product_mz.size(); i++) {
		double max = 0;
		for (auto it_exp = product_ions.begin(); it_exp != product_ions.end(); ++it_exp) {
			if (!it_exp->empty()) {
				for (auto jt = it_exp->at(i).intensity.begin(); jt != it_exp->at(i).intensity.end();
						++jt) {
					if (max < *jt)
						max = *jt;
				}
			}
		}
		intensity.push_back(max);
	}
}

void TransitionIdentified::get_product_type(vector<string>& type) const {
//	type = this->product_type;
	type.resize(product_mz.size());
	for (int i = 0; i < product_mz.size(); i++) {
		auto it = map_product_type.find(product_mz.at(i));
		if (it != map_product_type.end())
			type.at(i) = it->second;
		else
			type.at(i) = "-";
	}
}

bool TransitionIdentified::is_exist(int exp_num) const {
	if (exp_num >= precursor_ions.size()) {
		return false;
	}
	if (product_mz.size() <= 0) {
		return false;
	}
	return this->precursor_ions.at(exp_num).is_exist();
}

bool TransitionIdentified::is_exist() const {
	if (product_mz.size() <= 0) {
		return false;
	}
	bool pre_exist = false;
	for (auto it = precursor_ions.begin(); it != precursor_ions.end(); ++it)
		if (it->is_exist())
			pre_exist = true;
	return pre_exist;
}

void TransitionIdentified::get_product_sn(vector<double>& sn, int exp_num) const {
	vector<double> rt;
	this->product_ions.at(exp_num).get_rt(rt);

	sn.resize(get_product_num(), 0);
	for (int i = 0; i < sn.size(); i++) {
		sn.at(i) = ToolsScore::cal_sn_score(product_ions.at(exp_num).at(i).intensity, rt,
				precursor_ions.at(exp_num).get_peak_apex_rt());
	}
}

void TransitionIdentified::get_product_intensity(vector<double>& intensity, int exp_num) {
	if (final_product_intensity.empty())
		calc_product_intensity();

	if (exp_num >= final_product_intensity.size())
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, String(exp_num));

	intensity = final_product_intensity[exp_num];
}

void TransitionIdentified::set_product_by_seq(bool is_selected_by_sequence) {
	vector<double> ori_mz = product_mz;
	vector<double> remained_mz = ori_mz;
	product_mz.clear();
	vector<ProductIons> ori_product = this->product_ions;
	for (auto it = product_ions.begin(); it != product_ions.end(); ++it)
		it->clear();

	TheoreticalSpectrumGenerator spectrum_generator;
	RichPeakSpectrum spectrum;

	spectrum_generator.addPeaks(spectrum, aasequence, Residue::YIon, 1);
	add_product_ions(spectrum, "1_y_", remained_mz);

	spectrum.clear(true);
	spectrum_generator.addPeaks(spectrum, aasequence, Residue::YIon, 2);
	add_product_ions(spectrum, "2_y_", remained_mz);

	spectrum.clear(true);
	spectrum_generator.addPeaks(spectrum, aasequence, Residue::BIon, 1);
	add_product_ions(spectrum, "1_b_", remained_mz);

	spectrum.clear(true);
	spectrum_generator.addPeaks(spectrum, aasequence, Residue::BIon, 2);
	add_product_ions(spectrum, "2_b_", remained_mz);

	if (!is_selected_by_sequence) {
		//Random select product ions
		vector<double> decoy_mz;
		decoy_mz.reserve(remained_mz.size());
		for (auto it = remained_mz.begin(); it != remained_mz.end(); ++it)
			if (*it > 0)
				decoy_mz.push_back(*it);

		//generate uniformed decoy
		int seed = time(0);

		ToolsSpectrum::random(decoy_mz.begin(), decoy_mz.end());
		ToolsSpectrum::random(product_mz.begin(), product_mz.end());

		if (product_mz.size() > decoy_mz.size()) {
			product_mz.resize(decoy_mz.size());
		} else {
			decoy_mz.resize(product_mz.size());
		}

		map<double, double> decoy_map_product_theoretical_mz;
		map<double, string> decoy_map_product_type;
		auto itarget = product_mz.begin(), idecoy = decoy_mz.begin();
		for (; itarget != product_mz.end() and idecoy != decoy_mz.end(); ++itarget, ++idecoy) {
			auto i_theor_mz = map_product_theoretical_mz.find(*itarget);
			decoy_map_product_theoretical_mz[*idecoy] = i_theor_mz->second;

			auto i_type = map_product_type.find(*itarget);
			decoy_map_product_type[*idecoy] = i_type->second;

			*itarget = *idecoy;
		}
		map_product_theoretical_mz = decoy_map_product_theoretical_mz;
		map_product_type = decoy_map_product_type;
	}

	//reassign the product
	sort(product_mz.begin(), product_mz.end());
	for (int i = 0; i < product_mz.size(); i++) {
		//Find the origin place
		int place = 0;
		while (product_mz[i] != ori_mz[place])
			place++;

		//assign the product ion
		for (int cur_exp = 0; cur_exp < get_exp_size(); cur_exp++) {
			if (precursor_ions.size() <= cur_exp or product_ions.size() <= cur_exp)
				break;
			product_ions.at(cur_exp).push_back(ori_product.at(cur_exp)[place]);
		}
	}
}

int TransitionIdentified::add_product_ions(RichPeakSpectrum& spectrum, string add,
		vector<double>& mz) {
	spectrum.sortByPosition();
	int count = 0;
	int num = 0;
	for (auto it = spectrum.begin(); it != spectrum.end(); ++it) {
		num++;
		for (auto jt = mz.begin(); jt != mz.end(); ++jt) {
			if (abs(it->getMZ() - (*jt)) < MS2_WIN) {
				product_mz.push_back(*jt);
				map_product_type.insert(pair<double, string>(*jt, add + String(num)));
				map_product_theoretical_mz.insert(pair<double, double>(*jt, it->getMZ()));
				*jt = 0;
				count++;
				break;
			}
		}
	}
	return count;
}

void TransitionIdentified::get_product_theoretical_mz(vector<double>& mz) const {
	mz.clear();
	for (auto it = product_mz.begin(); it != product_mz.end(); ++it) {
		auto tt = map_product_theoretical_mz.find(*it);
		if (tt != map_product_theoretical_mz.end()) {
			mz.push_back(tt->second);
		}
	}
}

void TransitionIdentified::choose_best_flyer_by_percursor_pattern(vector<int>& remain) {
	remain.clear();

	int exp_num = precursor_ions.size();
	int product_ions_num = product_mz.size();
	if (exp_num <= 0 or product_ions_num < MIN_PRODUCT_IONS_NUM)
		return;
	bool error = false;
	for (auto it = product_ions.begin(); it != product_ions.end(); ++it) {
		if (it->empty())
			error = true;
	}
	if (error) {
		remain.clear();
		return;
	}

// Get the intensity of product ions
	vector<vector<double> > product_intensity(product_ions_num); //Ions, exp
	for (int cur_ions = 0; cur_ions < product_ions_num; cur_ions++) {
		for (int cur_exp = 0; cur_exp < exp_num; cur_exp++) {
			vector<double> cur_intensity;
			int peak_start, peak_length;
			this->precursor_ions[cur_exp].get_peak_rt_in_nearby_no(peak_start, peak_length);
			product_ions[cur_exp].get_intensity(cur_ions, peak_start, peak_length, cur_intensity);
			product_intensity[cur_ions].insert(product_intensity[cur_ions].end(),
					cur_intensity.begin(), cur_intensity.end());
		}
	}

// Get the intensity of precursor ion
	vector<double> precursor_intensity;
	for (int i = 0; i < exp_num; i++) {
		vector<double> ms1;
		this->precursor_ions.at(i).get_peak_intensity(ms1);
		precursor_intensity.insert(precursor_intensity.end(), ms1.begin(), ms1.end());
	}

//Get the intensity for every product ions
	calc_product_intensity();
	vector<double> product_final_intensity(product_ions_num, 0);
	for (int i = 0; i < product_ions_num; i++) {
		for (int j = 0; j < exp_num; j++)
			product_final_intensity[i] += final_product_intensity[j][i];
	}

	////Cal the correlation from each ions to mean
	vector<double> correlation(product_ions_num, 0);
	for (int i = 0; i < product_ions_num; i++) {
		if (product_final_intensity[i] <= 0) {
			correlation[i] = -1;
			continue;
		}
		correlation[i] = OpenMS::Math::pearsonCorrelationCoefficient(precursor_intensity.begin(),
				precursor_intensity.end(), product_intensity[i].begin(),
				product_intensity[i].end());
		if (correlation[i] != correlation[i])
			correlation[i] = -1;
	}

	map<double, int> map_cor_place;
	for (int i = 0; i < correlation.size(); i++) {
		double cor = correlation[i];
		if (cor > 0 and product_final_intensity[i] > 0)
			map_cor_place[cor] = i;
	}
	int select_num = BEST_FLYER_NUM;
	if (select_num > map_cor_place.size())
		select_num = map_cor_place.size();
	auto imap = map_cor_place.rbegin();

	//First find mz >= MIN_FLYER_MZ
	for (int i = 0; i < select_num and imap != map_cor_place.rend(); ++imap) {
		if (product_mz[imap->second] >= MIN_FLYER_MZ) {
			remain.push_back(imap->second);
			i++;
		}
	}

	//If not enough, select mz<MIN_FLYER_MZ
	if (remain.size() < MIN_PRODUCT_IONS_NUM) {
		imap = map_cor_place.rbegin();
		for (int i = remain.size(); i < MIN_PRODUCT_IONS_NUM and imap != map_cor_place.rend();
				++imap) {
			if (product_mz[imap->second] < MIN_FLYER_MZ) {
				remain.push_back(imap->second);
				i++;
			}
		}
	}

	sort(remain.begin(), remain.end());
}

GroupDIA::TransitionIdentified::TransitionIdentified(const Param& para) {
	identification_score = 0;
	accession = "";
	decoy = 0;
	BEST_FLYER_NUM = (int) para.getValue("best_flyer_num");
	MIN_PRODUCT_IONS_NUM_IN_INTENSITY = (int) para.getValue(
			"min_allowed_product_ions_in_calc_intensity");
	MIN_FLYER_MZ = (DoubleReal) para.getValue("wanted_min_product_mz");
	MAX_DIFF_IN_SCORE = (DoubleReal) para.getValue("max_similar_score");
	MIN_PRODUCT_IONS_NUM = (int) para.getValue("min_product_ions_num");
}

GroupDIA::TransitionIdentified::TransitionIdentified() {
	identification_score = 0;
	accession = "";
	decoy = 0;
}

bool GroupDIA::TransitionIdentified::load_store(fstream& f, int io_type, int type, int ion_num) {
	if (type == LoadSave::Type::Basic_Precursor_Iden) {
		LoadSave::io(f, type, io_type);
		if (f.eof())
			return false;
		LoadSave::io(f, this->feature_id, io_type);
		LoadSave::io(f, this->scan_num, io_type);
		LoadSave::io(f, this->feature_exp_num, io_type);
		LoadSave::io(f, this->feature_score, io_type);
		LoadSave::io(f, this->precursor_charge, io_type);
		LoadSave::io(f, this->precursor_mz, io_type);
		LoadSave::io(f, this->decoy, io_type);

		string sequence = "";
		if (io_type == LoadSave::IOType::Write)
			sequence = aasequence.toString();
		LoadSave::io(f, sequence, io_type);
		if (!sequence.empty() and io_type == LoadSave::IOType::Read)
			this->set_pep_seq(sequence);

		LoadSave::io(f, this->accession, io_type);
		LoadSave::io(f, this->identification_score, io_type);
		LoadSave::vector_io(f, this->product_mz, io_type);
		LoadSave::map_io(f, map_product_theoretical_mz, io_type);
		LoadSave::map_io(f, map_product_type, io_type);
		LoadSave::map_io(f, map_product_intensity, io_type);
		LoadSave::ion_io(f, feature_precursor, io_type);
		LoadSave::vector_ion_io(f, precursor_ions, io_type);
		return true;

	} else {
		return TransitionGroup::load_store(f, io_type, type, ion_num);
	}
}

void TransitionIdentified::reselect_product_ions(const Param& para) {
//Select best flyer by calc the correlation between precursor and product ions.
	vector<int> remain;
	if (!this->is_decoy()) {
		if (para.getValue("library_method") == "correlation") {
			//For target, select best correlation product ions.
			choose_best_flyer_by_percursor_pattern(remain);
		} else if (para.getValue("library_method") == "max") {
			choose_best_flyer_by_select_max_intensity(remain);
		} else {
			throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
					"library_method error:", para.getValue("library_method"));
		}
	} else
		//For decoy, random select product ions
		choose_best_flyer_by_random_select(remain);

//Remove not needed ions
	this->refine_product_ions(remain);

//Calculate product intensity
	calc_product_intensity();

//Generate library intensity
	double precursor_intensity;
	vector<double> product_intensity;
	bool is_exist = get_library_intensity(precursor_intensity, product_intensity);
	if (is_exist) {
		for (int i = 0; i < product_mz.size(); i++) {
			double mz = product_mz[i];
			map_product_intensity[mz] = product_intensity[i];
		}
	} else {
		this->product_mz.clear();
		for (auto it = this->product_ions.begin(); it != this->product_ions.end(); ++it)
			it->clear();
	}
}

string TransitionIdentified::get_pep_seq() const {
	return this->aasequence.toString();
}

string TransitionIdentified::get_accession() const {
	return this->accession;
}

void TransitionIdentified::generate_decoy_by_simple_copy(TransitionIdentified& decoy_trans) const {
//generate decoy transition
	decoy_trans = *this;
	decoy_trans.set_decoy(decoy_trans.Decoy_Rand);
	decoy_trans.accession = "DECOY_" + decoy_trans.accession;
}

void TransitionIdentified::choose_best_flyer_by_random_select(vector<int>& remain) const {
	remain.clear();

	vector<int> ori(product_mz.size());
	for (int i = 0; i < ori.size(); i++)
		ori[i] = i;
	ToolsSpectrum::random(ori.begin(), ori.end());
	int select_num = BEST_FLYER_NUM;
	if (select_num > ori.size())
		select_num = ori.size();
	for (int i = 0; i < select_num; i++) {
		remain.push_back(ori[i]);
	}
	sort(remain.begin(), remain.end());
}

void TransitionIdentified::choose_best_flyer_by_select_max_intensity(vector<int>& remain) {
	remain.clear();
// Calculate intensity for all ions
	calc_product_intensity();
	int product_ions_num = product_mz.size();
	int exp_num = precursor_ions.size();
	vector<double> product_intensity(product_ions_num, 0);
	for (int i = 0; i < product_ions_num; i++) {
		for (int j = 0; j < exp_num; j++)
			product_intensity[i] += final_product_intensity[j][i];
	}

	vector<pair<double, int> > intensity_place;
	for (int i = 0; i < product_ions_num; i++) {
		intensity_place.push_back(pair<double, int>(product_intensity[i], i));
	}

	sort(intensity_place.begin(), intensity_place.end(), large_first<double, int>);

	int select_num = BEST_FLYER_NUM;
	if (select_num > intensity_place.size())
		select_num = intensity_place.size();

	for (int i = 0; i < select_num; i++) {
		remain.push_back(intensity_place[i].second);
	}

	sort(remain.begin(), remain.end());
}

void TransitionIdentified::set_identification_result(const PeptideHit& pephit, bool is_decoy) {
	if (pephit.getProteinAccessions().size() < 1)
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Protein accession not found, Sequence is " + pephit.getSequence().toString());
	this->accession = pephit.getProteinAccessions()[0];
	for (auto it = pephit.getProteinAccessions().begin() + 1;
			it != pephit.getProteinAccessions().end(); ++it)
		this->accession = this->accession + ";" + *it;
	this->aasequence = pephit.getSequence();
	this->identification_score = pephit.getScore();
}

double TransitionIdentified::get_identification_score() const {
	return this->identification_score;
}

void TransitionIdentified::calc_product_intensity() {
	final_product_intensity.clear();
	final_product_intensity.resize(this->get_exp_size());

//If product ions nums are less than min required num, the intensity is 0.
	if (this->get_product_num() < MIN_PRODUCT_IONS_NUM_IN_INTENSITY or product_ions.empty()) {
		for (int exp_num = 0; exp_num < this->get_exp_size(); exp_num++) {
			final_product_intensity[exp_num] = vector<double>(this->get_product_num(), 0);
		}
		return;
	}
	if (product_ions[0].size() != product_mz.size())
		cout << "Error in TransitionIdentified::calc_product_intensity, product ions size not equal"
				<< endl;

	for (int exp_num = 0; exp_num < this->get_exp_size(); exp_num++) {
		final_product_intensity[exp_num] = vector<double>();
		//Get peak info
		int peak_start, peak_length;
		this->precursor_ions.at(exp_num).get_peak_rt_in_nearby_no(peak_start, peak_length);
		int peak_apex_no = this->precursor_ions.at(exp_num).get_peak_apex_rt_in_nearby_no();

		//Get production ion intensity
		vector<vector<double> > origin_intensity(this->get_product_num());
		for (int i = 0; i < this->get_product_num(); i++) {
			origin_intensity.at(i) = this->product_ions.at(exp_num).at(i).get_intensity();
		}

		//calc the intensity
		IntensityCalc int_calc;
		int_calc.get_intensity(origin_intensity, peak_start, peak_length,
				final_product_intensity[exp_num]);

		int not_zero_num = 0;
		for (auto it = final_product_intensity[exp_num].begin();
				it != final_product_intensity[exp_num].end(); ++it) {
			if (*it > 0)
				not_zero_num++;
		}
		if (not_zero_num < MIN_PRODUCT_IONS_NUM_IN_INTENSITY) {
			//Is empty
			for (auto& it : final_product_intensity[exp_num])
				it = 0;
//			precursor_ions[exp_num].not_exist();
		}
	}
}

double TransitionIdentified::calc_dist_score(const vector<double>& _ms1, vector<double>& ms2,
		double lib_ms1, double lib_ms2) const {
	if (_ms1.size() != ms2.size())
		cout << "Error in calc_dist_score\n";
//Normalize ms1 and ms2
	auto ms1 = _ms1;
	double sum_ms1 = (accumulate(ms1.begin(), ms1.end(), 0.0)) / lib_ms1;
	double sum_ms2 = (accumulate(ms2.begin(), ms2.end(), 0.0)) / lib_ms2;

	if (sum_ms1 == 0 or sum_ms2 == 0)
		return MAX_DIFF_IN_SCORE;

	double score = 0;
	if (sum_ms1 > sum_ms2) {
		auto i1 = ms1.begin();
		auto i2 = ms2.begin();
		auto end1 = ms1.end();
		for (; i1 != end1; ++i1, ++i2) {
			double delta = *i1 / *i2 / lib_ms1 * lib_ms2;
			if (delta > MAX_DIFF_IN_SCORE or delta != delta)
				delta = MAX_DIFF_IN_SCORE;
			score += delta * delta;
		}
	} else {
		auto i1 = ms1.begin();
		auto i2 = ms2.begin();
		auto end1 = ms1.end();
		for (; i1 != end1; ++i1, ++i2) {
			double delta = *i2 / *i1 / lib_ms2 * lib_ms1;
			if (delta > MAX_DIFF_IN_SCORE or delta != delta)
				delta = MAX_DIFF_IN_SCORE;
			score += delta * delta;
		}
	}

	return sqrt(score) / ms1.size();
}

void TransitionIdentified::generate_mrm(int exp_num, const Param& param_,
		MRMTransitionGroup<MSSpectrum<ChromatogramPeak>, OpenSwath::LightTransition>& mrm) const {
	OpenSwath::LightTargetedExperiment exp;

	if (product_mz.size() == 0)
		cout << "The size of product mz is zero.";

//Build transition
	for (int i = 0; i < product_mz.size(); i++) {
		double mz = product_mz[i];
		OpenSwath::LightTransition t;
		t.product_mz = mz;
		t.precursor_mz = precursor_mz;
		t.peptide_ref = this->get_pep_seq();

		auto iint = map_product_intensity.find(mz);
		if (iint != map_product_intensity.end())
			t.library_intensity = iint->second;
		else
			cout << "Not find product intensity" << this->get_pep_seq() << endl;

		t.transition_name = String(mz);
//		t.charge = 1;
		auto itype = map_product_type.find(mz);
		if (itype != map_product_type.end())
			t.transition_name = itype->second;
		else
			cout << "Not find product type" << this->get_pep_seq() << endl;
		t.charge = String(t.transition_name.substr(0, 1)).toInt();
		mrm.addTransition(t, t.getNativeID());

		//Build chromatrom
		auto& product = product_ions[exp_num];
		vector<double> intensity;
		product.get_intensity(i, intensity);
		MSSpectrum<ChromatogramPeak> spec;
		spec.reserve(intensity.size());
		vector<double> rt;
		product_ions[exp_num].get_rt(rt);
		if (rt.size() != intensity.size())
			cout << "RT size not equal to intensity size." << rt.size() << "\t" << intensity.size()
					<< endl;

		auto irt = rt.begin();
		for (auto it = intensity.begin(); it != intensity.end(); ++it, ++irt) {
			ChromatogramPeak peak;
			peak.setIntensity(*it);
			peak.setMZ(*irt);
			spec.push_back(peak);
		}

		spec.setMetaValue("product_mz", product_mz[i]);
		spec.setMetaValue("precursor_mz", this->precursor_mz);
		spec.setNativeID(t.getNativeID());
		mrm.addChromatogram(spec, t.getNativeID());
	}
	mrm.setTransitionGroupID(this->get_scan_num());

//	cout << "Start pick\n";
//Pick mrm
	MRMTransitionGroupPicker trgroup_picker;
	trgroup_picker.setParameters(param_);
	trgroup_picker.pickTransitionGroup(mrm);
}

void TransitionIdentified::get_pep_aaseq(AASequence& seq) const {
	seq = aasequence;
}

bool TransitionIdentified::is_decoy() const {
	if (decoy == this->Target)
		return false;
	else
		return true;
}
int TransitionIdentified::get_decoy_type() const {
	return this->decoy;
}
void TransitionIdentified::set_decoy(int decoy_type) {
	decoy = decoy_type;
}

void TransitionIdentified::get_information(map<string, string>& map_info) const {
//Add decoy information
//	int decoy = this->get_decoy_type();
//	if (decoy == this->Target) {
//		map_info.insert(pair<string, string>("decoy_all", String(0)));
//		map_info.insert(pair<string, string>("decoy_spec", String(0)));
//		map_info.insert(pair<string, string>("decoy_iden", String(0)));
//		map_info.insert(pair<string, string>("decoy", String(0)));
//	} else if (decoy == this->Decoy_Iden) {
//		map_info.insert(pair<string, string>("decoy_all", String(1)));
//		map_info.insert(pair<string, string>("decoy_spec", String(0)));
//		map_info.insert(pair<string, string>("decoy_iden", String(1)));
//		map_info.insert(pair<string, string>("decoy", String(0)));
//	} else if (decoy == this->Decoy_Rand) {
//		map_info.insert(pair<string, string>("decoy_all", String(1)));
//		map_info.insert(pair<string, string>("decoy_spec", String(0)));
//		map_info.insert(pair<string, string>("decoy_iden", String(0)));
//		map_info.insert(pair<string, string>("decoy", String(1)));
//	} else if (decoy == this->Decoy_Spec) {
//		map_info.insert(pair<string, string>("decoy_all", String(1)));
//		map_info.insert(pair<string, string>("decoy_spec", String(1)));
//		map_info.insert(pair<string, string>("decoy_iden", String(0)));
//		map_info.insert(pair<string, string>("decoy", String(0)));
//	}

	map_info.insert(pair<string, string>("ide_score", String(this->get_identification_score())));

//Add identification information
	map_info.insert(pair<string, string>("accessions", this->get_accession()));

	AASequence aaseq;
	this->get_pep_aaseq(aaseq);
	map_info.insert(pair<string, string>("full_sequence", aaseq.toString()));
	map_info.insert(pair<string, string>("sequence", aaseq.toUnmodifiedString()));

//Add other information
	map_info.insert(pair<string, string>("charge", String(this->get_precursor_charge())));
	map_info.insert(pair<string, string>("mz", String(this->get_precursor_mz())));
	int scan_num = this->get_scan_num();
	map_info.insert(pair<string, string>("scan_num", String(scan_num)));
	map_info.insert(pair<string, string>("id", String(this->feature_id)));
}

bool TransitionIdentified::get_library_intensity(double& lib_ms1_intensity,
		vector<double>& lib_ms2_intensity) {
	int product_size = this->get_product_num();
	vector<double> all_ms1_intensity;
	this->get_precursor_peak_intensity(all_ms1_intensity);
	lib_ms1_intensity = accumulate(all_ms1_intensity.begin(), all_ms1_intensity.end(), 0.0);
	if (lib_ms1_intensity == 0) {
		return false;
	}

//add ms2 intensity
	lib_ms2_intensity = vector<double>(product_size, 0);

#ifdef SINGLE
	int cur_exp = this->feature_exp_num;
	for (int i = 0; i < product_size; i++) {
		lib_ms2_intensity[i] += this->final_product_intensity[cur_exp][i];
	}
	lib_ms1_intensity = all_ms1_intensity[cur_exp];
	if ((accumulate(lib_ms2_intensity.begin(), lib_ms2_intensity.end(), 0.0) == 0)
			or (lib_ms1_intensity == 0))
	return false;
	else
	return true;
#endif

	for (int cur_exp = 0; cur_exp < this->get_exp_size(); cur_exp++) {
//		if (!this->precursor_ions[cur_exp].is_exist())
//			continue;

		for (int i = 0; i < product_size; i++) {
			lib_ms2_intensity[i] += this->final_product_intensity[cur_exp][i];
		}
	}
	if (accumulate(lib_ms2_intensity.begin(), lib_ms2_intensity.end(), 0.0) == 0)
		return false;
	return true;
}

void TransitionIdentified::reset_peak(int exp_num, double left_peak_rt, double right_peak_rt) {
	precursor_ions[exp_num].set_peak(left_peak_rt, right_peak_rt);
}

int TransitionIdentified::get_exist_exp_num() const {
	int exist_exp_num = 0;
	for (int i = 0; i < precursor_ions.size(); i++)
		if (precursor_ions[i].is_exist())
			exist_exp_num++;
	return exist_exp_num;
}

void TransitionIdentified::generate_decoy_by_openswath_method(TransitionIdentified& decoy_trans,
		bool& is_right, String method, double identity_threshold, int max_attempts,
		double mz_threshold, double mz_shift, double similarity_threshold) const {
	decoy_trans = *this;
	decoy_trans.set_decoy(decoy_trans.Decoy_Rand);
	decoy_trans.accession = "DECOY_" + decoy_trans.accession;

	is_right = true;

	MRMDecoy mrm;

// Generate decoy sequence
	OpenMS::String original_sequence = this->aasequence.toUnmodifiedString();
	OpenMS::TargetedExperiment::Peptide target_peptide, decoy_peptide;

#pragma omp critical
	ToolsSpectrum::convert_aasequence_to_peptide(this->aasequence, target_peptide);
	target_peptide.setChargeState(this->get_precursor_charge());

	if (method == "pseudo-reverse") {
		decoy_peptide = mrm.pseudoreversePeptide(target_peptide);
	} else if (method == "reverse") {
		decoy_peptide = mrm.reversePeptide(target_peptide);
	} else if (method == "shuffle") {
		decoy_peptide = mrm.shufflePeptide(target_peptide, identity_threshold, -1, max_attempts);
	}

	if (mrm.AASequenceIdentity(original_sequence, decoy_peptide.sequence) > identity_threshold) {
		is_right = false;
		return;
	}

	decoy_trans.aasequence = TargetedExperimentHelper::getAASequence(decoy_peptide);
	MRMDecoy::IonSeries decoy_ionseries = mrm.getIonSeries(decoy_trans.aasequence,
			decoy_peptide.getChargeState());
	MRMDecoy::IonSeries target_ionseries = mrm.getIonSeries(this->aasequence,
			target_peptide.getChargeState());

	for (int i = 0; i < this->product_mz.size(); i++) {
		std::pair<String, double> targetion = mrm.getTargetIon(this->product_mz[i], mz_threshold,
				target_ionseries);
		std::pair<String, double> decoyion = mrm.getDecoyIon(targetion.first, decoy_ionseries);
		if (targetion.second == -1 or decoyion.second == -1) {
			decoy_trans.product_mz[i] = -1;
			continue;
		}

		if (method == "shift") {
			decoy_trans.product_mz[i] = (decoyion.second + mz_shift);
		} else {
			decoy_trans.product_mz[i] = (decoyion.second);
		}

		if (decoyion.second > 0) {
			if (similarity_threshold >= 0) {
				if (std::fabs(this->product_mz[i] - decoy_trans.product_mz[i])
						< similarity_threshold) {
					decoy_trans.product_mz[i] = -1;
				}
			}
		}
	}

	vector<double> new_mz;
	for (int i = 0; i < decoy_trans.product_mz.size(); i++) {
		if (decoy_trans.product_mz[i] > 0) {
			new_mz.push_back(decoy_trans.product_mz[i]);
		}
	}

	if (new_mz.size() < MIN_PRODUCT_IONS_NUM) {
		is_right = false;
		return;
	} else {
		decoy_trans.product_mz = new_mz;
	}
}

void TransitionIdentified::add_product_mz(const vector<double>& product_mz) {
	this->product_mz.insert(this->product_mz.end(), product_mz.begin(), product_mz.end());
}
} /* namespace GroupDIA */
