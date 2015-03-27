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
#include "TransitionGroup.h"

//double MS2_WIN = 0.05;

namespace GroupDIA {
void TransitionGroup::set_pep_seq(const string& seq) {
#pragma omp critical
	aasequence.setStringSequence(seq);
//	this->sequence = seq;
}
void TransitionGroup::set_pep_seq(const AASequence& seq) {
	aasequence = seq;
}

bool TransitionGroup::set_transition(const SwathExperiment& exp, const Feature& f, int exp_num,
		int scan_num, int exact_scan_num) {
	this->scan_num = scan_num;
	this->feature_id = f.getUniqueId();
	this->feature_score = f.getOverallQuality();
	this->precursor_charge = f.getCharge();
	this->precursor_mz = f.getMZ();
	this->feature_exp_num = exp_num;
	bool is_exist = this->feature_precursor.set_ms1(exp, precursor_mz, f.getRT(), exact_scan_num);
	return is_exist;
}

void TransitionGroup::add_ms1(const SwathExperiment& exp, const RTNormalizer& rt, int exp_num,
		int exact_scan_num) {
	if (exp_num >= precursor_ions.size()) {
		precursor_ions.resize(exp_num + 1);
	}
	if (feature_exp_num == exp_num)
		this->precursor_ions[exp_num] = feature_precursor;
	else
		feature_precursor.add_ms1(exp, precursor_mz, rt, exact_scan_num,
				this->precursor_ions[exp_num]);
}

void TransitionGroup::set_product_ions_mz(const SwathExperiment& exp, double ms1_win_start,
		double ms1_win_end) {
	const auto& ori_spec = exp.get_spec(feature_precursor.get_peak_apex_rt_no());

	//Find peak
	MSSpectrum<Peak1D> filted_spectrum;
	ToolsSwathGaussFilter spec;
	spec.swath_filter(ori_spec, filted_spectrum);
	filted_spectrum.sortByPosition();

	MSSpectrum<Peak1D> picked_spectrum;
	ToolsSwathPickPeak pp;
	pp.pick_peak(filted_spectrum, picked_spectrum);
	picked_spectrum.sortByPosition();

//	MSSpectrum<Peak1D> refined_spectrum;
//	pp.refine_peak(picked_spectrum, exp.get_para(), refined_spectrum);

	MSSpectrum<Peak1D> result_spectrum;
	ToolsDeIsotoper deiso;
	deiso.de_isotoper(picked_spectrum, result_spectrum);

	product_mz.clear();
	for (auto it = result_spectrum.begin(); it != result_spectrum.end(); ++it) {
		if (it->getMZ() >= ms1_win_start and it->getMZ() <= ms1_win_end)
			continue;
		if (it->getIntensity() > 0) {
			product_mz.push_back(it->getMZ());
		}
	}

//#pragma omp critical (debug)
//	if (feature_precursor.get_peak_apex_rt_no() == 712) {
//		for (auto it = ori_spec.begin(); it != ori_spec.end(); ++it) {
//			cout << it->getMZ() << "\t";
//		}
//		cout << endl;
//		cout << feature_precursor.get_peak_apex_rt_no() << endl;
//		for (auto it = product_mz.begin(); it != product_mz.end(); ++it) {
//			cout << *it << "\t";
//		}
//		cout << endl;
//	}
}

void TransitionGroup::get_precursor_nearby_intensity(int exp_num, vector<double>& intensity) const {
	precursor_ions[exp_num].get_nearby_intensity(intensity);
}

void GroupDIA::TransitionGroup::set_scan_num(int num) {
	this->scan_num = num;
}

void TransitionGroup::set_precursor_rt_no(int exp_num, int& peak_start_no, int& peak_length) {
	precursor_ions[exp_num].set_peak_rt_no(peak_start_no, peak_length);
}

bool TransitionGroup::load_store(fstream& f, int io_type, int type, int ion_num) {
	LoadSave::io(f, type, io_type);
	if (f.eof())
		return false;

	if (type == LoadSave::Type::Basic) {
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
		LoadSave::ion_io(f, feature_precursor, io_type);
		return true;

	} else if (type == LoadSave::Type::Feature_Product) {
		int id = scan_num;
		LoadSave::io(f, id, io_type);
		this->feature_product.load_store(f, io_type);
		if (id != scan_num) {
			this->feature_product = ProductIons();
			return false;
		} else
			return true;

	} else if (type == LoadSave::Type::Basic_Precursor) {
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
		LoadSave::ion_io(f, feature_precursor, io_type);
		LoadSave::vector_ion_io(f, precursor_ions, io_type);
		return true;

	} else if (type == LoadSave::Type::Precursor) {
		int id = scan_num;
		LoadSave::io(f, id, io_type);
		if (ion_num >= this->precursor_ions.size())
			precursor_ions.resize(ion_num + 1);
		this->precursor_ions[ion_num].load_store(f, io_type);
		if (id != scan_num) {
			this->precursor_ions[ion_num] = PrecursorIon();
			return false;
		} else
			return true;

	} else if (type == LoadSave::Type::Product) {
		int id = scan_num;
		LoadSave::io(f, id, io_type);
		if (ion_num >= product_ions.size())
			product_ions.resize(ion_num + 1);
		this->product_ions[ion_num].load_store(f, io_type);
		if (id != scan_num) {
			this->product_ions[ion_num] = ProductIons();
			return false;
		} else
			return true;

	} else if (type == LoadSave::Type::All) {
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
		LoadSave::ion_io(f, feature_precursor, io_type);

		LoadSave::vector_ion_io(f, precursor_ions, io_type);
		LoadSave::vector_ion_io(f, product_ions, io_type);
		return true;
	} else {
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Error in load: " + String(type));
	}
	return false;
}

void TransitionGroup::get_product_mz(vector<double>& mz) const {
	mz = this->product_mz;
}

UInt64 TransitionGroup::get_uid() const {
	return this->feature_id;
}

int TransitionGroup::get_exp_size() const {
	return precursor_ions.size();
}

double TransitionGroup::get_precursor_mz() const {
	return this->precursor_mz;
}

int TransitionGroup::get_feature_exp_num() const {
	return this->feature_exp_num;
}

int TransitionGroup::get_product_num() const {
	return this->product_mz.size();
}

int TransitionGroup::get_precursor_charge() const {
	return this->precursor_charge;
}

double TransitionGroup::get_precursor_rt(int exp_num) const {
	return precursor_ions[exp_num].get_peak_apex_rt();
}

void TransitionGroup::refine_product_ions(const vector<int>& remain_place) {
	vector<double> remain_mz(remain_place.size());
	for (int i = 0; i < remain_place.size(); i++) {
		remain_mz.at(i) = product_mz.at(remain_place.at(i));
	}
	product_mz = remain_mz;
	for (int i = 0; i < this->product_ions.size(); i++) {
		this->product_ions[i].refine_product_ions(remain_place);
	}
}

int TransitionGroup::size() const {
	return precursor_ions.size();
}

int TransitionGroup::get_scan_num() const {
	return scan_num;
}

void TransitionGroup::record_ms2_nearby(SwathExperiment& exp, int exp_num) {
	int rt_start_no, rt_length;
	if (exp_num >= precursor_ions.size()) {
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"record_ms2_nearby error");
	}
	precursor_ions[exp_num].get_nearby_rt_no(rt_start_no, rt_length);
	if (rt_start_no + rt_length > exp.size()) {
		cout << "Error in record_ms2, ID: " << this->get_uid() << ", rt_start_no: " << rt_start_no
				<< ", rt_length: " << rt_length << "." << endl;
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"record_ms2_nearby error");
	}
	exp.record_mz(this->get_scan_num(), rt_start_no, rt_length, product_mz);
}

void TransitionGroup::record_feature_product_peak(SwathExperiment& exp) {
	int peak_rt_start_no, peak_rt_length;
	feature_precursor.get_peak_rt_no(peak_rt_start_no, peak_rt_length);
	exp.record_mz(this->scan_num, peak_rt_start_no, peak_rt_length, product_mz);
}

void TransitionGroup::add_feature_product(const SwathExperiment& exp) {
//	int peak_rt_start_no, peak_rt_length;
//	feature_precursor.get_peak_rt_no(peak_rt_start_no, peak_rt_length);
//	feature_product.add_ms2(exp, product_mz, peak_rt_start_no, peak_rt_length);

	int nearby_rt_start_no, nearby_rt_length;
	feature_precursor.get_nearby_rt_no(nearby_rt_start_no, nearby_rt_length);
	feature_product.add_ms2(exp, product_mz, nearby_rt_start_no, nearby_rt_length);
}

void TransitionGroup::add_ms2_nearby_in_record_mode(const SwathExperiment& exp, int exp_num) {
	if (exp_num >= product_ions.size()) {
		product_ions.resize(exp_num + 1);
	}

	int nearby_rt_start_no, nearby_rt_length;
	precursor_ions[exp_num].get_nearby_rt_no(nearby_rt_start_no, nearby_rt_length);
	product_ions[exp_num].add_ms2(exp, this->get_scan_num(), product_mz, nearby_rt_start_no,
			nearby_rt_length);

	if (nearby_rt_start_no + nearby_rt_length > exp.size())
		cout << "Error in record_ms2, ID: " << this->get_uid() << ", rt_start_no: "
				<< nearby_rt_start_no << ", rt_length: " << nearby_rt_length << "." << endl;
}

void TransitionGroup::add_ms2_nearby(const SwathExperiment& exp, int exp_num) {
//	if (product_mz.size() != feature_product.size())
//		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
//				"product_mz not equal to feature_product size, product_mz size:"
//						+ String(product_mz.size()) + ", feature_product size:"
//						+ String(feature_product.size()));

	if (exp_num >= precursor_ions.size()) {
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Error in add_ms2_nearby, size of precursor_ions is:"
						+ String(precursor_ions.size()) + ", exp_num is:" + String(exp_num));
	}

	if (exp_num >= product_ions.size()) {
		product_ions.resize(exp_num + 1);
	}

	int nearby_rt_start_no, nearby_rt_length;
	precursor_ions[exp_num].get_nearby_rt_no(nearby_rt_start_no, nearby_rt_length);
	product_ions[exp_num].add_ms2(exp, product_mz, nearby_rt_start_no, nearby_rt_length);

	if (nearby_rt_start_no + nearby_rt_length > exp.size())
		cout << "Error in record_ms2, ID: " << this->get_uid() << ", rt_start_no: "
				<< nearby_rt_start_no << ", rt_length: " << nearby_rt_length << "." << endl;
}

void TransitionGroup::add_ms1_nearby(const SwathExperiment& exp, int exp_num, int exact_scan_num) {
	if (exp_num >= precursor_ions.size()) {
		cout << "Error in TransitionGroup::add_ms1_nearby" << endl;
	}

	int rt_start, rt_length;
	precursor_ions[exp_num].get_peak_rt_no(rt_start, rt_length);
	this->precursor_ions[exp_num].set_peak(exp, precursor_mz, rt_start, rt_length, exact_scan_num);
}

void TransitionGroup::get_product_intensity_in_each_time(int cur_exp_num,
		vector<vector<double> >& intensity) const {
	intensity.resize(product_mz.size());
	auto iadd = intensity.begin();
	for (auto it = product_ions[cur_exp_num].begin(); it != product_ions[cur_exp_num].end();
			++iadd, ++it) {
		*iadd = it->get_intensity();
	}
}

void TransitionGroup::get_precursor_peak_intensity(vector<double>& intensity) const {
	intensity.clear();
	for (auto it = precursor_ions.begin(); it != precursor_ions.end(); ++it) {
		if (it->is_exist()) {
			vector<double> cur_int;
			it->get_peak_intensity(cur_int);
			intensity.insert(intensity.end(), cur_int.begin(), cur_int.end());
		} else
			intensity.push_back(0);
	}
}

void TransitionGroup::get_product_intensity_in_each_time_nearby_version(
		vector<vector<double> >& intensity) const {
	intensity.clear();
	intensity.resize(product_mz.size());
	for (int i = 0; i < get_exp_size(); i++) {
		int peak_start, peak_length;
		precursor_ions[i].get_peak_rt_in_nearby_no(peak_start, peak_length);
		for (int mz_no = 0; mz_no < product_mz.size(); mz_no++) {
			vector<double> cur_int;
			product_ions[i].get_intensity(mz_no, cur_int);
			if (peak_length > 0)
				intensity[mz_no].insert(intensity[mz_no].end(), cur_int.begin() + peak_start,
						cur_int.begin() + peak_start + peak_length);
			else
				intensity[mz_no].push_back(0);
		}
	}
}

int TransitionGroup::get_product_rt_no(int exp_num, double rt) const {
	return product_ions[exp_num].get_rt_no(rt);
}

void TransitionGroup::set_uid(UInt64 uid) {
	feature_id = uid;
}

int TransitionGroup::get_product_rt_size(int exp_num, double rt_start, double rt_end) const {
	return product_ions[exp_num].get_rt_size(rt_start, rt_end);
}

void TransitionGroup::remove_nearby_intensity(int exp_num) {
	int peak_start_in_nearby_no, peak_length;
	precursor_ions[exp_num].get_peak_rt_in_nearby_no(peak_start_in_nearby_no, peak_length);
	precursor_ions[exp_num].remove_nearby_intensity();
	product_ions[exp_num].remove_nearby_intensity(peak_start_in_nearby_no, peak_length);
}

void TransitionGroup::remove_feature_product() {
	ProductIons().swap(feature_product);
}

void TransitionGroup::remove_ions() {
	vector<PrecursorIon>().swap(precursor_ions);
	vector<ProductIons>().swap(product_ions);
}

void TransitionGroup::crude_refine(double MIN_ALLOWED_CORRELATION) {
//Get precursor and product ions intensity
	vector<double> precursor_intensity;
	vector<vector<double> > product_intensity;

	feature_precursor.get_peak_intensity(precursor_intensity);
	int peak_rt_from, peak_length;
	feature_precursor.get_peak_rt_in_nearby_no(peak_rt_from, peak_length);
	feature_product.get_intensity(peak_rt_from, peak_length, product_intensity);

//Refine
	int mz_size = this->product_mz.size();
	vector<double> new_product_mz;
	vector<int> remain;
	for (int i = 0; i < mz_size; i++) {
		double cor = Math::pearsonCorrelationCoefficient(precursor_intensity.begin(),
				precursor_intensity.end(), product_intensity.at(i).begin(),
				product_intensity.at(i).end());
		if (cor > MIN_ALLOWED_CORRELATION) {
			new_product_mz.push_back(product_mz.at(i));
			remain.push_back(i);
		}
	}
	product_mz = new_product_mz;
	product_mz.shrink_to_fit();
	feature_product.refine_product_ions(remain);

#pragma omp critical
	if (product_mz.size() != feature_product.size())
		cout << product_mz.size() << "\t" << feature_product.size() << endl;

}

void TransitionGroup::refined_by_crosscorrelation(int cur_exp_num, int MAX_ALLOWED_DELAY,
		int FINE_NEARBY, double MIN_INTENSITY) {
	int exp_size = this->get_exp_size();
	int mz_size = this->product_mz.size();

//Get feature precursor and product intensity
	vector<double> ori_feature_precursor_intensity; //intensity
	vector<vector<double> > feature_product_intensity; //mz, intensity
	int feature_peak_start_in_nearby_no, feature_peak_length;
	feature_precursor.get_peak_rt_in_nearby_no(feature_peak_start_in_nearby_no,
			feature_peak_length);
	feature_precursor.get_nearby_intensity(ori_feature_precursor_intensity);
	feature_product.get_intensity(feature_product_intensity);

//Get mz and intensity
	vector<pair<double, double> > all_possible_mz_intensity;
	all_possible_mz_intensity.reserve(mz_size);
	for (int i = 0; i < mz_size; i++) {
		const double& mz = product_mz.at(i);
		all_possible_mz_intensity.push_back(
				pair<double, double>(mz,
						accumulate(
								feature_product_intensity.at(i).begin()
										+ feature_peak_start_in_nearby_no,
								feature_product_intensity.at(i).begin()
										+ feature_peak_start_in_nearby_no + feature_peak_length,
								0.0)));
	}

	if (all_possible_mz_intensity.size() != mz_size)
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Error in TransitionSimilar: " + this->feature_id);
//Refine peak position
	vector<double> feature_precursor_intensity = ori_feature_precursor_intensity;

	vector<double> cur_precursor_intensity;
	vector<vector<double> > cur_product_intensity;

	precursor_ions.at(cur_exp_num).get_nearby_intensity(cur_precursor_intensity);
	product_ions.at(cur_exp_num).get_intensity(cur_product_intensity);

	double corrected_offset = 0;
	double all_intensity = 0;

	vector<int> mz_used_for_calc_delay;
	for (int cur_mz_num = 0; cur_mz_num < all_possible_mz_intensity.size(); cur_mz_num++) {
		if (all_possible_mz_intensity.at(cur_mz_num).second > MIN_INTENSITY) {
			mz_used_for_calc_delay.push_back(cur_mz_num);
			all_intensity += all_possible_mz_intensity.at(cur_mz_num).second;
		}
	}

	if (all_intensity > 0) {
////////////////////////// Calc the delay first /////////////////////////////////
		map<int, double> delay_score;
		for (int cur_mz_num : mz_used_for_calc_delay) {
			int min_size = cur_product_intensity.at(cur_mz_num).size();
			if (min_size > feature_product_intensity.at(cur_mz_num).size())
				min_size = feature_product_intensity.at(cur_mz_num).size();
			vector<double> feature_product_intensity_cur_mz(
					feature_product_intensity.at(cur_mz_num).begin(),
					feature_product_intensity.at(cur_mz_num).begin() + min_size);
			vector<double> cur_product_intensity_cur_mz(
					cur_product_intensity.at(cur_mz_num).begin(),
					cur_product_intensity.at(cur_mz_num).begin() + min_size);
			//offset=b-a
			vector<pair<int, double> > cur_delay_score;
			ToolsSpectrum::get_best_delay(feature_product_intensity_cur_mz,
					cur_product_intensity_cur_mz, MAX_ALLOWED_DELAY, cur_delay_score);

			for (auto& t : cur_delay_score) {
				auto id = delay_score.find(t.first);
				if (id == delay_score.end()) {
					delay_score.insert(
							pair<int, double>(t.first,
									t.second * all_possible_mz_intensity.at(cur_mz_num).second));
				} else {
					id->second += t.second * all_possible_mz_intensity.at(cur_mz_num).second;
				}
			}
		}
		//Find max delay
		int max_delay = 0;
		double max_score = 0;
		for (auto& t : delay_score) {
			if (t.second > max_score) {
				max_delay = t.first;
				max_score = t.second;
			}
		}
		corrected_offset = max_delay;

////////////////////////// Calc the delay Again /////////////////////////////////
		delay_score.clear();
		//Calc feature start and length
		int feature_start = feature_peak_start_in_nearby_no - FINE_NEARBY;
		int feature_length = feature_peak_length + 2 * FINE_NEARBY;
		if (feature_start < 0)
			feature_start = 0;
		if (feature_start >= feature_precursor_intensity.size() - 1)
			feature_start = feature_precursor_intensity.size();
		if (feature_start + feature_length > feature_precursor_intensity.size())
			feature_length = feature_precursor_intensity.size() - feature_start;

		//Calc cur start and length
		int cur_start = feature_start;
		int cur_lenth = feature_length;
		cur_start += corrected_offset;
		if (cur_start < 0)
			cur_start = 0;
		if (cur_start >= cur_precursor_intensity.size())
			cur_start = cur_precursor_intensity.size() - 1;
		if (cur_start + cur_lenth > cur_precursor_intensity.size()) {
			cur_lenth = cur_precursor_intensity.size() - cur_start;
			feature_length = cur_lenth;
		}

		for (int cur_mz_num : mz_used_for_calc_delay) {
			vector<double> feature_product_intensity_cur_mz(
					feature_product_intensity.at(cur_mz_num).begin() + feature_start,
					feature_product_intensity.at(cur_mz_num).begin() + feature_start
							+ feature_length);
			vector<double> cur_product_intensity_cur_mz(
					cur_product_intensity.at(cur_mz_num).begin() + cur_start,
					cur_product_intensity.at(cur_mz_num).begin() + cur_start + cur_lenth);
			//offset=b-a
			vector<pair<int, double> > cur_delay_score;
			ToolsSpectrum::get_best_delay(feature_product_intensity_cur_mz,
					cur_product_intensity_cur_mz, MAX_ALLOWED_DELAY, cur_delay_score);

			for (auto& t : cur_delay_score) {
				auto id = delay_score.find(t.first);
				if (id == delay_score.end()) {
					delay_score.insert(
							pair<int, double>(t.first,
									t.second * all_possible_mz_intensity.at(cur_mz_num).second));
				} else {
					id->second += t.second * all_possible_mz_intensity.at(cur_mz_num).second;
				}
			}
		}

		//Find max delay
		int new_corrected_delay = 0;
		double new_max_score = 0;
		for (auto& t : delay_score) {
			if (t.second > new_max_score) {
				new_corrected_delay = t.first;
				new_max_score = t.second;
			}
		}
		corrected_offset = new_corrected_delay + (cur_start - feature_start);
	}

////////////////////////// Correct peak /////////////////////////////////
//Correct peak
	if (corrected_offset >= feature_precursor_intensity.size()
			or 0 - corrected_offset >= cur_precursor_intensity.size()) {
		precursor_ions.at(cur_exp_num).set_peak_rt_in_nearby_no(feature_peak_start_in_nearby_no, 0);
		return;
	}

//Make the start equal
	int start = feature_peak_start_in_nearby_no;
	if (corrected_offset > 0) {
		vector<double> temp = vector<double>(cur_precursor_intensity.begin() + corrected_offset,
				cur_precursor_intensity.end());
		cur_precursor_intensity.swap(temp);
	} else {
		vector<double> temp = vector<double>(feature_precursor_intensity.begin() - corrected_offset,
				feature_precursor_intensity.end());
		feature_precursor_intensity.swap(temp);
		start += corrected_offset;
	}

//Make the size equal
	int length = feature_peak_length;
	if (cur_precursor_intensity.size() < feature_precursor_intensity.size()) {
		feature_precursor_intensity.resize(cur_precursor_intensity.size());
	} else {
		cur_precursor_intensity.resize(feature_precursor_intensity.size());
	}
	if (start + length > cur_precursor_intensity.size()) {
		length = cur_precursor_intensity.size() - start;
	}

	if (feature_precursor_intensity.size() != cur_precursor_intensity.size()) {
		cout << "intensity not equal" << endl;
	}

/////////////////////// Determint peak ///////////////////////////////////////
//Get the range to calc
	int cur_peak_start_no = start, cur_peak_length = length;
	if (length <= 0 or start < 0 or length + start > cur_precursor_intensity.size()) {
		if (corrected_offset > 0)
			cur_peak_start_no += corrected_offset;
		if (cur_peak_start_no < 0)
			cur_peak_start_no = 0;
		else if (cur_peak_start_no >= cur_precursor_intensity.size())
			cur_peak_start_no = cur_precursor_intensity.size() - 1;
		cur_peak_length = 0;
	} else {
		if (corrected_offset > 0)
			cur_peak_start_no += corrected_offset;
	}
	precursor_ions.at(cur_exp_num).set_peak_rt_in_nearby_no(cur_peak_start_no, cur_peak_length);
}

void TransitionGroup::get_precursor_rt_in_nearby_no(int exp_num, int& peak_start_no,
		int& peak_length) const {
	precursor_ions[exp_num].get_peak_rt_in_nearby_no(peak_start_no, peak_length);
}

void TransitionGroup::get_precursor_rt_no(int exp_num, int& peak_start_no, int& peak_length) const {
	precursor_ions[exp_num].get_peak_rt_no(peak_start_no, peak_length);
}
}/* namespace GroupDIA */
