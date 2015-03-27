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

#include "TransitionSimilar.h"

//#define DEBUG

namespace GroupDIA {

} /* namespace GroupDIA */

template<class T1, class T2> bool large_first(const pair<T1, T2>& i, const pair<T1, T2>& j) {
	return (i.first > j.first);
}

template<class T1, class T2> bool small_first(const pair<T1, T2>& i, const pair<T1, T2>& j) {
	return (i.first < j.first);
}

void GroupDIA::TransitionSimilar::refined_by_crosscorrelation(int MAX_ALLOWED_DELAY, int FINE_NEARBY,
		double MIN_INTENSITY) {
	int exp_size = this->get_exp_size();
	int mz_size = this->product_mz.size();

//Get feature precursor and product intensity
	vector<double> ori_feature_precursor_intensity; //intensity
	vector<vector<double> > feature_product_intensity; //mz, intensity
	int feature_peak_start_in_nearby_no, feature_peak_length;
	feature_precursor.get_peak_rt_in_nearby_no(feature_peak_start_in_nearby_no,
			feature_peak_length);
	feature_precursor.get_nearby_intensity(ori_feature_precursor_intensity);
	product_ions.at(feature_exp_num).get_intensity(feature_product_intensity);

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
	for (int cur_exp_num = 0; cur_exp_num < exp_size; cur_exp_num++) {
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
			precursor_ions.at(cur_exp_num).set_peak_rt_in_nearby_no(feature_peak_start_in_nearby_no,
					0);
			continue;
		}

		//Make the start equal
		int start = feature_peak_start_in_nearby_no;
		if (corrected_offset > 0) {
			vector<double> temp = vector<double>(cur_precursor_intensity.begin() + corrected_offset,
					cur_precursor_intensity.end());
			cur_precursor_intensity.swap(temp);
		} else {
			vector<double> temp = vector<double>(
					feature_precursor_intensity.begin() - corrected_offset,
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
			//Determine which peak is the wanted peak
//			determine_peak(feature_precursor_intensity.begin(), feature_precursor_intensity.end(),
//					cur_precursor_intensity.begin(), cur_precursor_intensity.end(), start, length, cur_peak_start_no,
//					cur_peak_length);
			if (corrected_offset > 0)
				cur_peak_start_no += corrected_offset;
		}
		precursor_ions.at(cur_exp_num).set_peak_rt_in_nearby_no(cur_peak_start_no, cur_peak_length);
	}
}

void GroupDIA::TransitionSimilar::crude_refine(double MIN_ALLOWED_CORRELATION) {
//Get precursor and product ions intensity
	vector<double> precursor_intensity;
	vector<vector<double> > product_intensity;

	feature_precursor.get_peak_intensity(precursor_intensity);
	int peak_rt_from, peak_length;
	feature_precursor.get_peak_rt_in_nearby_no(peak_rt_from, peak_length);
	feature_product.get_intensity(0, peak_length, product_intensity);

//Refine
	int mz_size = this->product_mz.size();
	vector<double> new_product_mz;
	for (int i = 0; i < mz_size; i++) {
		double cor = Math::pearsonCorrelationCoefficient(precursor_intensity.begin(),
				precursor_intensity.end(), product_intensity.at(i).begin(),
				product_intensity.at(i).end());
		if (cor > MIN_ALLOWED_CORRELATION)
			new_product_mz.push_back(product_mz.at(i));
	}
	product_mz = new_product_mz;
	product_mz.shrink_to_fit();
	feature_product = ProductIons();
}

void GroupDIA::TransitionSimilar::determine_peak(const VecIt ref_intensity_start,
		const VecIt ref_intensity_end, const VecIt sam_intensity_start,
		const VecIt sam_intensity_end, int ref_peak_start, int ref_peak_length, int& sam_peak_start,
		int& sam_peak_length) const {
//Find peak
	vector<PeakInfo> peak_info;
	vector<double> ori_intensity(sam_intensity_start, sam_intensity_end);
//	ToolsSpectrum::smooth_intensity(ori_intensity);
	ToolsSpectrum::peak_finder(ori_intensity, sam_peak_start, sam_peak_length, peak_info);

//Test peak one by one
	double max_cor = -2;
	PeakInfo candidate_peak(sam_peak_start, 0);
	peak_info.push_back(candidate_peak);
	for (auto ipeak = peak_info.begin(); ipeak != peak_info.end(); ++ipeak) {
		//Generate ref intensity
		vector<double> ref(ref_intensity_start, ref_intensity_end);
		for (auto it = ref.begin() + ref_peak_start;
				it != ref.begin() + ref_peak_start + ref_peak_length; ++it)
			*it = 0;
		vector<double> sam(sam_intensity_start, sam_intensity_end);
		for (auto it = sam.begin() + ipeak->first; it != sam.begin() + ipeak->first + ipeak->second;
				++it)
			*it = 0;
		double cur_cor = Math::pearsonCorrelationCoefficient(ref.begin(), ref.end(), sam.begin(),
				sam.end());
		if (cur_cor > max_cor) {
			max_cor = cur_cor;
			candidate_peak = *ipeak;
		}
	}

	if (candidate_peak.second > 0) {
		sam_peak_start = candidate_peak.first;
		sam_peak_length = candidate_peak.second;
	} else {
		vector<int> rt(ref_peak_length);
		for (int i = 0; i < ref_peak_length; i++)
			rt.at(i) = i;
		//Get the apex of the ref peak
		double apex = ToolsSpectrum::cal_peak_apex(ref_intensity_start + ref_peak_start,
				ref_intensity_start + ref_peak_start + ref_peak_length, rt.begin(), rt.end());
		sam_peak_start = ref_peak_start + apex;
		sam_peak_length = 1;
	}
}
