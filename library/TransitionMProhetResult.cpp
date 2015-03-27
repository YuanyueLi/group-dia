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

#include "TransitionMProhetResult.h"

namespace GroupDIA {

} /* namespace GroupDIA */

void GroupDIA::MProhetResult::add_intensity(const String& ori_value) {
	vector<String> s;
	ori_value.split(",", s);
	intensity.clear();
	intensity.reserve(s.size());
	for (auto it = s.begin(); it != s.end(); ++it) {
		if (*it != "") {
			intensity.push_back(it->toDouble());
		}
	}
}

void GroupDIA::TransitionMProhetResult::add_result(int exp_num, MProhetResult& result) {
	int exist_num = 0;
	for (auto it = result.intensity.begin(); it != result.intensity.end(); ++it) {
		if (*it > 0)
			exist_num++;
	}
	if (exist_num < 2)
		return;

	auto it = tran.find(exp_num);
	if (it != tran.end())
		cout << "Error in add result: " << exp_num << endl;

	tran.insert(pair<int, MProhetResult>(exp_num, result));
}

double GroupDIA::TransitionMProhetResult::get_intensity(int exp_num) const {
	auto it = tran.find(exp_num);
	if (it == tran.end())
		return 0;
	return accumulate(it->second.intensity.begin(), it->second.intensity.end(), 0.0);
}

double GroupDIA::TransitionMProhetResult::get_average_d_score() const {
	return 0;
//	if (size() == 0)
//		return 0;
//	double d_score = 0;
//	for (auto it = tran.begin(); it != tran.end(); ++it) {
//		d_score += it->second.d_score;
//	}
//	return d_score / size();
}

double GroupDIA::TransitionMProhetResult::get_ide_score() const {
	return tran.begin()->second.ide_score;
}

int GroupDIA::TransitionMProhetResult::size() const {
	return tran.size();
}

void GroupDIA::TransitionMProhetResult::get_experiment_dist_score(map<string, double>& score) const {
	int product_size = tran.begin()->second.intensity.size();
	vector<vector<double> > exp_product_intensity;
	//Get product intensity
	for (auto it = tran.begin(); it != tran.end(); ++it) {
		vector<double> product_intensity;
		product_intensity = it->second.intensity;
		exp_product_intensity.push_back(product_intensity);
	}
	int exist_exp_size = exp_product_intensity.size();

	const vector<vector<double> >& intensity = exp_product_intensity;

	//Calc exp_intensity
	vector<double> exp_intensity(exist_exp_size, 0);
	for (int cur_exp = 0; cur_exp < exist_exp_size; cur_exp++) {
		for (int cur_product = 0; cur_product < product_size; cur_product++) {
			exp_intensity[cur_exp] += intensity[cur_exp][cur_product];
		}
	}

	//cal product_intensity
	vector<double> product_intensity(product_size, 0);
	for (int cur_product = 0; cur_product < product_size; cur_product++) {
		for (int cur_exp = 0; cur_exp < exist_exp_size; cur_exp++) {
			product_intensity[cur_product] += intensity[cur_exp][cur_product];
		}
	}

	//Calc exp_distance
	ToolsSpectrum::normalize_max(exp_intensity);
	double dist_exp_intensity = 0;
	for (int cur_product = 0; cur_product < product_size; cur_product++) {
		vector<double> temp(exist_exp_size);
		for (int cur_exp = 0; cur_exp < exist_exp_size; cur_exp++) {
			temp[cur_exp] = intensity[cur_exp][cur_product];
		}
		ToolsSpectrum::normalize_max(temp);
		for (int cur_exp = 0; cur_exp < exist_exp_size; cur_exp++)
			dist_exp_intensity += abs(temp[cur_exp] - exp_intensity[cur_exp]);
	}

	//Calc product_distance
	ToolsSpectrum::normalize_max(product_intensity);
	double dist_product_intensity = 0;
	for (int cur_exp = 0; cur_exp < exist_exp_size; cur_exp++) {
		vector<double> temp(product_size);
		for (int cur_product = 0; cur_product < product_size; cur_product++) {
			temp[cur_product] = intensity[cur_exp][cur_product];
		}
		ToolsSpectrum::normalize_max(temp);
		for (int cur_product = 0; cur_product < product_size; cur_product++)
			dist_product_intensity += abs(temp[cur_product] - product_intensity[cur_product]);
	}

	//Record
	double var_mean_product_intensity_dist = 0, var_mean_exp_intensity_dist = 0;
	if (dist_product_intensity != 0)
		var_mean_product_intensity_dist = log(
				dist_product_intensity / (exist_exp_size * product_size));
	else
		var_mean_product_intensity_dist = -10;
	if (dist_exp_intensity != 0)
		dist_exp_intensity = log(dist_exp_intensity / (exist_exp_size * product_size));
	else
		dist_exp_intensity = -10;

	score.insert(pair<string, double>("var_mean_product_intensity_dist", dist_product_intensity));
	score.insert(pair<string, double>("var_mean_exp_intensity_dist", dist_exp_intensity));
}

void GroupDIA::TransitionMProhetResult::remove_interference_ions() {
	//Build intensity matrix
	int tran_size = tran.size();
	if (tran_size < 3)
		return;
	int ion_size = tran.begin()->second.intensity.size();
	bool exist_interferenct = true;
	int remove_cycle = 0;

	while (exist_interferenct and remove_cycle <= tran_size) {
		exist_interferenct = false;
		remove_cycle++;
		vector<vector<double> > matrix(tran_size);
		auto tt = matrix.begin();
		for (auto it = tran.begin(); it != tran.end(); ++it, ++tt) {
			*tt = it->second.intensity;
		}

		//Normalize
		for (auto it = matrix.begin(); it != matrix.end(); ++it) {
			double sum = accumulate(it->begin(), it->end(), 0.0);
			if (sum != 0) {
				for (auto& jt : *it) {
					if (jt < sum) {
						jt = jt / (sum - jt);
					} else {
						jt = INT_MAX;
					}
				}
			} else {
				for (auto& jt : *it) {
					jt = 0;
				}
			}
		}

		//Find outlier by median + 3 * mad
		for (int j = 0; j < ion_size; j++) {
			double median = 0, sad = 0;
			vector<double> value;
			for (int i = 0; i < tran_size; i++) {
				if (matrix[i][j] > 0) {
					value.push_back(matrix[i][j]);
				}
			}
			if (value.size() < 3)
				continue;
			ToolsSpectrum::median_and_sad(value, median, sad);
			if (sad == 0)
				continue;
			double cut_off = median + 3 * sad;
			for (int i = 0; i < tran_size; i++) {
				double& x = matrix[i][j];
				if (x >= cut_off) {
					//Is Outlier
					exist_interferenct = true;
					//Replace values with mean.
					auto cur_tran = tran.begin();
					advance(cur_tran, i);
					double& origin_value = cur_tran->second.intensity[j];
					double origin_sum = accumulate(cur_tran->second.intensity.begin(),
							cur_tran->second.intensity.end(), 0.0);
					sort(value.begin(), value.end());
					double median = value[value.size() / 2];
					origin_value = (origin_sum - origin_value) * median;
				}
			}
		}
	}
}

void GroupDIA::TransitionMProhetResult::get_distance_score(map<string, double>& score) const {
	int product_size = tran.begin()->second.intensity.size();
	vector<vector<double> > product_ions_intensity_in_exp;
	vector<double> exp_weight;
//Get product intensity
	for (auto it = tran.begin(); it != tran.end(); ++it) {
		vector<double> product_intensity;
		product_intensity = it->second.intensity;
		product_ions_intensity_in_exp.push_back(product_intensity);
		exp_weight.push_back(accumulate(product_intensity.begin(), product_intensity.end(), 0.0));
	}
	int exp_size = product_ions_intensity_in_exp.size();

//Find the max exp
	int max_exp = 0;
	double max_value = 0;
	for (int i = 0; i < exp_size; i++) {
		if (exp_weight.at(i) > max_value) {
			max_value = exp_weight.at(i);
			max_exp = i;
		}
	}
//Normalize
	for (int i = 0; i < product_size; i++) {
		double max = product_ions_intensity_in_exp.at(max_exp).at(i);
		for (int j = 0; j < exp_size; j++) {
			if (max > 0) {
				product_ions_intensity_in_exp.at(j).at(i) /= max;
			} else {
				product_ions_intensity_in_exp.at(j).at(i) = 0;
			}
		}
	}

//Cal the distance
	double product_ions_distance = 0;
	double product_ions_weighted_distance = 0;
	for (int cur_exp = 0; cur_exp < exp_size; cur_exp++) {
		for (int start = 0; start < product_size; start++) {
			for (int end = start + 1; end < product_size; end++) {
				double dist = abs(
						product_ions_intensity_in_exp.at(cur_exp).at(start)
								- product_ions_intensity_in_exp.at(cur_exp).at(end));
				product_ions_distance += dist * dist;
				product_ions_weighted_distance += dist * dist * exp_weight.at(cur_exp)
						* exp_weight.at(cur_exp);
			}
		}
	}
	product_ions_distance = sqrt(product_ions_distance);
	product_ions_weighted_distance = sqrt(product_ions_weighted_distance);

	product_ions_distance /= exp_size * (product_size * (product_size - 1));
	product_ions_weighted_distance /= exp_size * (product_size * (product_size - 1));

	if (product_ions_distance != 0)
		product_ions_distance = log(product_ions_distance);
	else
		product_ions_distance = -10;

	if (product_ions_weighted_distance != 0)
		product_ions_weighted_distance = log(product_ions_weighted_distance);
	else
		product_ions_weighted_distance = -10;

	score.insert(pair<string, double>("var_product_dist", product_ions_distance));
	score.insert(pair<string, double>("var_product_dist_weight", product_ions_weighted_distance));
}

double GroupDIA::TransitionMProhetResult::get_max_p() const {
	if (size() == 0)
		return 0;
	double max_p = 0;
	for (auto it = tran.begin(); it != tran.end(); ++it) {
		if (max_p < it->second.p)
			max_p = it->second.p;
	}

	return max_p;
}
