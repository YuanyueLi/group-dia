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

#include "TransitionSwathScore.h"

namespace GroupDIA {

TransitionSwathScore::TransitionSwathScore(TransitionIdentified& t) :
		trans(t) {
}

void TransitionSwathScore::generate_score_tabfile(ToolsTabFile& score) {
	map<string, string> map_info;
	map<string, double> map_score;
	get_transition_info(map_info);
	get_transition_score(map_score);

	for (auto& it : map_info) {
		auto jt = map_score.find(it.first);
		if (jt != map_score.end())
			map_score.erase(jt);
	}

//#pragma omp critical
	{
		for (auto it = map_info.begin(); it != map_info.end(); ++it)
			score.add_item(it->first, it->second);
		for (auto it = map_score.begin(); it != map_score.end(); ++it)
			score.add_item(it->first, it->second);
	}
}

void TransitionSwathScore::generate_intensity_tabfile(ToolsTabFile& intensity) {
	vector<string> name;
	string acc = trans.get_accession();
	AASequence aaseq;
	trans.get_pep_aaseq(aaseq);

	//Cal the intensity
	vector<vector<double> > all_intensity;
	vector<string> product_type;
	vector<double> product_mz;
	vector<double> product_th_mz;
	trans.get_product_type(product_type);
	trans.get_product_theoretical_mz(product_th_mz);
	trans.get_product_mz(product_mz);
	for (int cur_exp_num = 0; cur_exp_num < trans.get_exp_size(); cur_exp_num++) {
		vector<double> product_intensity;
		trans.get_product_intensity(product_intensity, cur_exp_num);
		all_intensity.push_back(product_intensity);
	}

	vector<string> cur_value;

	map<string, string> map_info;
	trans.get_information(map_info);
	for (auto& it : map_info) {
		name.push_back(it.first);
		cur_value.push_back(it.second);
	}

	for (int i = 0; i < trans.get_exp_size(); i++) {
		name.push_back("exp_" + String(i));
	}

	for (int cur_exp_num = 0; cur_exp_num < trans.get_exp_size(); cur_exp_num++) {
		cur_value.push_back(String(accumulate(all_intensity[cur_exp_num].begin(), all_intensity[cur_exp_num].end(), 0.0)));
	}
//#pragma omp critical
	for (int i = 0; i < name.size(); i++)
		intensity.add_item(name[i], cur_value[i]);
}

void TransitionSwathScore::cal_distance_score() {
	vector<vector<double> > product_ions_intensity_in_exp(exp_size); //exp, ions
	for (int i = 0; i < exp_size; i++) {
		trans.get_product_intensity(product_ions_intensity_in_exp.at(i), i);
	}
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
				double dist = abs(product_ions_intensity_in_exp.at(cur_exp).at(start) - product_ions_intensity_in_exp.at(cur_exp).at(end));
				product_ions_distance += dist * dist;
				product_ions_weighted_distance += dist * dist * exp_weight.at(cur_exp) * exp_weight.at(cur_exp);
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

	add_score("var_product_dist", product_ions_distance);
	add_score("var_product_dist_weight", product_ions_weighted_distance);
}

void TransitionSwathScore::cal_max_intensity_score(vector<double>& max_intensity) {
	max_intensity.clear();
	trans.get_product_max_intensity(max_intensity);
	double max_intensity_sum = accumulate(max_intensity.begin(), max_intensity.end(), 0.0);
	add_score("var_max_intensity", log10(max_intensity_sum));
}

void TransitionSwathScore::cal_sn_intensity_score(vector<double>& exp_sum_intensity) {
	exp_sum_intensity.clear();
	vector<double> exp_mean_sn;
	for (int i = 0; i < exp_size; i++) {
		if (!trans.is_exist(i)) {
			exp_sum_intensity.push_back(0);
			exp_mean_sn.push_back(0);
			exist_exp_size--;
			continue;
		}

		vector<double> cur_exp_intensity;
		vector<double> cur_exp_sn;
		trans.get_product_sn(cur_exp_sn, i);
		trans.get_product_intensity(cur_exp_intensity, i);

		double cur_intensity = accumulate(cur_exp_intensity.begin(), cur_exp_intensity.end(), 0.0);
		exp_sum_intensity.push_back(cur_intensity);

		double cur_sn = accumulate(cur_exp_sn.begin(), cur_exp_sn.end(), 0.0) / (product_size * 1.0);
		exp_mean_sn.push_back(cur_sn);
	}

	add_score("exp_size", exist_exp_size);

	/*
	 if (exist_exp_size == 0) {
	 //		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
	 //				"exist_exp_size is zero. ID:" + String(this->trans.get_uid()));
	 //		cout << "exist_exp_size is zero. ID:" + String(this->trans.get_uid()) << endl;
	 add_score("var_log_mean_intensity", 0);
	 add_score("var_mean_sn", 0);
	 return;
	 }

	 double sum_intensity = accumulate(exp_sum_intensity.begin(), exp_sum_intensity.end(), 0.0) / (exist_exp_size * 1.0);
	 add_score("var_log_mean_intensity", log(sum_intensity) / exist_exp_size);
	 double mean_sn = accumulate(exp_mean_sn.begin(), exp_mean_sn.end(), 0.0) / (exist_exp_size * 1.0);
	 double var_mean_sn = mean_sn / exist_exp_size;
	 add_score("var_mean_sn", var_mean_sn);
	 */
}

void TransitionSwathScore::cal_similar_score() {
	vector<vector<double> > exp_product_intensity;
	//Get product intensity
	for (int cur_exp = 0; cur_exp < exp_size; cur_exp++) {
		vector<double> product_intensity(product_size);
		if (trans.is_exist(cur_exp)) {
			trans.get_product_intensity(product_intensity, cur_exp);
		}
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
		var_mean_product_intensity_dist = log(dist_product_intensity / (exist_exp_size * product_size));
	else
		var_mean_product_intensity_dist = -10;
	if (dist_exp_intensity != 0)
		dist_exp_intensity = log(dist_exp_intensity / (exist_exp_size * product_size));
	else
		dist_exp_intensity = -10;

	add_score("var_mean_product_intensity_dist", dist_product_intensity);
	add_score("var_mean_exp_intensity_dist", dist_exp_intensity);
}

void inline TransitionSwathScore::add_score(string name, double value) {
//	score.add_score(name, value);
	map_score[name] = value;
}

void TransitionSwathScore::get_transition_score(map<string, double>& _map_score) {
	this->map_score.clear();
	add_score("ide_score", trans.get_identification_score());
//	add_score("var_peak_quality", trans.get_quality());

	exp_size = trans.get_exp_size();
	product_size = trans.get_product_num();
	exist_exp_size = exp_size;
//	factor_exp_size = exp_size;

//	exist_exp_size = 1;
//	factor_exp_size = 1;
	vector<double> max_intensity;
	vector<double> exp_sum_intensity;

//****************************************************************************
//cal the max_intensity score
	cal_max_intensity_score(max_intensity);

//****************************************************************************
//Cal the intensity and SN
	cal_sn_intensity_score(exp_sum_intensity);

//****************************************************************************
	//Cal the weight for rest cal
	exp_weight = exp_sum_intensity;
	ToolsSpectrum::normalize_sum(exp_weight);
	ion_weight = max_intensity;
	ToolsSpectrum::normalize_sum(ion_weight);

//****************************************************************************
//cal the co-elution score
//	cal_co_elution_score();

//****************************************************************************
//cal the similar of the intensity
	cal_similar_score();

//****************************************************************************
//Cal the distance between product ions
	cal_distance_score();

//****************************************************************************
//Cal the b y ion numbers
//	cal_ions_score();

//****************************************************************************
//Cal the mw difference
//	cal_mz_difference_score();

	_map_score = this->map_score;
}

void TransitionSwathScore::get_transition_info(map<string, string>& map_info) {
//	trans.get_information(map_info);
	int decoy = trans.get_decoy_type();
	if (decoy == TransitionIdentified::Decoy_Rand)
		map_info.insert(pair<string, string>("decoy", "1"));
	else
		map_info.insert(pair<string, string>("decoy", "0"));

	map_info.insert(pair<string, string>("ide_score", String(trans.get_identification_score())));
//	map_info.insert(pair<string, string>("transition_group_id", trans.get_pep_seq() + "_" + String(trans.get_precursor_charge())));
	map_info.insert(pair<string, string>("id", String(trans.get_uid())));
	map_info.insert(pair<string, string>("transition_group_id", String(trans.get_uid())));
}

}/* namespace GroupDIA */
