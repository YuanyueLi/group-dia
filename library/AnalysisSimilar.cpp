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

#include "AnalysisSimilar.h"

//const double MIN_INTENSITY = 50;
//const double MIN_ALLOWED_COR = 0.5;
//const int MAX_IONS_NUM = 25;
//const int AVERGAE_NUM_IN_A_CLUSTER = 5;

bool large_first(const Candidate& i, const Candidate& j) {
	return (i.cor > j.cor);
}

template<class T1, class T2> bool small_first(const pair<T1, T2>& i, const pair<T1, T2>& j) {
	return (i.first < j.first);
}

namespace GroupDIA {
void AnalysisSimilar::get_intensity_score(const vector<vector<double> >& product_ms2_intensity,
		vector<double>& score) const {
	score.clear();
	for (int i = 0; i < product_ms2_intensity.size(); i++) {
		if (product_ms2_intensity.at(i).size() <= 0)
			score.push_back(0);
		else
			score.push_back(
					*max_element(product_ms2_intensity.at(i).begin(),
							product_ms2_intensity.at(i).end()));
	}
}

void AnalysisSimilar::generate_result_exp(const TransitionGroup& trans,
		const vector<pair<double, double> >& select_info, MSSpectrum<Peak1D>& spectrum) const {
	OpenMS::Precursor pre;
	pre.setCharge(trans.get_precursor_charge());
	pre.setMZ(trans.get_precursor_mz());
	vector<OpenMS::Precursor> v_pre;
	v_pre.push_back(pre);
	spectrum.setPrecursors(v_pre);

	vector<double> mz;
	vector<double> intensity;
	mz.reserve(select_info.size()), intensity.reserve(select_info.size());
	for (auto it = select_info.begin(); it != select_info.end(); ++it) {
		Peak1D ion;
		ion.setMZ(it->first);
		ion.setIntensity(it->second);
		spectrum.push_back(ion);
	}

	spectrum.setMSLevel(2);
	spectrum.sortByPosition();
}

void AnalysisSimilar::cal_similar_by_kmeans(const TransitionGroup& transition,
		MSSpectrum<Peak1D>& target_spectrum, MSSpectrum<Peak1D>& decoy_spectrum) const {
	vector<double> precurosr_intensity;
	vector<vector<double> > product_ms2_intensity;
	int peak_start, peak_length;
	transition.get_product_intensity_in_each_time_nearby_version(product_ms2_intensity);
	transition.get_precursor_peak_intensity(precurosr_intensity);

	//Get the intensity
	vector<double> intensity;
	get_intensity_score(product_ms2_intensity, intensity);

	//Get score
	double MIN_SCORE = -2;
	vector<double> score;
	ToolsSpectrum::get_pearson_score(precurosr_intensity, product_ms2_intensity, score);

	for (int i = 0; i < intensity.size(); i++) {
		if (intensity[i] < MIN_INTENSITY)
			score[i] = MIN_SCORE;
	}

//Get product mz
	vector<double> mz;
	transition.get_product_mz(mz);

	/////////////////////////
	//Calculate target product
	/////////////////////////

	vector<pair<double, double> > target_product;
	//Generate result by calc the distance

	vector<Candidate> all_candi;
	vector<Candidate> selected_candi;
	for (int i = 0; i < mz.size(); i++) {
		Candidate c;
		c.ion_mz = mz[i];
		c.ion_intensity = intensity[i];
		c.cor = score[i];
		c.intensity = &(product_ms2_intensity[i]);
		all_candi.push_back(c);
	}
	refine_similar_by_k_means(precurosr_intensity, all_candi, selected_candi);

	target_product.clear();
	for (auto& it : selected_candi) {
		target_product.push_back(pair<double, double>(it.ion_mz, it.ion_intensity));
	}

	target_spectrum.clear(true);
	generate_result_exp(transition, target_product, target_spectrum);

	/////////////////////////
	//Calculate decoy product
	/////////////////////////
	vector<pair<double, double> > decoy_product;
	//Get candidate
	vector<Candidate> decoy_candi;
	for (auto cur_candi : all_candi) {
		bool is_target = false;
		for (auto it : selected_candi) {
			if (cur_candi.ion_mz == it.ion_mz) {
				is_target = true;
				break;
			}
		}
		if (!is_target)
			decoy_candi.push_back(cur_candi);
	}

	//Random select
	ToolsSpectrum::random(decoy_candi.begin(), decoy_candi.end());
	decoy_product.clear();
	for (int i = 0; i < selected_candi.size() and i < decoy_candi.size(); i++) {
		decoy_product.push_back(
				pair<double, double>(decoy_candi[i].ion_mz, decoy_candi[i].ion_intensity));
	}

	generate_result_exp(transition, decoy_product, decoy_spectrum);

	target_spectrum.setRT(transition.get_precursor_rt(transition.get_feature_exp_num()));
	target_spectrum.setNativeID(transition.get_scan_num() * 2 + 0);
	decoy_spectrum.setRT(transition.get_precursor_rt(transition.get_feature_exp_num()));
	decoy_spectrum.setNativeID(transition.get_scan_num() * 2 + 1);
}

void AnalysisSimilar::cal_similar_by_simply_cor(const TransitionSimilar& target,
		MSSpectrum<Peak1D>& target_spectrum, MSSpectrum<Peak1D>& decoy_spectrum) const {
	vector<double> precurosr_intensity;
	vector<vector<double> > product_ms2_intensity;
	target.get_product_intensity_in_each_time_nearby_version(product_ms2_intensity);
	target.get_precursor_peak_intensity(precurosr_intensity);

	//Get the intensity
	vector<double> intensity;
	get_intensity_score(product_ms2_intensity, intensity);

	//Get score
	double MIN_SCORE = -2;
	vector<double> score;
	ToolsSpectrum::get_pearson_score(precurosr_intensity, product_ms2_intensity, score);

	for (int i = 0; i < intensity.size(); i++) {
		if (intensity[i] < MIN_INTENSITY)
			score[i] = MIN_SCORE;
	}

	//Select score cutoff
	vector<double> sorted_score = score;
	double target_cutoff = ToolsSpectrum::get_nth_large(sorted_score, MAX_IONS_NUM);
	double decoy_cutoff = ToolsSpectrum::get_nth_small(sorted_score, MAX_IONS_NUM);
	if (target_cutoff < MIN_ALLOWED_COR)
		target_cutoff = MIN_ALLOWED_COR;
	if (decoy_cutoff > target_cutoff)
		decoy_cutoff = target_cutoff;

	//Get product mz
	vector<double> mz;
	target.get_product_mz(mz);

	//Calc target product
	vector<pair<double, double> > target_product;
	target_product.reserve(MAX_IONS_NUM);
	for (int i = 0; i < score.size(); i++) {
		if (score[i] >= target_cutoff) {
			target_product.push_back(pair<double, double>(mz[i], intensity[i]));
		}
	}

	generate_result_exp(target, target_product, target_spectrum);

	//Calc decoy product
	vector<pair<double, double> > decoy_product;
	decoy_product.reserve(MAX_IONS_NUM);
	for (int i = 0; i < score.size(); i++) {
		if (score[i] <= decoy_cutoff) {
			decoy_product.push_back(pair<double, double>(mz[i], intensity[i]));
		}
	}

	generate_result_exp(target, decoy_product, decoy_spectrum);

	target_spectrum.setRT(target.get_scan_num() * 2 + 0);
	decoy_spectrum.setRT(target.get_scan_num() * 2 + 1);
}

void AnalysisSimilar::k_means(const vector<vector<double> >& data,
		vector<vector<double> >& means_data, vector<int>& cluster) const {
	int cluster_size = data.size() / AVERGAE_NUM_IN_A_CLUSTER;
	if (cluster_size <= 1) {
		means_data.resize(1);
		calc_k_means_center(data, means_data[0]);

		cluster.clear();
		for (int i = 0; i < data.size(); i++)
			cluster.push_back(1);
		return;
	}

	int data_size = data.size();
	const int MAX_IT = 100;

	//First select cluster
	cluster = vector<int>(data_size, -1);
	means_data.clear();
	means_data.reserve(cluster_size);
	for (int i = 0; i < cluster_size; i++)
		means_data.push_back(data[i]);

	for (int cur_it = 0; cur_it < MAX_IT; cur_it++) {
		bool changed = false;
		//Assign to cluster
		for (int i = 0; i < data_size; i++) {
			double max_score = -2;
			int max_cluster = -1;
			for (int c = 0; c < cluster_size; c++) {
				double cur_score = cor(means_data[c], data[i]);
				if (cur_score > max_score) {
					max_score = cur_score;
					max_cluster = c;
				}
			}
			if (cluster[i] != max_cluster) {
				cluster[i] = max_cluster;
				changed = true;
			}
		}

		//Generate means
		for (int c = 0; c < cluster_size; c++) {
			vector<vector<double> > cur_cluster_data;
			for (int i = 0; i < data_size; i++) {
				if (cluster[i] == c) {
					cur_cluster_data.push_back(data[i]);
				}
			}
			if (cur_cluster_data.size() > 0)
				calc_k_means_center(cur_cluster_data, means_data[c]);
		}

		//Remove not exist cluster
		for (int c = 0; c < cluster_size; c++) {
			bool is_exist = false;
			for (int i = 0; i < data_size; i++) {
				if (cluster[i] == c) {
					is_exist = true;
					break;
				}
			}
			if (!is_exist) {
				int to_swap = cluster_size - 1;
				cluster_size--;
				for (int i = 0; i < data_size; i++)
					if (cluster[i] == to_swap)
						cluster[i] = c;
				swap(means_data[c], means_data[to_swap]);

				means_data.resize(cluster_size);
				c--;
			}
		}

		if (!changed) {
			break;
		}
	}
}

void AnalysisSimilar::calc_k_means_center(const vector<vector<double> >& data,
		vector<double>& result) const {
	int data_size = data.size();
	int col_size = data[0].size();
	result.resize(col_size, 0);
	for (int i = 0; i < col_size; i++) {
		for (int j = 0; j < data_size; j++)
			result[i] += data[j][i];
		result[i] /= col_size;
	}
}

void AnalysisSimilar::refine_similar_by_k_means(const vector<double>& precursor_intensity,
		vector<Candidate>& all_ions, vector<Candidate>& selected_result) const {
	selected_result.clear();
	for (Candidate& c : all_ions) {
		if (c.cor > MIN_ALLOWED_COR)
			selected_result.push_back(c);
	}

	if (selected_result.size() < MAX_IONS_NUM) {
		return;
	}

	sort(selected_result.begin(), selected_result.end(), large_first);

	//Cluster by kmeans
	vector<vector<double> > data, means;
	vector<int> clusters;
	for (int i = 0; i < selected_result.size(); i++) {
		data.push_back(*(selected_result[i].intensity));
	}
	k_means(data, means, clusters);

	//Select cluster
	int cluster_size = means.size();
	int data_size = data.size();

	map<double, int, greater<double> > cor_cluster;
	for (int i = 0; i < cluster_size; i++) {
		cor_cluster.insert(pair<double, int>(cor(precursor_intensity, means[i]), i));
	}
	set<int> result;
	for (auto it = cor_cluster.begin(); it != cor_cluster.end(); ++it) {
		if (it->first < MIN_ALLOWED_COR or result.size() >= MAX_IONS_NUM)
			break;
		//Get cur_cluster num
		vector<int> cur_cluster_item;
		for (int i = 0; i < data_size; i++)
			if (clusters[i] == it->second)
				result.insert(i);
	}

	//Generate result;
	vector<Candidate> final_result;
	for (int i = 0; i < data_size; i++)
		if (result.find(i) != result.end())
			final_result.push_back(selected_result[i]);

	selected_result = final_result;
}

AnalysisSimilar::AnalysisSimilar(const Param& para) {
	MIN_INTENSITY = para.getValue("min_allowed_intensity");
	MIN_ALLOWED_COR = para.getValue("min_allowed_correlation");
	MAX_IONS_NUM = para.getValue("min_ion_numbers_in_spectrum");
	AVERGAE_NUM_IN_A_CLUSTER = para.getValue("number_in_a_cluster");
}

inline double AnalysisSimilar::cor(const vector<double>& a, const vector<double>& b) const {
	return Math::pearsonCorrelationCoefficient(a.begin(), a.end(), b.begin(), b.end());
}

}/* namespace GroupDIA */
