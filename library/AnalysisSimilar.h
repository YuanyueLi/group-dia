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

#ifndef ANALYSISSIMILAR_H_
#define ANALYSISSIMILAR_H_
#include "TransitionGroup.h"
#include "ToolsTabFile.h"
#include "TransitionSimilar.h"
#include "ToolsSpectrum.h"
#include "TransitionIdentified.h"

namespace GroupDIA {
struct Candidate {
	double ion_mz;
	double ion_intensity;
	double cor;
	double score;
	vector<double>* intensity;
};

class AnalysisSimilar {
public:
	AnalysisSimilar(const Param& para);

	void cal_similar_by_kmeans(const TransitionGroup& target, MSSpectrum<Peak1D>& target_spectrum,
			MSSpectrum<Peak1D>& decoy_spectrum) const;

	void cal_similar_by_simply_cor(const TransitionSimilar& target,
			MSSpectrum<Peak1D>& target_spectrum, MSSpectrum<Peak1D>& decoy_spectrum) const;
private:
	void get_intensity_score(const vector<vector<double> >& product_ms2_intensity,
			vector<double>& score) const;
	void generate_result_exp(const TransitionGroup& trans,
			const vector<pair<double, double> >& select_info, MSSpectrum<Peak1D>& spectrum) const;

	void refine_similar_by_k_means(const vector<double>& precursor_intensity,
			vector<Candidate>& all_ions, vector<Candidate>& selected_result) const;

	void k_means(const vector<vector<double> >& data, vector<vector<double> >& means_data,
			vector<int>& cluster) const;

	void calc_k_means_center(const vector<vector<double> >& data, vector<double>& result) const;

	inline double cor(const vector<double>& a, const vector<double>& b) const;

	double MIN_INTENSITY;
	double MIN_ALLOWED_COR;
	int MAX_IONS_NUM;
	int AVERGAE_NUM_IN_A_CLUSTER;
};

} /* namespace GroupDIA */
#endif /* ANALYSISSIMILAR_H_ */
