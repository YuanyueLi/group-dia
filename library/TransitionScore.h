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

#ifndef TRANSITIONSCORE_H_
#define TRANSITIONSCORE_H_
#include "TransitionIdentified.h"
#include "ToolsSpectrum.h"
#include "ToolsTabFileWrite.h"
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/MRMScoring.h>
#define run_identifier "unique_run_identifier"

namespace GroupDIA {

class TransitionScore {
	typedef MSSpectrum<ChromatogramPeak> SpectrumT;
	typedef OpenSwath::LightTransition TransitionT;
public:
	TransitionScore(TransitionIdentified& t, Param p);

	void score(int cur_exp, OpenSwath::SpectrumAccessPtr swath_ptr, ToolsTabFileWrite& output);

private:

	void score(const TransitionIdentified& tr, int exp_num, OpenSwath::SpectrumAccessPtr swath_ptr,
			FeatureMap<Feature>& output);

	void calculateSwathScores_(MRMTransitionGroup<SpectrumT, TransitionT> & transition_group,
			MRMFeature & mrmfeature_, OpenSwath::SpectrumAccessPtr swath_map,
			std::vector<double>& normalized_library_intensity, OpenSwath_Scores scores);

	OpenSwath::SpectrumPtr getAddedSpectra_(OpenSwath::SpectrumAccessPtr swath_map, double RT,
			int nr_spectra_to_add);

	void updateMembers_(Param param);

	void scorePeakgroups_(MRMTransitionGroup<SpectrumT, TransitionT> & transition_group,
			int exp_num, OpenSwath::SpectrumAccessPtr swath_map, FeatureMap<Feature>& output);

	void cal_similar_score(vector<MRMFeature>::iterator& i_cur_feature, int exp_num);

	void cal_cor_score(vector<MRMFeature>::iterator& i_cur_feature, int exp_num);

	void decide_group_place(vector<MRMFeature>::iterator& i_cur_feature, int exp_num);

	OpenMS::DIAScoring diascoring_;
	OpenMS::EmgScoring emgscoring_;

	// Which scores to use
	// Variables

	// Which scores to use
	bool use_coelution_score_;bool use_shape_score_;bool use_rt_score_;bool use_library_score_;bool use_elution_model_score_;bool use_intensity_score_;bool use_total_xic_score_;bool use_nr_peaks_score_;bool use_sn_score_;

	double MAX_DIFF_IN_SCORE;
	double max_allowed_delta_rt;
	DoubleReal spacing_for_spectra_resampling_;
	int add_up_spectra_;

	// All the filters expect MSSpectrum<PeakT>, thus we give it an "MSSpectrum"
	// but filled with Chromatogram Peaks.
	typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram; // this is the type in which we store the chromatograms for this analysis
	typedef OpenSwath::LightTransition TransitionType;
	typedef OpenSwath::LightTargetedExperiment TargetedExpType;
	typedef OpenSwath::LightPeptide PeptideType;
	typedef OpenSwath::LightProtein ProteinType;
	typedef OpenSwath::LightModification ModificationType;
	typedef MRMTransitionGroup<MSSpectrum<ChromatogramPeak>, TransitionType> MRMTransitionGroupType; // a transition group holds the MSSpectra with the Chromatogram peaks from above
	typedef std::map<String, MRMTransitionGroupType> TransitionGroupMapType;

	TransitionIdentified& tr;

	Param param_;
};

}
/* namespace GroupDIA */

#endif /* TRANSITIONSCORE_H_ */
