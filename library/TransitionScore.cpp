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

#include "TransitionScore.h"
namespace GroupDIA {
/// Score all peak groups
void TransitionScore::scorePeakgroups_(
		MRMTransitionGroup<SpectrumT, TransitionT> & transition_group, int exp_num,
		OpenSwath::SpectrumAccessPtr swath_map, FeatureMap<Feature>& output) {
	//std::vector<SignalToNoiseEstimatorMedian<RichPeakChromatogram> > signal_noise_estimators;
	typedef typename MRMTransitionGroup<SpectrumT, TransitionT>::PeakType PeakT;
	std::vector<OpenSwath::ISignalToNoisePtr> signal_noise_estimators;
	std::vector<MRMFeature> feature_list;

	DoubleReal sn_win_len_ = (DoubleReal) param_.getValue("TransitionGroupPicker:sn_win_len");
	DoubleReal sn_bin_count_ = (DoubleReal) param_.getValue("TransitionGroupPicker:sn_bin_count");
	for (Size k = 0; k < transition_group.getChromatograms().size(); k++) {
		OpenSwath::ISignalToNoisePtr snptr(
				new OpenMS::SignalToNoiseOpenMS<PeakT>(transition_group.getChromatograms()[k],
						sn_win_len_, sn_bin_count_));
		signal_noise_estimators.push_back(snptr);
	}

	// get the expected rt value for this peptide
	double expected_rt = tr.get_precursor_rt(exp_num);

// Go through all peak groups (found MRM features) and score them
	for (std::vector<MRMFeature>::iterator mrmfeature =
			transition_group.getFeaturesMuteable().begin();
			mrmfeature != transition_group.getFeaturesMuteable().end(); mrmfeature++) {
		int group_size = boost::numeric_cast<int>(transition_group.size());
		if (group_size == 0) {
			throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
					"Error: Transition group " + transition_group.getTransitionGroupID()
							+ " has no chromatograms.");
		}
		if (group_size < 2) {
			std::cerr << "Error: transition group " << transition_group.getTransitionGroupID()
					<< " has less than 2 chromatograms. It has " << group_size << std::endl;
			continue;
		}
		if (tr.get_product_rt_size(exp_num, (double) mrmfeature->getMetaValue("leftWidth"),
				(double) mrmfeature->getMetaValue("rightWidth")) < 2)
			continue;
		if (abs(mrmfeature->getRT() - expected_rt) > max_allowed_delta_rt)
			continue;

		OpenSwath::IMRMFeature* imrmfeature;
		imrmfeature = new MRMFeatureOpenMS(*mrmfeature);

		OpenSwath::ITransitionGroup* itransition_group;
		itransition_group = new TransitionGroupOpenMS<SpectrumT, TransitionT>(transition_group);

		// calculate the normalized library intensity (expected value of the intensities)
		std::vector<double> normalized_library_intensity;
		transition_group.getLibraryIntensity(normalized_library_intensity);
		OpenSwath::Scoring::normalize_sum(&normalized_library_intensity[0],
				boost::numeric_cast<int>(normalized_library_intensity.size()));

		// calcxcorr -> for each lag do the correlation, normally use lag 0
		// xcorr_matrix  => correlate chromatogram i with chromatogram j

		OpenSwath::MRMScoring mrmscore_;
		bool normalize = true;
		mrmscore_.initializeXCorrMatrix(imrmfeature, itransition_group, normalize);

		// XCorr score (coelution)
		double xcorr_coelution_score = 0;
		if (use_coelution_score_) {
			xcorr_coelution_score = mrmscore_.calcXcorrCoelutionScore();
			mrmfeature->addScore("var_xcorr_coelution", xcorr_coelution_score);
		}

		double weighted_coelution_score = 0;
		if (use_coelution_score_) {
			weighted_coelution_score = mrmscore_.calcXcorrCoelutionScore_weighted(
					normalized_library_intensity);
			mrmfeature->addScore("var_xcorr_coelution_weighted ", weighted_coelution_score);
		}

		// XCorr score (shape)
		// mean over the intensities at the max of the crosscorrelation
		// FEATURE : weigh by the intensity as done by mQuest
		// FEATURE : normalize with the intensity at the peak group apex?
		double xcorr_shape_score = 0;
		if (use_shape_score_) {
			xcorr_shape_score = mrmscore_.calcXcorrShape_score();
			mrmfeature->addScore("var_xcorr_shape", xcorr_shape_score);
		}

		double weighted_xcorr_shape = 0;
		if (use_shape_score_) {
			weighted_xcorr_shape = mrmscore_.calcXcorrShape_score_weighted(
					normalized_library_intensity);
			mrmfeature->addScore("var_xcorr_shape_weighted", weighted_xcorr_shape);
		}

		// FEATURE : how should we best calculate correlation between library and experiment?
		// FEATURE : spectral angle
		double library_corr = 0, library_rmsd = 0;
		double library_manhattan, library_dotprod;
		if (use_library_score_) {
			mrmscore_.calcLibraryScore(imrmfeature, transition_group.getTransitions(), library_corr,
					library_rmsd, library_manhattan, library_dotprod);
			mrmfeature->addScore("var_library_corr", library_corr);
			mrmfeature->addScore("var_library_rmsd", library_rmsd);
			mrmfeature->addScore("var_library_manhattan", library_manhattan); // new score
			mrmfeature->addScore("var_library_dotprod", library_dotprod); // new score
		}

		// Retention time score
		double rt_score = 0, norm_rt_score = 0;
		if (use_rt_score_) {
			// get the id, then get the expected and the experimental retention time
			String native_id = transition_group.getChromatograms()[0].getNativeID();
			TransitionType tr = transition_group.getTransition(native_id);
			double experimental_rt = mrmfeature->getFeature(native_id).getRT();
			rt_score = abs(experimental_rt - expected_rt);
			norm_rt_score = rt_score / max_allowed_delta_rt;
			mrmfeature->addScore("delta_rt", mrmfeature->getRT() - expected_rt);
			mrmfeature->addScore("assay_rt", expected_rt);
			mrmfeature->addScore("norm_RT", experimental_rt);
			mrmfeature->addScore("rt_score", rt_score);
			mrmfeature->addScore("var_norm_rt_score", norm_rt_score);
		}
		mrmfeature->addScore("rt", mrmfeature->getRT());

		// Intensity score
		double intensity_score = 0;
		if (use_intensity_score_) {
			intensity_score = mrmfeature->getIntensity()
					/ (double) mrmfeature->getMetaValue("total_xic");
			mrmfeature->addScore("var_intensity_score", intensity_score);
		}

		double total_xic_score = 0;
		if (use_total_xic_score_) {
			total_xic_score = (double) mrmfeature->getMetaValue("total_xic");
			mrmfeature->addScore("total_xic", total_xic_score);
		}

		double nr_peaks_score = 0;
		if (use_nr_peaks_score_) {
			nr_peaks_score = group_size;
			mrmfeature->addScore("nr_peaks", nr_peaks_score);
		}

		double sn_score = 0, log_sn_score = 0;
		if (use_sn_score_) {
			sn_score = mrmscore_.calcSNScore(imrmfeature, signal_noise_estimators);
			if (sn_score < 1) { // fix to make sure, that log(sn_score = 0) = -inf does not occur
				log_sn_score = 0;
			} else {
				log_sn_score = std::log(sn_score);
			}
			mrmfeature->addScore("sn_ratio", sn_score);
			mrmfeature->addScore("var_log_sn_score", log_sn_score);
		}

		OpenSwath_Scores scores;
		double quick_lda_dismiss = 0;
		double lda_quick_score = -scores.get_quick_lda_score(library_corr, library_rmsd,
				norm_rt_score, xcorr_coelution_score, xcorr_shape_score, log_sn_score);

		if (lda_quick_score < quick_lda_dismiss) {
			// continue;
		}

		double elution_model_fit_score = 0;
		if (use_elution_model_score_) {
#pragma omp critical
			elution_model_fit_score = emgscoring_.calcElutionFitScore((*mrmfeature),
					transition_group);
			mrmfeature->addScore("var_elution_model_fit_score", elution_model_fit_score);
		}

		//Add the similar score
		cal_similar_score(mrmfeature, exp_num);
		cal_cor_score(mrmfeature, exp_num);
		decide_group_place(mrmfeature, exp_num);

		double xx_lda_prescore;
		scores.library_corr = library_corr;
		scores.library_rmsd = library_rmsd;
		scores.norm_rt_score = norm_rt_score;
		scores.elution_model_fit_score = elution_model_fit_score;
		scores.log_sn_score = log_sn_score;
		scores.xcorr_coelution_score = xcorr_coelution_score;
		scores.xcorr_shape_score = xcorr_shape_score;
		xx_lda_prescore = -scores.calculate_lda_prescore(scores);

		bool swath_present = (swath_map->getNrSpectra() > 0);
		if (!swath_present) {
			mrmfeature->addScore("main_var_xx_lda_prelim_score", xx_lda_prescore);
			mrmfeature->setOverallQuality(xx_lda_prescore);
		} else {
			mrmfeature->addScore("xx_lda_prelim_score", xx_lda_prescore);
		}

		if (swath_present) {
//#pragma omp critical
			calculateSwathScores_(transition_group, *mrmfeature, swath_map,
					normalized_library_intensity, scores);
		}

		///////////////////////////////////////////////////////////////////////////
		// add the peptide hit information to the feature
		///////////////////////////////////////////////////////////////////////////

		PeptideIdentification pep_id_ = PeptideIdentification();
		PeptideHit pep_hit_ = PeptideHit();

		pep_hit_.setCharge(tr.get_precursor_charge());
		pep_hit_.setScore(xx_lda_prescore);
		AASequence aas;
		tr.get_pep_aaseq(aas);
		pep_hit_.setSequence(aas);
		pep_hit_.addProteinAccession(tr.get_accession());
		pep_id_.insertHit(pep_hit_);
		pep_id_.setIdentifier(run_identifier);

		mrmfeature->getPeptideIdentifications().push_back(pep_id_);
		mrmfeature->setMetaValue("PrecursorMZ",
				transition_group.getTransitions()[0].getPrecursorMZ());
		mrmfeature->setSubordinates(mrmfeature->getFeatures()); // add all the subfeatures as subordinates
		double total_intensity = 0, total_peak_apices = 0;
		for (std::vector<Feature>::iterator sub_it = mrmfeature->getSubordinates().begin();
				sub_it != mrmfeature->getSubordinates().end(); sub_it++) {
			total_intensity += sub_it->getIntensity();
			total_peak_apices += (DoubleReal) sub_it->getMetaValue("peak_apex_int");
		}
		// overwrite the reported intensities with those above the m/z cutoff
		mrmfeature->setIntensity(total_intensity);
		mrmfeature->setMetaValue("peak_apices_sum", total_peak_apices);
		feature_list.push_back((*mrmfeature));

		delete imrmfeature;
		delete itransition_group;
	}

	// Order by quality
	std::sort(feature_list.begin(), feature_list.end(), OpenMS::Feature::OverallQualityLess());
	std::reverse(feature_list.begin(), feature_list.end());

	for (Size i = 0; i < feature_list.size(); i++) {
		output.push_back(feature_list[i]);
	}
}

OpenSwath::SpectrumPtr TransitionScore::getAddedSpectra_(OpenSwath::SpectrumAccessPtr swath_map,
		double RT, int nr_spectra_to_add) {
	std::vector<std::size_t> indices = swath_map->getSpectraByRT(RT, 0.0);
	int closest_idx = boost::numeric_cast<int>(indices[0]);
	if (indices[0] != 0
			&& std::fabs(
					swath_map->getSpectrumMetaById(boost::numeric_cast<int>(indices[0]) - 1).RT
							- RT)
					< std::fabs(
							swath_map->getSpectrumMetaById(boost::numeric_cast<int>(indices[0])).RT
									- RT)) {
		closest_idx--;
	}

	if (nr_spectra_to_add == 1) {
		OpenSwath::SpectrumPtr spectrum_ = swath_map->getSpectrumById(closest_idx);
		return spectrum_;
	} else {
		std::vector<OpenSwath::SpectrumPtr> all_spectra;
		// always add the spectrum 0, then add those right and left
		all_spectra.push_back(swath_map->getSpectrumById(closest_idx));
		for (int i = 1; i <= nr_spectra_to_add / 2; i++) { // cast to int is intended!
			all_spectra.push_back(swath_map->getSpectrumById(closest_idx - i));
			all_spectra.push_back(swath_map->getSpectrumById(closest_idx + i));
		}
		OpenSwath::SpectrumPtr spectrum_ = SpectrumAddition::addUpSpectra(all_spectra,
				spacing_for_spectra_resampling_, true);
		return spectrum_;
	}
}

void TransitionScore::calculateSwathScores_(
		MRMTransitionGroup<SpectrumT, TransitionT> & transition_group, MRMFeature & mrmfeature_,
		OpenSwath::SpectrumAccessPtr swath_map, std::vector<double>& normalized_library_intensity,
		OpenSwath_Scores scores) {
	MRMFeature* mrmfeature = &mrmfeature_;

	// parameters
	int by_charge_state = 1; // for which charge states should we check b/y series

	// find spectrum that is closest to the apex of the peak using binary search
	OpenSwath::SpectrumPtr spectrum_ = getAddedSpectra_(swath_map, mrmfeature->getRT(),
			add_up_spectra_);
	OpenSwath::SpectrumPtr* spectrum = &spectrum_;

	// Isotope correlation / overlap score: Is this peak part of an
	// isotopic pattern or is it the monoisotopic peak in an isotopic
	// pattern?
	OpenSwath::IMRMFeature* imrmfeature = new MRMFeatureOpenMS(*mrmfeature);
	double isotope_corr = 0, isotope_overlap = 0;
	diascoring_.dia_isotope_scores(transition_group.getTransitions(), (*spectrum), imrmfeature,
			isotope_corr, isotope_overlap);
	mrmfeature->addScore("var_isotope_correlation_score", isotope_corr);
	mrmfeature->addScore("var_isotope_overlap_score", isotope_overlap);

	double bseries_score = 0, yseries_score = 0, massdev_score = 0;
	if (param_.getValue("method") != "random") {
		// Mass deviation score
		double ppm_score = 0, ppm_score_weighted = 0;
		diascoring_.dia_massdiff_score(transition_group.getTransitions(), (*spectrum),
				normalized_library_intensity, ppm_score, ppm_score_weighted);
		// FEATURE we should not punish so much when one transition is missing!
		massdev_score = ppm_score / transition_group.size();
		double massdev_score_weighted = ppm_score_weighted;
		mrmfeature->addScore("var_massdev_score", massdev_score);
		mrmfeature->addScore("var_massdev_score_weighted", massdev_score_weighted);

		// Presence of b/y series score
		OpenMS::AASequence aas;
//	OpenSwathDataAccessHelper::convertPeptideToAASequence(*PeptideRefMap_[transition_group.getTransitions()[0].getPeptideRef()], aas);
		tr.get_pep_aaseq(aas);
		diascoring_.dia_by_ion_score((*spectrum), aas, by_charge_state, bseries_score,
				yseries_score);
		mrmfeature->addScore("var_bseries_score", bseries_score);
		mrmfeature->addScore("var_yseries_score", yseries_score);
	}

	double dotprod_score_dia;
	double manhatt_score_dia;
	diascoring_.score_with_isotopes((*spectrum), transition_group.getTransitions(),
			dotprod_score_dia, manhatt_score_dia);
	mrmfeature->addScore("var_dotprod_score", dotprod_score_dia);
	mrmfeature->addScore("var_manhatt_score", manhatt_score_dia);

	scores.yseries_score = yseries_score;
	scores.isotope_correlation = isotope_corr;
	scores.isotope_overlap = isotope_overlap;
	scores.massdev_score = massdev_score;
	double xx_swath_prescore = -scores.calculate_swath_lda_prescore(scores);
	mrmfeature->addScore("main_var_xx_swath_prelim_score", xx_swath_prescore);
	mrmfeature->setOverallQuality(xx_swath_prescore);
	delete imrmfeature;
}

void TransitionScore::cal_cor_score(vector<MRMFeature>::iterator& i_cur_feature, int exp_num) {
	int product_num = tr.get_product_num();
	if (product_num == 0)
		cout << "Product_num is zero" << endl;

	//Get origin intensity
	vector<double> origin_precursor_intensity;
	tr.get_precursor_nearby_intensity(exp_num, origin_precursor_intensity);

	vector<vector<double> > origin_product_intensity;
	tr.get_product_intensity_in_each_time(exp_num, origin_product_intensity);
	if (product_num != origin_product_intensity.size())
		cout << "Product_num not equal!" << endl;
	if (origin_product_intensity[0].size() != origin_precursor_intensity.size())
		cout << "Precursor num is not equal to product num" << endl;

	//Get peak information
	double left_rt = i_cur_feature->getMetaValue("leftWidth");
	double right_rt = i_cur_feature->getMetaValue("rightWidth");
	int peak_left_rt_no = tr.get_product_rt_no(exp_num, left_rt);
	int peak_right_rt_no = tr.get_product_rt_no(exp_num, right_rt) + 1;
	int peak_length = peak_right_rt_no - peak_left_rt_no;

	//Get precurosr intensity
	vector<double> precursor_intensity(origin_precursor_intensity.begin() + peak_left_rt_no,
			origin_precursor_intensity.begin() + peak_left_rt_no + peak_length);

	//Get product intensity
	vector<vector<double> > product_intensity(product_num);
	for (int i = 0; i < product_num; i++) {
		product_intensity[i] = vector<double>(origin_product_intensity[i].begin() + peak_left_rt_no,
				origin_product_intensity[i].begin() + peak_left_rt_no + peak_length);
	}

	//Calc cor score
	double score = 0;
	for (int i = 0; i < product_num; i++) {
		double cur_score = Math::pearsonCorrelationCoefficient(precursor_intensity.begin(),
				precursor_intensity.end(), product_intensity[i].begin(),
				product_intensity[i].end());
		if (cur_score != cur_score)
			cur_score = -1;
		score += cur_score;
	}
	i_cur_feature->addScore("var_cor_score", score / product_num);

}

void TransitionScore::cal_similar_score(vector<MRMFeature>::iterator& i_cur_feature, int exp_num) {
	int product_num = tr.get_product_num();
	if (product_num == 0)
		cout << "Product_num is zero" << endl;

	//Get origin intensity
	vector<double> origin_precursor_intensity;
	tr.get_precursor_nearby_intensity(exp_num, origin_precursor_intensity);

	vector<vector<double> > origin_product_intensity;
	tr.get_product_intensity_in_each_time(exp_num, origin_product_intensity);
	if (product_num != origin_product_intensity.size())
		cout << "Product_num not equal!" << endl;
	if (origin_product_intensity[0].size() != origin_precursor_intensity.size())
		cout << "Precursor num is not equal to product num" << endl;

	double similar_score = MAX_DIFF_IN_SCORE;
	//Get library intensity
	double lib_precursor_intensity;
	vector<double> lib_product_intensity;
	bool is_lib_exist = tr.get_library_intensity(lib_precursor_intensity, lib_product_intensity);
	double& sum_lib_precursor_intensity = lib_precursor_intensity;
	double sum_lib_product_intensity = accumulate(lib_product_intensity.begin(),
			lib_product_intensity.end(), 0.0);
	if (!is_lib_exist) {
		i_cur_feature->addScore("var_similar_score", similar_score);
		return;
	}

	//Get peak information
	double left_rt = i_cur_feature->getMetaValue("leftWidth");
	double right_rt = i_cur_feature->getMetaValue("rightWidth");
	int peak_left_rt_no = tr.get_product_rt_no(exp_num, left_rt);
	int peak_right_rt_no = tr.get_product_rt_no(exp_num, right_rt) + 1;
	int peak_length = peak_right_rt_no - peak_left_rt_no;

	//Get precurosr intensity
	double sum_precursor_intensity = accumulate(
			origin_precursor_intensity.begin() + peak_left_rt_no,
			origin_precursor_intensity.begin() + peak_left_rt_no + peak_length, 0.0);

	//Get product intensity
	double sum_product_intensity = 0;
	vector<vector<double> > product_intensity(product_num);
	for (int i = 0; i < product_num; i++) {
		sum_product_intensity = accumulate(origin_product_intensity[i].begin() + peak_left_rt_no,
				origin_product_intensity[i].begin() + peak_left_rt_no + peak_length,
				sum_product_intensity);
	}

	//Calc similar score
	if (!(sum_precursor_intensity == 0 or sum_product_intensity == 0
			or sum_lib_precursor_intensity == 0 or sum_lib_product_intensity == 0))
		similar_score = abs(
				log(sum_precursor_intensity / sum_product_intensity)
						- log(sum_lib_precursor_intensity / sum_lib_product_intensity));
	i_cur_feature->addScore("var_similar_score", similar_score);
}

TransitionScore::TransitionScore(TransitionIdentified& t, Param p) :
		tr(t) {
	updateMembers_(p);
}

void TransitionScore::updateMembers_(Param param) {
	param_ = param;

	param_.copy("TransitionGroupPicker:", true);

	max_allowed_delta_rt = (DoubleReal) param_.getValue("max_allowed_delta_rt");
	MAX_DIFF_IN_SCORE = (DoubleReal) param_.getValue("identification:max_similar_score");

	diascoring_.setParameters(param_.copy("DIAScoring:", true));
	emgscoring_.setFitterParam(param_.copy("EmgScoring:", true));

	add_up_spectra_ = param_.getValue("add_up_spectra");
	spacing_for_spectra_resampling_ = param_.getValue("spacing_for_spectra_resampling");
	max_allowed_delta_rt = (DoubleReal) param_.getValue("max_allowed_delta_rt");

	use_coelution_score_ = param_.getValue("Scores:use_coelution_score").toBool();
	use_shape_score_ = param_.getValue("Scores:use_shape_score").toBool();
	use_rt_score_ = param_.getValue("Scores:use_rt_score").toBool();
	use_library_score_ = param_.getValue("Scores:use_library_score").toBool();
	use_elution_model_score_ = param_.getValue("Scores:use_elution_model_score").toBool();
	use_intensity_score_ = param_.getValue("Scores:use_intensity_score").toBool();
	use_total_xic_score_ = param_.getValue("Scores:use_total_xic_score").toBool();
	use_nr_peaks_score_ = param_.getValue("Scores:use_nr_peaks_score").toBool();
	use_sn_score_ = param_.getValue("Scores:use_sn_score").toBool();

//	rt_normalization_factor_ = 600;

//	use_coelution_score_ = param_.getValue("Scores:use_coelution_score").toBool();
//	use_shape_score_ = param_.getValue("Scores:use_shape_score").toBool();
//	use_library_score_ = param_.getValue("Scores:use_library_score").toBool();
//	use_elution_model_score_ = param_.getValue("Scores:use_elution_model_score").toBool();
//	use_intensity_score_ = param_.getValue("Scores:use_intensity_score").toBool();
//	use_total_xic_score_ = param_.getValue("Scores:use_total_xic_score").toBool();
//	use_nr_peaks_score_ = param_.getValue("Scores:use_nr_peaks_score").toBool();
//	use_sn_score_ = param_.getValue("Scores:use_sn_score").toBool();
}

void TransitionScore::score(int cur_exp, OpenSwath::SpectrumAccessPtr swath_ptr,
		ToolsTabFileWrite& output) {
	map<string, string> info;
	info.insert(
			pair<string, string>("decoy",
					(tr.get_decoy_type() == TransitionIdentified::Decoy_Rand) ? "1" : "0"));
	info.insert(pair<string, string>("decoy_type", String(tr.get_decoy_type())));
	tr.get_information(info);

	FeatureMap<Feature> cur_result;
	score(tr, cur_exp, swath_ptr, cur_result);

	if (cur_result.empty())
		return;

	//Rentent time
	double tr_rt = tr.get_precursor_rt(cur_exp);
	double min_delta_rt = abs(cur_result[0].getMetaValue("rt").toString().toDouble() - tr_rt);
	auto min_peak = cur_result.begin();
	for (auto it = cur_result.begin(); it != cur_result.end(); ++it) {
		double rt = it->getMetaValue("rt").toString().toDouble();
		if (abs(tr_rt - rt) < min_delta_rt) {
			min_delta_rt = abs(tr_rt - rt);
			min_peak = it;
		}
	}

	vector<String> names;
	cur_result[0].getKeys(names);
#pragma omp critical
	for (auto it = cur_result.begin(); it != cur_result.end(); ++it) {
		//Add openswath score
		for (auto iname = names.begin(); iname != names.end(); ++iname) {
			output.add_item(*iname, it->getMetaValue(*iname).toString());
		}
		//Add transition intensity
		string aggr_intensity = "";
		for (auto jt = it->getSubordinates().begin(); jt != it->getSubordinates().end(); ++jt) {
			aggr_intensity += String(jt->getIntensity()) + ",";
		}
		output.add_item("aggr_Peak_Area", aggr_intensity);

		//Add transition group id
		output.add_item("transition_exp_id", String(tr.get_uid()));
		output.add_item("transition_group_id",
				String(tr.get_uid()) + String("_") + String(cur_exp));

		//Add other information
		output.add_item("intensity", it->getIntensity());

		//Add transition information
		output.add_item("exp_num", cur_exp);
		for (auto& it : info) {
			output.add_item(it.first, it.second);
		}
		output.end_add_item();
	}
}

void TransitionScore::score(const TransitionIdentified& tr, int exp_num,
		OpenSwath::SpectrumAccessPtr swath_ptr, FeatureMap<Feature>& output) {
	MRMTransitionGroup<SpectrumT, TransitionT> transition_group;

	tr.generate_mrm(exp_num, param_.copy("TransitionGroupPicker:", true), transition_group);

//	MapType swath_map;
//	OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);

	scorePeakgroups_(transition_group, exp_num, swath_ptr, output);

}

void TransitionScore::decide_group_place(vector<MRMFeature>::iterator& i_cur_feature, int exp_num) {
	int ingroup = 0;
	//Get peak information
	double left_rt = i_cur_feature->getMetaValue("leftWidth");
	double right_rt = i_cur_feature->getMetaValue("rightWidth");
	double apex_rt = i_cur_feature->getMetaValue("norm_RT");
	int peak_left_rt_no = tr.get_product_rt_no(exp_num, left_rt);
	int peak_right_rt_no = tr.get_product_rt_no(exp_num, right_rt) + 1;
	int apex_rt_no = tr.get_product_rt_no(exp_num, apex_rt);

	int group_peak_start_no, group_peak_length;
	tr.get_precursor_rt_in_nearby_no(exp_num, group_peak_start_no, group_peak_length);
	int group_peak_end_no = group_peak_start_no + group_peak_length;

	if (apex_rt_no > group_peak_start_no and apex_rt_no < group_peak_end_no)
		ingroup = 1;
	if (peak_right_rt_no > group_peak_start_no + group_peak_length / 2
			and peak_left_rt_no <= group_peak_start_no)
		ingroup = 1;
	if (peak_left_rt_no <= group_peak_start_no + group_peak_length / 2
			and peak_right_rt_no > group_peak_end_no)
		ingroup = 1;

	i_cur_feature->addScore("ingroup", ingroup);
	//Add selected information
	if (ingroup == 1 and exp_num == tr.get_feature_exp_num()) {
		i_cur_feature->addScore("selected", 1);
	} else
		i_cur_feature->addScore("selected", 0);
}
} /* namespace GroupDIA */
