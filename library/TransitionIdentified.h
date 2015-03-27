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

#ifndef TRANSITIONIDENTIFIED_H_
#define TRANSITIONIDENTIFIED_H_

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include "TransitionGroup.h"
#include "IntensityCalc.h"
#include "ToolsScore.h"
//#define DEBUG

namespace GroupDIA {
class TransitionIdentified: public TransitionGroup {
public:
	enum {
		Target = 0, Decoy_Spec = 1, Decoy_Iden = 2, Decoy_Rand = 3,
	} DecoyType;

	TransitionIdentified(const Param& para);
	TransitionIdentified();

//Set identified result
	void set_identification_result(const PeptideHit& pephit, bool is_decoy = false);
	void set_product_by_seq(bool is_selected_by_sequence);
	void add_product_mz(const vector<double>& product_mz);

	void reselect_product_ions(const Param& para);
	void calc_product_intensity();
	void choose_best_flyer_by_percursor_pattern(vector<int>& remain);
	void choose_best_flyer_by_random_select(vector<int>& remain) const;
	void choose_best_flyer_by_select_max_intensity(vector<int>& remain);

//Get identified result
	string get_pep_seq() const;
	void get_pep_aaseq(AASequence& seq) const;
	string get_accession() const;
	double get_identification_score() const;

//Decoy
	void generate_decoy_by_simple_copy(TransitionIdentified& decoy) const;
	void generate_decoy_by_openswath_method(TransitionIdentified& decoy_trans, bool& is_right,
			String method, double identity_threshold, int max_attempts, double mz_threshold,
			double mz_shift, double similarity_threshold) const;

	bool is_decoy() const;
	int get_decoy_type() const;
	void set_decoy(int decoy_type);

//Refine the intensity.
	bool get_library_intensity(double& precursor_intensity, vector<double>& product_intensity);
	double calc_dist_score(const vector<double>& ms1, vector<double>& ms2, double lib_ms1,
			double lib_ms2) const;

//Exist or not
	bool is_exist(int exp_num) const;
	int get_exist_exp_num() const;
	bool is_exist() const;

//Generate mrm for openswath score.
	void generate_mrm(int exp_num, const Param& param_,
			MRMTransitionGroup<MSSpectrum<ChromatogramPeak>, OpenSwath::LightTransition>& mrm) const;

//Get basic information
	void get_information(map<string, string>& info) const;

//Get information for score
	void get_product_max_intensity(vector<double>& intensity) const;
	void get_product_type(vector<string>&type) const;
	void get_product_sn(vector<double>& sn, int exp_num) const;

	void get_product_intensity(vector<double>& intensity, int exp_num);
	void get_product_theoretical_mz(vector<double>&mz) const;

//Reset peak
	void reset_peak(int exp_num, double left_peak_rt, double right_peak_rt);
	virtual bool load_store(fstream& f, int io_type, int type = LoadSave::Type::Null, int ion_num =
			0);

private:
	vector<vector<double> > final_product_intensity; //exp_num, product_ions

	map<double, double> map_product_theoretical_mz;
	map<double, string> map_product_type;
	map<double, double> map_product_intensity;

	int BEST_FLYER_NUM;
	int MIN_PRODUCT_IONS_NUM_IN_INTENSITY;
	double MIN_FLYER_MZ;
	int MIN_PRODUCT_IONS_NUM;
	double MAX_DIFF_IN_SCORE;

//Set identified result
	int add_product_ions(RichPeakSpectrum& spectrum, string add, vector<double>& mz);
}
;

} /* namespace GroupDIA */
#endif /* TRANSITIONIDENTIFIED_H_ */
