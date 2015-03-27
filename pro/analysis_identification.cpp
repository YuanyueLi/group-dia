/*************************************************************************
Group-DIA
Copyright (C) 2015 .

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*************************************************************************/

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/KERNEL/Peak2D.h>
#include <AnalysisIdentification.h>
#include <ToolsFilename.h>

using namespace std;
using namespace OpenMS;
using namespace GroupDIA;
class Analysis_Identification: public TOPPBase {
public:
	Analysis_Identification() :
			TOPPBase("Analysis Identification", "Analysis database searching result.", false) {
	}
protected:
	void registerOptionsAndFlags_() {
		registerInputFile_("in", "<file>", "", "Input parameter list");
		setValidFormats_("in", StringList::create("txt"));

		registerInputFile_("id", "<file>", "", "Input idenfitication result");
		setValidFormats_("id", StringList::create("pepxml"));

		registerIntOption_("win", "", 0, "Cur windows number", false);

//Parameter
		registerSubsection_("algorithm", "Algorithm parameters section");
	}

	Param getSubsectionDefaults_(const String &) const {
		Param defaults_;

		Param spectrum_ms1;
		spectrum_ms1.setValue("use_ppm", "true");
		spectrum_ms1.setValue("delta_mz", 0.05);
		spectrum_ms1.setValue("delta_ppm", 50);
		defaults_.insert("spectrum_ms1:", spectrum_ms1);

		Param spectrum;
		spectrum.setValue("use_ppm", "false");
		spectrum.setValue("delta_mz", 0.05);
		spectrum.setValue("delta_ppm", 50);
		defaults_.insert("spectrum:", spectrum);

		Param identified;
		identified.setValue("min_product_ions_num", 4);
		identified.setValue("best_flyer_num", 6);
		identified.setValue("min_allowed_product_ions_in_calc_intensity", 2);
		identified.setValue("wanted_min_product_mz", 400.0);
		identified.setValue("max_similar_score", 10.0);
		defaults_.insert("identification:", identified);

		//Para method for decoy
		defaults_.setValue("library_method", "correlation",
				"method for select best flyer ('max','correlation')");
		defaults_.setValue("method", "shuffle",
				"decoy generation method ('random','shuffle','pseudo-reverse','reverse','shift')");
		defaults_.setValue("mz_threshold", 0.8,
				"MZ threshold in Thomson for fragment ion annotation");
		defaults_.setValue("similarity_threshold", 0.05,
				"Similarity threshold for absolute difference of the product mz of target and decoy assays for exclusion in Dalton. Suggested value: 0.05");
		defaults_.setValue("identity_threshold", 0.7,
				"shuffle: identity threshold for the shuffle algorithm");
		defaults_.setValue("max_attempts", 10,
				"shuffle: maximum attempts to lower the sequence identity between target and decoy for the shuffle algorithm");
		defaults_.setValue("mz_shift", 20, "shift: MZ shift in Thomson for shift decoy method");

		defaults_.setValue("exact_scan_num", 30, "");
		defaults_.setValue("stop_report_after_feature", -1,
				"Stop reporting after feature (ordered by quality; -1 means do not stop).");
		defaults_.setValue("max_allowed_delta_rt", 60,
				"a value of 150 means to extract around +/- 150 s of the expected elution. For this to work, the TraML input file needs to contain normalized RT values.");
		defaults_.setValue("rt_normalization_factor", 1.0,
				"The normalized RT is expected to be between 0 and 1. If your normalized RT has a different range, pass this here (e.g. it goes from 0 to 100, set this value to 100)");
		defaults_.setValue("quantification_cutoff", 0.0,
				"Cutoff below which peaks should not be used for quantification any more",
				StringList::create("advanced"));
		defaults_.setMinFloat("quantification_cutoff", 0.0);
		defaults_.setValue("write_convex_hull", "false",
				"Whether to write out all points of all features into the featureXML",
				StringList::create("advanced"));
		defaults_.setValidStrings("write_convex_hull", StringList::create("true,false"));
		defaults_.setValue("add_up_spectra", 1,
				"Add up spectra around the peak apex (needs to be a non-even integer)",
				StringList::create("advanced"));
		defaults_.setMinInt("add_up_spectra", 1);
		defaults_.setValue("spacing_for_spectra_resampling", 0.005,
				"If spectra are to be added, use this spacing to add them up",
				StringList::create("advanced"));
		defaults_.setMinFloat("spacing_for_spectra_resampling", 0.0);

		defaults_.insert("TransitionGroupPicker:", MRMTransitionGroupPicker().getDefaults());
		defaults_.insert("DIAScoring:", DIAScoring().getDefaults());
		defaults_.insert("EMGScoring:", EmgScoring().getDefaults());

		// One can turn on / off each score individually
		Param scores_to_use;
		scores_to_use.setValue("use_shape_score", "true", "Use the shape score",
				StringList::create("advanced"));
		scores_to_use.setValidStrings("use_shape_score", StringList::create("true,false"));
		scores_to_use.setValue("use_coelution_score", "true", "Use the coelution score",
				StringList::create("advanced"));
		scores_to_use.setValidStrings("use_coelution_score", StringList::create("true,false"));
		scores_to_use.setValue("use_rt_score", "true", "Use the retention time score",
				StringList::create("advanced"));
		scores_to_use.setValidStrings("use_rt_score", StringList::create("true,false"));
		scores_to_use.setValue("use_library_score", "true", "Use the library score",
				StringList::create("advanced"));
		scores_to_use.setValidStrings("use_library_score", StringList::create("true,false"));
		scores_to_use.setValue("use_elution_model_score", "true",
				"Use the elution model (EMG) score", StringList::create("advanced"));
		scores_to_use.setValidStrings("use_elution_model_score", StringList::create("true,false"));
		scores_to_use.setValue("use_intensity_score", "true", "Use the intensity score",
				StringList::create("advanced"));
		scores_to_use.setValidStrings("use_intensity_score", StringList::create("true,false"));
		scores_to_use.setValue("use_nr_peaks_score", "true", "Use the number of peaks score",
				StringList::create("advanced"));
		scores_to_use.setValidStrings("use_nr_peaks_score", StringList::create("true,false"));
		scores_to_use.setValue("use_total_xic_score", "true", "Use the total XIC score",
				StringList::create("advanced"));
		scores_to_use.setValidStrings("use_total_xic_score", StringList::create("true,false"));
		scores_to_use.setValue("use_sn_score", "true", "Use the SN (signal to noise) score",
				StringList::create("advanced"));
		scores_to_use.setValidStrings("use_sn_score", StringList::create("true,false"));
		defaults_.insert("Scores:", scores_to_use);

		return defaults_;
	}

	typedef OpenMS::MSExperiment<Peak1D> MapType;
	ExitCodes main_(int, const char **) {
		typedef OpenMS::MSExperiment<Peak1D> MapType;
//******************************* Deal filename *******************************
		String in = getStringOption_("in");
		String in_id = getStringOption_("id");
		int cur_windows_num = getIntOption_("win");

		Param para = getParam_();

		ToolsFilename project;
		int exp_size = project.read_filelist(in);
		int windows_size = project.get_swath_windows_size();
		string project_name = project.get_project_name();

//******************************* Start analysis *******************************
		AnalysisIdentification score(para.copy("algorithm:", true), windows_size, cur_windows_num,
				exp_size);
		vector<String> ms1_filename(exp_size), ms2_filename(exp_size), ms2_part_chro_filename(
				exp_size);
		string input_chro_filename = project.get_ms2_chro_filename(-1, cur_windows_num);
		if (!File::exists(input_chro_filename)) {
			// If file not exist, merge the small files
			fstream fo(input_chro_filename, fstream::out | fstream::binary);
			for (int i = 0; i < exp_size; i++) {
				if (!File::exists(project.get_ms2_chro_filename(i, cur_windows_num))) {
					throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__,
							project.get_ms2_chro_filename(i, cur_windows_num));
				}
				fstream fi(project.get_ms2_chro_filename(i, cur_windows_num),
						fstream::in | fstream::binary);
				char e;
				while (true) {
					fi.read((char*) &e, sizeof(e));
					if (fi.eof())
						break;

					fo.write((char*) &e, sizeof(e));
				}
			}
			fo.close();
		}

		string ms1_part_chro_filename = project.get_identification_chro_filename(cur_windows_num);
		for (int i = 0; i < exp_size; i++) {
			ms2_part_chro_filename[i] = project.get_identification_splited_chro_filename(i,
					cur_windows_num);
			ms1_filename.at(i) = project.get_ms1_mzML_filename(i);
			ms2_filename.at(i) = project.get_ms2_mzML_filename(i, cur_windows_num);
		}
		score.reexact_intensity(in_id, input_chro_filename, ms1_part_chro_filename, ms1_filename,
				ms2_part_chro_filename, ms2_filename,
				project.get_identification_result_chro_filename(cur_windows_num));
		score.score_transition(input_chro_filename, ms1_part_chro_filename, ms2_part_chro_filename,
				ms2_filename, project.get_identification_result_chro_filename(cur_windows_num),
				project.get_identification_result(cur_windows_num));

		return EXECUTION_OK;
	}
};

int main(int argc, const char ** argv) {
	Analysis_Identification tool;
	return tool.main(argc, argv);
}

/// @endcond
