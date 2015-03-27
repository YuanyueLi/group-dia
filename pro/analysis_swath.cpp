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
#include <omp.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <AnalysisSWATH.h>
#include <ToolsFilename.h>
#include <stdio.h>
#include <MGFFile.h>
#include <MzMLFile.h>
#include <thread>
using namespace std;
using namespace OpenMS;
using namespace GroupDIA;
class Analysis_SWATH: public TOPPBase {
public:
	Analysis_SWATH() :
			TOPPBase("Analysis SWATH", "Analysis SWATH-MS data files to generate pseudo-spectra.", false) {
	}

protected:
	void registerOptionsAndFlags_() {
		registerInputFile_("in", "<file>", "", "Input parameter list");
		setValidFormats_("in", StringList::create("txt"));

		registerIntOption_("win", "", 0, "Cur windows number", false);
		registerIntOption_("run", "", -1,
				"The run number to analysis, -1 means analysis all runs at a time", false);

		registerDoubleOption_("max", "", 100, "Max deal kilo features in one time", false);

		registerTOPPSubsection_("ms1", "");
		registerStringOption_("ms1:use_ppm", "", "true", "use ppm in analysis ms1", false);
		registerDoubleOption_("ms1:delta_ppm", "", 50, "the ppm window", false);
		registerDoubleOption_("ms1:delta_mz", "", 0.05, "the mz window", false);

		registerTOPPSubsection_("ms2", "");
		registerStringOption_("ms2:use_ppm", "", "false", "use ppm in analysis ms2", false);
		registerDoubleOption_("ms2:delta_ppm", "", 50, "the ppm window", false);
		registerDoubleOption_("ms2:delta_mz", "", 0.05, "the mz window", false);

		registerTOPPSubsection_("similar", "");
		registerDoubleOption_("similar:min_allowed_intensity", "", 50.0, "", false);
		registerDoubleOption_("similar:min_allowed_correlation", "", 0.5, "", false);
		registerIntOption_("similar:min_ion_numbers_in_spectrum", "", 25, "", false);
		registerIntOption_("similar:number_in_a_cluster", "", 5, "", false);

		registerIntOption_("similar:min_allowed_product_ion_numbers", "", 4, "", false);
		registerIntOption_("similar:retention_window", "", 100, "", false);
		registerDoubleOption_("similar:min_allowed_correlation_in_crude_refine", "", 0.3, "",
				false);

		registerTOPPSubsection_("delay", "");
		registerIntOption_("delay:max_allowd_scan_num_delay_in_normalized", "", 60, "", false);
		registerIntOption_("delay:max_allowd_scan_num_delay_in_second_normalized", "", 10, "",
				false);
		registerDoubleOption_("delay:min_intensity_in_calc_delay", "", 300.0, "", false);
		//Param:
		//AnalysisSWATH
		//AnalysisSimilar
		//TransitionSimilar
	}

	ExitCodes main_(int, const char **) {
//******************************* Deal para *******************************
		String in = getStringOption_("in");
		int win_num = getIntOption_("win");
		int run_num = getIntOption_("run");
		int max_analysis_num = getDoubleOption_("max") * 1000;

		Param para = getParam_();
		Param ms1_para = para.copy("ms1:", true);
		Param ms2_para = para.copy("ms2:", true);

		ToolsFilename project;
		int project_size = project.read_filelist(in);
		int windows_size = project.get_swath_windows_size();
		string project_name = project.get_project_name();

		double mz_min, mz_max;
		project.get_mz(win_num, mz_min, mz_max);

		if (project_size < 1) {
			cout << "Filenum error:" << endl;
			return EXECUTION_OK;
		}
//******************************* Start analysis *******************************
		//Set the ms1
		int ref_exp_num = 0;
		int max_ref_exp_num = project_size;
		if (run_num != -1) {
			ref_exp_num = run_num;
			max_ref_exp_num = run_num + 1;
		}

		File::remove(project.get_ms2_chro_filename(run_num, win_num));
		File::remove(project.get_feature_chro_filename(win_num));
		int cur_feature_num = 0;
		int split_num = 0;
		static int scan_num = 0;
		vector<MSSpectrum<Peak1D> > result;

		while (true) {
			AnalysisSWATH swath(para);

			while (ref_exp_num < max_ref_exp_num) {
				FeatureXMLFile feature_file;
				FeatureMap<> features;
				if (File::exists(project.get_ms1_feature_filename(ref_exp_num))) {
					feature_file.load(project.get_ms1_feature_filename(ref_exp_num), features);
				} else {
					cout << project.get_ms1_feature_filename(ref_exp_num) << " is not exist."
							<< endl;
					ref_exp_num++;
					continue;
				}

				//******************************* Start set ms1 *******************************
				cout << "Start set feature" << endl;
				SwathExperiment ms1_exp(ms1_para);
				ms1_exp.set_exp(project.get_ms1_mzML_filename(ref_exp_num));

//				for (; cur_feature_num < features.size(); cur_feature_num++) {
#pragma omp parallel
				{
					while (cur_feature_num < features.size() and swath.size() < max_analysis_num) {
						int cur_num;
#pragma omp critical
						{
							cur_num = cur_feature_num;
							cur_feature_num++;
							scan_num++;
						}
//						if (features[cur_num].getRT() < 2930 or features[cur_num].getRT() > 2950)
//							continue;

						if (cur_num < features.size() and features[cur_num].getMZ() > mz_min
								and features[cur_num].getMZ() < mz_max) {
							swath.set_ms1(features[cur_num], ref_exp_num,
									(cur_num * project_size + ref_exp_num) * windows_size + win_num,
									ms1_exp);
						}
					}
				}
				if (swath.size() >= max_analysis_num)
					break;
				scan_num = 0;
				cur_feature_num = 0;
				ref_exp_num++;
			}

			if (swath.size() == 0)
				break;

			//******************************* Start set ms2 *******************************
			cout << "Start get m/z of product ions" << endl;
#ifdef QUICK
			SwathExperiment* temp_peak_exp;
			for (int ref_exp_num = 0; ref_exp_num < project_size; ref_exp_num++) {
				SwathExperiment* exp;
				SwathExperiment* peak_exp;
				thread *t, *t_peak;

				if (ref_exp_num == 0) {
					exp = new SwathExperiment(ms2_para);
					exp->set_exp(project.get_ms2_mzML_filename(ref_exp_num, win_num));
					peak_exp = new SwathExperiment(ms2_para);
					peak_exp->set_exp(project.get_ms2_peak_filename(ref_exp_num, win_num));
				} else {
					exp = temp_ms2_exp;
					peak_exp = temp_peak_exp;
				}

				if (ref_exp_num < project_size - 1) {
					//Load next experiment
					temp_ms2_exp = new SwathExperiment(ms2_para);
					t = new thread([temp_ms2_exp,&project,ref_exp_num, win_num]() {
								temp_ms2_exp->set_exp(project.get_ms2_mzML_filename(ref_exp_num+1, win_num));
							});

					temp_peak_exp = new SwathExperiment(ms2_para);
					t_peak = new thread([temp_peak_exp,&project,ref_exp_num, win_num]() {
								temp_peak_exp->set_exp(project.get_ms2_peak_filename(ref_exp_num+1, win_num));
							});
				}

				swath.set_ms2(ref_exp_num, *peak_exp, *exp);

				if (ref_exp_num < project_size - 1) {
					t->join();
					t_peak->join();
					delete t;
					delete t_peak;
				}
				//Delete cur_exp
				delete exp;
				delete peak_exp;
			}
#else
			for (int ref_exp_num = 0; ref_exp_num < project_size; ref_exp_num++) {
				SwathExperiment* exp;

				exp = new SwathExperiment(ms2_para);
				exp->set_exp(project.get_ms2_mzML_filename(ref_exp_num, win_num));

				swath.set_ms2(ref_exp_num, *exp, mz_min, mz_max,
						project.get_feature_chro_filename(win_num));

				//Delete cur_exp
				delete exp;
			}
#endif
			swath.refine();

			cout << "Start get intensity" << endl;
			//Add the ms1
			SwathExperiment* temp_ms1_exp;
			for (int add_exp_num = 0; add_exp_num < project_size; add_exp_num++) {
				SwathExperiment* ms1_exp;
				thread* load_ms1;
				if (add_exp_num == 0) {
					ms1_exp = new SwathExperiment(ms1_para);
					ms1_exp->set_exp(project.get_ms1_mzML_filename(add_exp_num));
				} else {
					load_ms1->join();
					delete load_ms1;
					ms1_exp = temp_ms1_exp;
				}

				SwathExperiment* ms2_exp = new SwathExperiment(ms2_para);
				thread* load_ms2 = new thread([ms2_exp,&project,add_exp_num, win_num]() {
					ms2_exp->set_exp(project.get_ms2_mzML_filename(add_exp_num, win_num));
				});

				for (int cur_ref_exp_num = 0; cur_ref_exp_num < project_size; cur_ref_exp_num++) {
//					cout << cur_ref_exp_num << "\t" << add_exp_num << endl;
					swath.add_ms1(cur_ref_exp_num, add_exp_num, *ms1_exp,
							project.get_rt_normalizer_filename(cur_ref_exp_num, add_exp_num));
				}

				//Delete cur_exp
				delete ms1_exp;
				load_ms2->join();
				if (add_exp_num < project_size - 1) {
					//Load next experiment
					temp_ms1_exp = new SwathExperiment(ms1_para);
					load_ms1 = new thread([temp_ms1_exp,&project,add_exp_num, win_num]() {
						temp_ms1_exp->set_exp(project.get_ms1_mzML_filename(add_exp_num+1));
					});
				}

				swath.add_ms2(add_exp_num, *ms2_exp, project.get_feature_chro_filename(win_num),
						project.get_ms2_splited_chro_filename(run_num, add_exp_num, win_num));

				delete load_ms2;
				delete ms2_exp;
			}
			File::remove(project.get_feature_chro_filename(win_num));

//******************************* Start select ms2 *******************************
			cout << "Start select product ions" << endl;
			vector<MSSpectrum<Peak1D> > target_result, decoy_result;
			vector<string> ms2_splited_chro_filename;
			for (int cur_exp_num = 0; cur_exp_num < project_size; cur_exp_num++)
				ms2_splited_chro_filename.push_back(
						project.get_ms2_splited_chro_filename(run_num, cur_exp_num, win_num));

			swath.select_ms2(ms2_splited_chro_filename,
					project.get_ms2_chro_filename(run_num, win_num), target_result, decoy_result);

			for (int cur_exp_num = 0; cur_exp_num < project_size; cur_exp_num++)
				File::remove(project.get_ms2_splited_chro_filename(run_num, cur_exp_num, win_num));

			result.insert(result.end(), target_result.begin(), target_result.end());
			result.insert(result.end(), decoy_result.begin(), decoy_result.end());

			split_num++;

			cout << ref_exp_num * 100.0 / (project_size * 1.0) << "% done" << endl;
		}
		//Output result
		MSExperiment<Peak1D> result_exp;
		result_exp.setSpectra(result);

		GroupDIA::MzMLFile mzml;
		mzml.store(project.get_mzML_filename(run_num, win_num), result_exp);

		GroupDIA::MGFFile mgf;
		mgf.store_simple(project.get_mgf_filename(run_num, win_num), project.get_project_name(),
				result_exp);

		return EXECUTION_OK;
	}

}
;
int main(int argc, const char ** argv) {
	Analysis_SWATH tool;
	return tool.main(argc, argv);
}

/// @endcond
