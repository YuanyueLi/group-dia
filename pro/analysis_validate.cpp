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
#include <AnalysisValidate.h>
#include <ToolsFilename.h>
#include <stdlib.h>
#define DEBUG
using namespace std;
using namespace OpenMS;
using namespace GroupDIA;
class Analysis_Validate: public TOPPBase {
public:
	Analysis_Validate() :
			TOPPBase("Analysis Validate", "Generate the final result.", false) {
	}
protected:
	void registerOptionsAndFlags_() {
		registerInputFile_("in", "<file>", "", "Input parameter list");
		setValidFormats_("in", StringList::create("txt"));

		registerSubsection_("algorithm", "Algorithm parameters section");
	}

	Param getSubsectionDefaults_(const String &) const {
		Param defaults_;
		defaults_.setValue("max_accepted_fdr", 0.01, "");
		defaults_.insert("mProhet:", GroupProphet().getParameters());
		return defaults_;
	}

	typedef OpenMS::MSExperiment<Peak1D> MapType;
	ExitCodes main_(int, const char **) {
		typedef OpenMS::MSExperiment<Peak1D> MapType;
		String in = getStringOption_("in");
		String threads_num = getParam_().getValue("threads").toString();
		Param para = getParam_().copy("algorithm:", true);

		ToolsFilename project;
		int project_size = project.read_filelist(in);
		string out = project.get_project_name();
		int windows_size = project.get_swath_windows_size();

		AnalysisValidate analysis(para);
		vector<string> input_filenames, output_filenames;
		for (int i = 0; i < windows_size; i++) {
			String cur_filename = project.get_identification_result(i);
			if (File::exists(cur_filename)) {
				input_filenames.push_back(cur_filename);
			} else {
				cout << "The " << i << " file is not exist" << endl;
				input_filenames.push_back("");
			}
		}
		for (int i = 0; i < project_size; i++) {
			output_filenames.push_back(project.get_mprohet_input_filename(i));
		}

		cout << "Start merge GroupProphet input files" << endl;
		analysis.merge_gprophet_input_file(input_filenames, output_filenames[0]);
		cout << "Start load GroupProphet input file" << endl;
		analysis.load(output_filenames[0]);
		analysis.run_mprohet();
		analysis.output_mprohet_result(output_filenames[0], project.get_mprohet_output_filename(0));

//		for (int i = 0; i < windows_size; i++) {
//			if (!File::exists(project.get_identification_result(i)))
//				continue;
//			AnalysisValidate cur_analysis(para);
//			cur_analysis.load(project.get_identification_result(i));
//			bool is_existed = cur_analysis.run_mprohet();
//			if (is_existed)
//				cur_analysis.output_mprohet_result(project.get_identification_result(i),
//						project.get_mprohet_output_filename(i));
//			else {
//				cout << "Error in GroupProphet" << endl;
//			}
//		}

		cout << "Start score transition" << endl;
//		for (int i = 0; i < windows_size; i++) {
//			analysis.load_mprohet_result(project_size, project.get_mprohet_output_filename(i));
//		}
		analysis.load_mprohet_result(project_size, project.get_mprohet_output_filename(0));
		analysis.score(project_size, project.get_mprohet_2nd_input_filename());
//		analysis.calc_fdr();
		analysis.output_intensity(out + ".intensity.csv");

		return EXECUTION_OK;
	}
};

int main(int argc, const char ** argv) {
	Analysis_Validate tool;
	return tool.main(argc, argv);
}

/// @endcond
