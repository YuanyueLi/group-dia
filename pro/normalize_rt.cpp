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
#include <ToolsSpectrum.h>
#include <RTNormalizer.h>
#include <ToolsFilename.h>
#include <stdio.h>
using namespace std;
using namespace OpenMS;
using namespace GroupDIA;
class NormalizeRT: public TOPPBase {
public:
	NormalizeRT() :
			TOPPBase("Normalize RT", "Normalize the retention time", false) {
	}

protected:
	void registerOptionsAndFlags_() {
		registerInputFile_("in", "<file>", "", "Input parameter list");
		setValidFormats_("in", StringList::create("txt"));

		registerSubsection_("algorithm", "Algorithm parameters section");
	}

	Param getSubsectionDefaults_(const String &) const {
		Param defaults_;
		defaults_.setValue("mz_start", 400.0, "");
		defaults_.setValue("mz_end", 1200.0, "");
		defaults_.setValue("mz_space", 0.1, "");
		return defaults_;
	}

	ExitCodes main_(int, const char **) {
		String in = getStringOption_("in");
		fstream fi(in, fstream::in);

		Param para = getParam_().copy("algorithm:", true);

		ToolsFilename project;
		int exp_size = project.read_filelist(in);

		vector<Ms1Intensity> ms1_exp(exp_size);

		cout << "Start load ms1" << endl;
		for (int i = 0; i < exp_size; i++) {
			ms1_exp[i].load_mzml(para, project.get_ms1_mzML_filename(i));
		}

		cout << "Start align" << endl;
		//Normalize RT
#pragma omp parallel for schedule(dynamic) collapse(2)
		for (int i = 0; i < exp_size; i++) {
			for (int j = 0; j < exp_size; j++) {
				RTNormalizer rt;
				rt.align(ms1_exp[i], ms1_exp[j]);
				rt.store(project.get_rt_normalizer_filename(i, j));
			}
		}

		return EXECUTION_OK;
	}

};
int main(int argc, const char ** argv) {
	NormalizeRT tool;
	return tool.main(argc, argv);
}

/// @endcond
