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

#ifndef ANALYSISIDENTIFICATION_H_
#define ANALYSISIDENTIFICATION_H_
#include "ToolsScore.h"
#include "ToolsTabFile.h"
#include "TransitionIdentified.h"
#include "TransitionSwathScore.h"
#include "ScanPepXMLFile.h"
#include "TransitionScore.h"
#include "TestTransition.h"
#include "MergeTransitions.h"
#include <OpenMS/FORMAT/MzIdentMLFile.h>
//#define DEBUG
namespace GroupDIA {

class AnalysisIdentification {
public:
	AnalysisIdentification(const Param& _para, int windows_size, int cur_windows_num,int exp_size);

	void reexact_intensity(string pepxml_filename, string input_chro_filename, string temp_chro_filename,
			vector<String> ms1_mzML_filename, vector<String> ms2_chro_filename,
			vector<String> ms2_mzML_filename, string output_chro_filename);

	void score_transition(string chro_filename, string temp_chro_filename,
			vector<String> ms2_chro_filename, vector<String> ms2_mzML_filename,
			string out_chro_filename, string output_score_filename);

protected:
	Param para;
	Param para_ide;

	int exp_size;
	int windows_size;
	int cur_windows_num;

	int MIN_PRODUCT_IONS_NUM;

	multimap<int, PeptideHit> map_scan_pep;
};

} /* namespace GroupDIA */
#endif /* ANALYSISIDENTIFICATION_H_ */
