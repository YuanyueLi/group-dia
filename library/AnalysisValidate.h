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


#ifndef ANALYSISVALIDATE_H_
#define ANALYSISVALIDATE_H_
#include "ToolsSpectrum.h"
#include "ToolsTabFile.h"
#include "TransitionIdentified.h"
#include "ToolsLoadSave.h"
#include "TransitionSwathScore.h"
#include "TransitionMProhetResult.h"
#include "ToolsTabFileWrite.h"
#include <ToolsSpectrum.h>
#include <GroupProphet.h>
#include <boost/math/distributions/normal.hpp>

namespace GroupDIA {

struct TransitionInfo {
	int decoy_type;
	string accessions;
	int charge;
	string sequence;
	double mz;
	int scan_num;
};

class AnalysisValidate {
public:
	AnalysisValidate(const Param& para);
	void merge_gprophet_input_file(const vector<string>& input_filename,
			string output_filename) const;

	void prepare_for_pyprophet(const vector<string>& input_filename,
			const vector<string>& output_filename) const;
	void load(vector<string> filename);
	void load(string filename);
	bool run_mprohet();
	void output_mprohet_result(string input_filename, string output_filename);

	void load_mprohet_result(int project_size, string filename);
	double get_gprophet_min_tran_p(string filename) const;
	double get_gprophet_min_spec_p(string filename, double tran_p_cutoff) const;
	void score(int exp_size, string score_filename);
	void calc_fdr();
	void output_intensity(string filename);

private:
	GroupProphet mprophet;
	double ACCEPT_FDR;

	vector<TransitionMProhetResult> trans_score;
	vector<TransitionInfo> trans_info;
	map<int, int> id_to_no;
	ToolsTabFile result_intensity;

};

} /* namespace GroupDIA */

#endif /* ANALYSISVALIDATE_H_ */
