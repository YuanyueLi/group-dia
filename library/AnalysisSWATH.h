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


#ifndef ANALYSISSWATH_H_
#define ANALYSISSWATH_H_
#include "TransitionGroup.h"
#include "AnalysisSimilar.h"
#include "ToolsSpectrum.h"

namespace GroupDIA {

class AnalysisSWATH {
public:
	AnalysisSWATH(const Param& para);

	int size() const {
		return trans.size();
	}

	bool set_ms1(const Feature& feature, int exp_num, int scan_num, const SwathExperiment& exp);

	void add_ms1(int ref_exp_num, int add_exp_num, const SwathExperiment& exp, string rt_filename);

	void set_ms2(int exp_num, SwathExperiment& ms2_exp, double ms1_win_start, double ms1_win_end,
			string output_chro_filename);

	void refine();

	void add_ms2(int exp_num, SwathExperiment& ms2_exp, string input_chro_filename,
			string output_chro_filename);

	void select_ms2(vector<string> intput_chro_filename, string out_chro_filename,
			vector<MSSpectrum<Peak1D> >& target_spectrum,
			vector<MSSpectrum<Peak1D> >& decoy_spectrum);
protected:
	vector<TransitionGroup*> trans;

	Param para;

	int MIN_ALLOWED_PRODUCT_NUM;
	int RT_WINDOWS;

	double MIN_ALLOWED_CORRELATION;
	int MAX_ALLOWED_DELAY;
	int FINE_NEARBY;
	double MIN_INTENSITY;

//	double MIN_ALLOWED_CORRELATION = 0.3;
//	int MAX_ALLOWED_DELAY = 60;
//	int FINE_NEARBY = 10;
//	double MIN_INTENSITY = 300;
//	const int MIN_ALLOWED_PRODUCT_NUM = 4;
//	const int RT_WINDOWS = 100;
};

} /* namespace GroupDIA */

#endif /* ANALYSISSWATH_H_ */
