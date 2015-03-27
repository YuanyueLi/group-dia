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

#ifndef TRANSITIONSWATHSCORE_H_
#define TRANSITIONSWATHSCORE_H_
#include "ToolsScore.h"
#include "TransitionIdentified.h"
#include "ToolsTabFile.h"

namespace GroupDIA {

class TransitionSwathScore {
public:
	TransitionSwathScore(TransitionIdentified& trans);
	void generate_score_tabfile(ToolsTabFile& score);
	void generate_intensity_tabfile(ToolsTabFile& intensity);
	void get_transition_score(map<string, double>& map_score);
	void get_transition_info(map<string, string>& map_info);

private:
	TransitionIdentified& trans;

	vector<string> info_name;
	vector<string> info_value;
	map<string, double> map_score;

	void inline add_score(string name, double value);

	void cal_distance_score();
	void cal_max_intensity_score(vector<double>& max_intensity);
	void cal_sn_intensity_score(vector<double>& exp_sum_intensity);
	void cal_similar_score();

	int exp_size;
	int product_size;
	int exist_exp_size;
//	double factor_exp_size;
	vector<double> exp_weight;
	vector<double> ion_weight;
};

} /* namespace GroupDIA */
#endif /* TRANSITIONSWATHSCORE_H_ */
