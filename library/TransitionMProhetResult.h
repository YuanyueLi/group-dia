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

#ifndef TRANSITIONMPROHETRESULT_H_
#define TRANSITIONMPROHETRESULT_H_
#include <climits>
#include <boost/math/distributions/students_t.hpp>
#include "ToolsSpectrum.h"
namespace GroupDIA {

struct MProhetResult {
	double left_peak;
	double right_peak;
	double p;
	double ide_score;
	vector<double> intensity;

	void add_intensity(const String& ori_value);
};

class TransitionMProhetResult {
public:
	void remove_interference_ions();

	void add_result(int exp_num, MProhetResult& result);

	int size() const;

	double get_intensity(int exp_num) const;

	double get_average_d_score() const;

	double get_max_p() const;

	double get_ide_score() const;

	void get_experiment_dist_score(map<string, double>& score) const;

	void get_distance_score(map<string, double>& score) const;

private:
	map<int, MProhetResult> tran;
};

} /* namespace GroupDIA */

#endif /* TRANSITIONMPROHETRESULT_H_ */
