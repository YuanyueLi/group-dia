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

#ifndef RTNORMALIZER_H_
#define RTNORMALIZER_H_
#include <OpenMS/FILTERING/SMOOTHING/LowessSmoothing.h>
#include "ToolsSpectrum.h"
#include "ToolsLoadSave.h"
namespace GroupDIA {

class Ms1Intensity {
	friend class RTNormalizer;
public:
	void load_mzml(const Param& para, string ms1_mzml_filename);

private:
	int size() const;
	double get_origin_time(int i) const;

	vector<double> origin_time;
	vector<vector<double> > intensity;
	vector<double> sum;
};

class RTNormalizer {
public:
	vector<double> align(const Ms1Intensity& exp_ref, const Ms1Intensity& exp_sam);
	double get_sample_rt(double reference_rt) const;

	void store(string rt_filename) const;
	void load(string rt_filename);

private:
	vector<double> find_anchor(const Ms1Intensity& exp_ref, const Ms1Intensity& exp_sam,
			vector<vector<int> >& path);
	vector<double> get_best_point(const Ms1Intensity& exp_ref, const Ms1Intensity& exp_sam,
			vector<vector<int> >& path);
	double get_dist(const Ms1Intensity& exp_ref, int time_ref, const Ms1Intensity& exp_sam,
			int time_sam);

	map<double, double> anchor; //reference rt,smaple rt
	vector<vector<double> > dist;
};

} /* namespace GroupDIA */
#endif /* RTNORMALIZER_H_ */
