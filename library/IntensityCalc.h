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


#ifndef INTENSITYCALC_H_
#define INTENSITYCALC_H_
#include "ToolsSpectrum.h"
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>
namespace GroupDIA {
const string _background_subtraction_ = "smoothed";
class IntensityCalc: private MRMTransitionGroupPicker {
public:
	void get_intensity(const vector<vector<double> >& intensity, int peak_start, int peak_length, vector<double>& result_intensity);

};

} /* namespace GroupDIA */
#endif /* INTENSITYCALC_H_ */
