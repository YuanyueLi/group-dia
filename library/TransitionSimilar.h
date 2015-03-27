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

#ifndef TRANSITIONSIMILAR_H_
#define TRANSITIONSIMILAR_H_

#include "TransitionGroup.h"

namespace GroupDIA {
typedef vector<double>::iterator VecIt;

class TransitionSimilar: public TransitionGroup {
public:
	void crude_refine(double MIN_ALLOWED_CORRELATION);

	void refined_by_crosscorrelation(int MAX_ALLOWED_DELAY, int FINE_NEARBY, double MIN_INTENSITY);

private:
	//ref_intensity_end - ref_intensity_start should be the same with sam_intensity_end-sam_intensity_start
	void determine_peak(const VecIt ref_intensity_start, const VecIt ref_intensity_end,
			const VecIt sam_intensity_start, const VecIt sam_intensity_end, int ref_peak_start,
			int ref_peak_length, int& sam_peak_start, int& sam_peak_length) const;

//	const double MIN_ALLOWED_CORRELATION = 0.3;
//	const int MAX_ALLOWED_DELAY = 60;
//	const int FINE_NEARBY = 10;
//	const double MIN_INTENSITY = 300;
};

} /* namespace GroupDIA */

#endif /* TRANSITIONSIMILAR_H_ */
