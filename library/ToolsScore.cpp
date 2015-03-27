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

#include "ToolsScore.h"
namespace GroupDIA {
double ToolsScore::cal_sn_score(const vector<double>& intensity, const vector<double>&rt, double apex_rt) {
	OpenMS::MSSpectrum<ChromatogramPeak> origin_chro;
	if (intensity.size() != rt.size()) {
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, String(intensity.size()) + " " + String(rt.size()));
	}
	for (auto it_i = intensity.begin(), it_rt = rt.begin(); it_i != intensity.end(); ++it_i, ++it_rt) {
		ChromatogramPeak p;
		p.setPos(*it_rt);
		p.setIntensity(*it_i);
		origin_chro.push_back(p);
	}

	DoubleReal sn_win_len_ = 1000.0;
	DoubleReal sn_bin_count_ = 30.0;
	OpenSwath::ISignalToNoisePtr snptr(new OpenMS::SignalToNoiseOpenMS<ChromatogramPeak>(origin_chro, sn_win_len_, sn_bin_count_));
	double sn = snptr->getValueAtRT(apex_rt);
	if (sn <= 0.001)
		sn = 0.001;
	sn = log10(sn);
	return sn;
}

} /* namespace GroupDIA */
