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
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "IntensityCalc.h"

namespace GroupDIA {
void IntensityCalc::get_intensity(const vector<vector<double> >& product_intensity, int peak_start, int peak_length,
		vector<double>& calc_intensity) {
	if (product_intensity.empty()) {
		calc_intensity.clear();
		cout << "Intensity is empty, may have bugs.!" << endl;
		return;
	}

	int product_size = product_intensity.size();
	calc_intensity = vector<double>(product_size, 0);

	if (product_intensity[0].size() < 2) {
		cout << "Intensity is less than 2, may have bugs.!" << endl;
		return;
	}

	int intensity_size = product_intensity[0].size();
	for (auto it = product_intensity.begin(); it != product_intensity.end(); ++it)
		if (it->size() != intensity_size) {
			cout << "Intensity is not equal, may have bugs.!" << endl;
			return;
		}

	peak_start--;
	peak_length++;
	if (peak_start < 0) {
		peak_length += peak_start;
		peak_start = 0;
	}
	if (peak_length < 0)
		peak_length = 0;
	if (intensity_size - 2 < peak_start)
		peak_start = intensity_size - 2;
	if (peak_start + peak_length > intensity_size - 2)
		peak_length = intensity_size - peak_start - 2;

	if (peak_length == 0) {
		return;
	}
	int peak_end = peak_start + peak_length;

	//Construct RichPeakChromatogram
	typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram;

	vector<RichPeakChromatogram> chromatograms, picked_chroms, smoothed_chroms;
	for (int product_num = 0; product_num < product_size; product_num++) {
		RichPeakChromatogram chromatogram, smoothed_chrom, picked_chrom;
		chromatogram.resize(product_intensity[product_num].size());
		for (int i = 0; i < product_intensity[product_num].size(); i++) {
			chromatogram[i].setIntensity(product_intensity[product_num][i]);
			chromatogram[i].setMZ(i);
		}

		//Find by openswath
		pickChromatogram(chromatogram, smoothed_chrom, picked_chrom);

		chromatograms.push_back(chromatogram);
		picked_chroms.push_back(picked_chrom);
		smoothed_chroms.push_back(smoothed_chrom);
	}

	//Calc by openswath

	// Remove other, overlapping, picked peaks (in this and other
	// chromatograms) and then ensure that at least one peak is set to zero
	// (the currently best peak).
	remove_overlapping_features(picked_chroms, peak_start, peak_end);

	// Prepare linear resampling of all the chromatograms, here creating the
	// empty master_peak_container with the same RT (m/z) values as the reference
	// chromatogram.
	RichPeakChromatogram master_peak_container;
	master_peak_container.resize(product_intensity[0].size());
	int i = 0;
	for (auto it = master_peak_container.begin(); it != master_peak_container.end(); ++it) {
		it->setMZ(i);
		i++;
	}

	double total_intensity = 0;
	double total_peak_apices = 0;
	double total_xic = 0;
	for (Size k = 0; k < chromatograms.size(); k++) {
		const RichPeakChromatogram& chromatogram = chromatograms[k];
		for (auto it = chromatogram.begin(); it != chromatogram.end(); it++) {
			total_xic += it->getIntensity();
		}

		// resample the current chromatogram
		const RichPeakChromatogram used_chromatogram = resampleChromatogram_(chromatogram, master_peak_container, peak_start, peak_end);
		// const SpectrumT& used_chromatogram = chromatogram; // instead of resampling

		DoubleReal intensity_sum(0.0), rt_sum(0.0);
		// FEATURE : use RTBegin / MZBegin -> for this we need to know whether the template param is a real chromatogram or a spectrum!
		for (auto it = used_chromatogram.begin(); it != used_chromatogram.end(); it++) {
			if (it->getMZ() > peak_start && it->getMZ() < peak_end) {
				intensity_sum += it->getIntensity();
			}
		}

		//Remove background
		if (_background_subtraction_ != "none") {
			double background = 0;
			// we use the smoothed chromatogram here to have a more accurate estimatation of the noise at the flanks of the peak
			if (_background_subtraction_ == "smoothed") {
				if (smoothed_chroms.size() <= k) {
					background = 0;
				} else {
					background = calculateBgEstimation_(smoothed_chroms[k], peak_start, peak_end);
				}
			} else if (_background_subtraction_ == "original") {
				background = calculateBgEstimation_(used_chromatogram, peak_start, peak_end);
			}
			intensity_sum -= background;
			if (intensity_sum < 0)
				intensity_sum = 0;
		}

		calc_intensity[k] = intensity_sum;
	}
}
}/* namespace GroupDIA */
