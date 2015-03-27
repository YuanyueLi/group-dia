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

#ifndef TOOLSSPECTRUMANALYSIS_H_
#define TOOLSSPECTRUMANALYSIS_H_

#include "ToolsSpectrum.h"
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/StatsHelpers.h>

const int MAX_PEAK_SPACE = 3;
namespace GroupDIA {
class ToolsSwathGaussFilter: public OpenMS::GaussFilter {
public:
	template<typename PeakType>
	void swath_filter(const MSSpectrum<PeakType> & ori_spectrum,
			MSSpectrum<PeakType> & result_spectrum) {
		double min_space = this->spacing_;
		for (auto it = ori_spectrum.begin(); it != ori_spectrum.end(); ++it) {
			if (it->getIntensity() < min_space)
				min_space = it->getIntensity();
		}
		this->spacing_ = min_space;

		double gauss_width_ = 0.1;
		Param filter_parameters = getParameters();
		filter_parameters.setValue("gaussian_width", gauss_width_);
		setParameters(filter_parameters);
		result_spectrum = ori_spectrum;
		GaussFilter::filter(result_spectrum);
	}

};

class ToolsSwathPickPeak {
public:
	template<typename PeakType>
	void pick_peak(const MSSpectrum<PeakType>& input, MSSpectrum<PeakType>& output) const {
		PeakPickerHiRes pp;
		Param pepi_param = PeakPickerHiRes().getDefaults();
		pp.setParameters(pepi_param);
		pp.pick(input, output);
	}

	template<typename PeakType>
	void refine_peak(const MSSpectrum<PeakType>& input, const Param& para,
			MSSpectrum<PeakType>& output) const {
		//Get para
		bool use_ppm = para.getValue("use_ppm").toBool();
		double delta_ppm, delta_mz;
		if (use_ppm)
			delta_ppm = para.getValue("delta_ppm");
		else
			delta_mz = para.getValue("delta_mz");

		//Start refine
		output = input;
		output.clear(false);
		if (input.begin() == input.end())
			return;

		auto icur = input.begin();
		auto inext = icur + 1;
		for (; inext != input.end(); ++icur, ++inext) {
			double delta;
			if (use_ppm)
				delta = icur->getMZ() * delta_ppm / 1.0e6;
			else
				delta = delta_mz;

			//the space is larger than resoltion
			if (inext->getMZ() - icur->getMZ() >= delta) {
				output.push_back(*icur);
				continue;
			}

			//the space is smaller than resoltion
			auto istart = icur;
			for (; inext != input.end(); ++icur, ++inext) {
				if (inext->getMZ() - icur->getMZ() >= delta)
					break;
			}
			auto iend = inext;
			double c = 0, total_intensity = 0, max_intensity = 0;
			for (auto it = istart; it != iend; ++it) {
				c += it->getIntensity() * it->getMZ();
				total_intensity += it->getIntensity();
				if (it->getIntensity() > max_intensity)
					max_intensity = it->getIntensity();
			}
			PeakType new_peak;
			new_peak.setMZ(c / total_intensity);
			new_peak.setIntensity(max_intensity);
			output.push_back(new_peak);
			if (inext == input.end())
				break;
		}
		output.push_back(*icur);
	}

	template<typename PeakType>
	void pick_peak_apex(const MSSpectrum<PeakType>& input, MSSpectrum<PeakType>& output,
			double min_space = 0) const {
		output = input;
		output.clear(false);

		if (input.size() <= 3)
			return;
		if (min_space == 0) {
			min_space = 100000;
			for (auto it = input.begin(); it != input.end(); ++it) {
				if (it->getIntensity() < min_space)
					min_space = it->getIntensity();
			}
		}

		auto ipre = input.begin();
		auto ic = ipre + 1;
		auto inext = ic + 1;
		const auto ie = input.end();
		int no_space_num = 0;
		while (inext != ie) {
			//The peak is the highest.
			int pre_space_num = (ic->getMZ() - ipre->getMZ()) / min_space + 0.5;
			int next_space_num = (inext->getMZ() - ic->getMZ()) / min_space + 0.5;
			if ((ic->getIntensity() > ipre->getIntensity() or pre_space_num > MAX_PEAK_SPACE)
					and (ic->getIntensity() > inext->getIntensity()
							or next_space_num > MAX_PEAK_SPACE)) {
				Peak1D p;
				p.setIntensity(ic->getIntensity());
				p.setMZ(ic->getMZ());
				output.push_back(p);
			}

			++ipre, ++ic, ++inext;
		}
	}

};

} /* namespace GroupDIA */
#endif /* TOOLSSPECTRUMANALYSIS_H_ */
