/*
 * ToolsDeIsotoper.h
 *
 *  Created on: Mar 27, 2015
 *      Author: lyy
 */

#ifndef TOOLSDEISOTOPER_H_
#define TOOLSDEISOTOPER_H_
#include "ToolsSpectrum.h"
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/StatsHelpers.h>
namespace GroupDIA {

class ToolsDeIsotoper {
public:
	ToolsDeIsotoper(Param para);
	ToolsDeIsotoper();
	void de_isotoper(const MSSpectrum<Peak1D>& ori_spec, MSSpectrum<Peak1D>& result_spec) const;

private:
	bool find_peak_group(const MSSpectrum<Peak1D>::const_iterator& begin,
			const MSSpectrum<Peak1D>::const_iterator& end, int charge,
			MSSpectrum<Peak1D>::const_iterator& next_begin, vector<Peak1D>& result) const;

	bool use_ppm;
	double delta_ppm;
	double delta_mz;
};

} /* namespace GroupDIA */

#endif /* TOOLSDEISOTOPER_H_ */
