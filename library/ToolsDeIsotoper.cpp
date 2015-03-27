/*
 * ToolsDeIsotoper.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: lyy
 */

#include "ToolsDeIsotoper.h"

const int MAX_ISOTOPE_NUM = 8;
const int MAX_CHARGE = 5, MIN_CHARGE = 1;
const double MIN_ISOTOPE_CORRELATION = 0.9;
namespace GroupDIA {

ToolsDeIsotoper::ToolsDeIsotoper(Param para) {
	use_ppm = para.getValue("use_ppm").toBool();
	if (use_ppm)
		delta_ppm = para.getValue("delta_ppm");
	else
		delta_mz = para.getValue("delta_mz");
}

ToolsDeIsotoper::ToolsDeIsotoper() {
	use_ppm = false;
	delta_mz = 0.05;
}

void ToolsDeIsotoper::de_isotoper(const MSSpectrum<Peak1D>& ori_spec, MSSpectrum<Peak1D>& result_spec) const {
	//////////////////////////////////////////////////////////
	//Get matched peaks
	//////////////////////////////////////////////////////////
	//Generate origin peaks
	vector<double> mz_removed;

	for (int charge = MAX_CHARGE; charge >= MIN_CHARGE; --charge) {
		MSSpectrum<Peak1D>::const_iterator i_ori_spec_begin = ori_spec.begin();
		MSSpectrum<Peak1D>::const_iterator i_ori_spec_end = ori_spec.end();
		vector<Peak1D> peak_spec;
		while (find_peak_group(i_ori_spec_begin, i_ori_spec_end, charge, i_ori_spec_begin, peak_spec)) {
			//Get distribution of isotope
			IsotopeDistribution dis;
			dis.setMaxIsotope(peak_spec.size());
			dis.estimateFromPeptideWeight(peak_spec[0].getMZ() * charge - OpenMS::Constants::C13C12_MASSDIFF_U);

			vector<double> isotope_intensity;
			isotope_intensity.reserve(peak_spec.size() + 1);
			for (auto it = dis.begin(); it != dis.end(); ++it)
				isotope_intensity.push_back(it->second);
			isotope_intensity.push_back(0);

			//Calc the correlation
			vector<double> real_intensity;
			real_intensity.reserve(peak_spec.size() + 1);
			for (auto it = peak_spec.begin(); it != peak_spec.end(); ++it)
				real_intensity.push_back(it->getIntensity());
			real_intensity.push_back(0);

			double cor = OpenMS::Math::pearsonCorrelationCoefficient(isotope_intensity.begin(), isotope_intensity.end(),
					real_intensity.begin(), real_intensity.end());

			//If the correlation is larger than thresthold, record the peak.
			if (cor >= MIN_ISOTOPE_CORRELATION) {
				auto ip = peak_spec.begin();
				++ip;
				for (; ip != peak_spec.end(); ++ip) {
					mz_removed.push_back(ip->getMZ());
				}
			}
		}
	}

//Generate output spectrum
	result_spec = ori_spec;
	result_spec.clear(false);
	MSSpectrum<Peak1D> temp_spec = ori_spec;
	temp_spec.sortByPosition();
	sort(mz_removed.begin(), mz_removed.end());

	auto iwalker = temp_spec.begin();
	auto iremove = mz_removed.begin();
	while (iwalker != temp_spec.end() and iremove != mz_removed.end()) {
		if (iwalker->getMZ() < *iremove) {
			result_spec.push_back(*iwalker);
			++iwalker;
		} else if (iwalker->getMZ() > *iremove) {
			++iremove;
		} else {
			++iwalker;
			++iremove;
		}
	}
	while (iwalker != temp_spec.end()) {
		result_spec.push_back(*iwalker);
		++iwalker;
	}
}

bool ToolsDeIsotoper::find_peak_group(const MSSpectrum<Peak1D>::const_iterator& begin, const MSSpectrum<Peak1D>::const_iterator& end,
		int charge, MSSpectrum<Peak1D>::const_iterator& next_begin, vector<Peak1D>& result) const {
	next_begin = begin;
	while (next_begin != end) {
		result.clear();

		int num = 1;
		bool exist = true;
		auto iwalker = next_begin;
		while (exist) {
			//Add the peak
			result.push_back(*iwalker);
			if (result.size() >= MAX_ISOTOPE_NUM)
				break;

			double wanted_mz = next_begin->getMZ() + OpenMS::Constants::C13C12_MASSDIFF_U * num / charge;

			double eps;
			if (use_ppm)
				eps = delta_ppm * iwalker->getMZ() / 1.0e6;
			else
				eps = delta_mz;
			double mzl = wanted_mz - eps;
			double mzr = wanted_mz + eps;

			//Move iwalker to the top peak in [mzl,mzr]
			while (iwalker->getMZ() <= mzl and iwalker != end)
				++iwalker;
			if (iwalker == end or iwalker->getMZ() > mzr)
				exist = false;
			else {
				auto itemp = iwalker;
				auto imax = iwalker;
				while (itemp != end and itemp->getMZ() <= mzr) {
					if (itemp->getIntensity() > imax->getIntensity())
						imax = itemp;
					++itemp;
				}
				iwalker = imax;
			}
			num++;
		}
		++next_begin;
		if (result.size() > 1) {
			for (auto it = result.begin(); it != result.end(); ++it) {
//				cout << it->getMass() << ":" << it->getIntensity() << endl;
			}
			return true;
		}
	}
	return false;
}
} /* namespace GroupDIA */
