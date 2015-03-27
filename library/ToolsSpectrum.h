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
#ifndef TOOLSSPECTRUM_H_
#define TOOLSSPECTRUM_H_
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
//#include <boost/numeric/conversion/cast.hpp>
#include <random>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <limits>
#include <sys/time.h>
#include <math.h>
#include <tuple>
#include <fstream>
#include <cmath>
#include <iostream>
#include <omp.h>

using namespace std;
using namespace OpenMS;
typedef OpenMS::MSExperiment<Peak1D> MapType;
typedef MSSpectrum<Peak1D> SpectrumT;
typedef MSSpectrum<RichPeak1D> RichPeakSpectrum;
typedef pair<int, int> PeakInfo;

namespace GroupDIA {
const double MS2_WIN = 0.05;

class ToolsSpectrum {
public:
	ToolsSpectrum();
	virtual ~ToolsSpectrum();

	static void mean_and_std_dev(const vector<double>& values, double& mean, double& sd);
	static void median_and_sad(vector<double> values, double& median, double& sad);

	static bool is_target_identification_result(const vector<String>& acc);

	static bool is_same_in_ppm(const double a, const double b, const double ppm);

	template<typename IteratorType1, typename IteratorType2>
	static DoubleReal cal_peak_apex(const IteratorType1 intensity_begin,
			const IteratorType1 intensity_end, const IteratorType2 rt_begin,
			const IteratorType2 rt_end) {
		double sum = 0.0;
		auto i = intensity_begin;
		auto r = rt_begin;
		for (; i != intensity_end and r != rt_end; ++i, ++r) {
			sum += *i * *r;
		}
		double sum_intensity = accumulate(intensity_begin, intensity_end, 0.0);
		if (sum_intensity > 0)
			return sum / sum_intensity;
		else
			return 0;
	}

	template<typename IteratorType1>
	static int calc_max_peak(const IteratorType1 intensity_begin,
			const IteratorType1 intensity_end) {
		auto imax = max_element(intensity_begin, intensity_end);
		return distance(intensity_begin, imax);
	}

	static int get_rt_no(const vector<double>& rt, double r);
	static double get_rt(const vector<double>& rt_no, int r);

	static void standardize_data(std::vector<double> & data);

	//Cal cross correlation
	static void get_best_delay(vector<double> a, vector<double> b, int max_delay,
			vector<pair<int, double> >& dela_score);

	static void normalize_sum(vector<double>& x);
	static void normalize_max(vector<double>& x);

	static void convert_aasequence_to_peptide(const AASequence& aa_sequence,
			TargetedExperiment::Peptide& result_peptide);

	template<typename T> static T get_nth_small(vector<T>& input, int n) {
		if (input.empty())
			throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, "");
		if (n > input.size())
			n = input.size();
		nth_element(input.begin(), input.begin() + n - 1, input.end());
		return input[n - 1];
	}

	template<typename T> static T get_nth_large(vector<T>& input, int n) {
		if (input.empty())
			throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, "");
		if (n > input.size())
			n = input.size();
		nth_element(input.begin(), input.begin() + n - 1, input.end(), greater<T>());
		return input[n - 1];
	}

	static void get_pearson_score(const vector<double>& precursor_intensity,
			const vector<vector<double> >& product_intensity, vector<double>&score);

	static void peak_finder(const vector<double>& intensity, int peak_finder_start,
			int peak_finder_length, vector<PeakInfo>& peak_info);

	template<typename IteratorType> static void random(IteratorType v_begin, IteratorType v_end) {
		if (v_begin == v_end)
			return;
		int size = distance(v_begin, v_end);
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dis(0, size - 1);
		auto generator = bind(dis, gen);
		random_shuffle(v_begin, v_end, generator);
	}

private:
	static void add_modification_(std::vector<TargetedExperiment::Peptide::Modification>& mods,
			int location, ResidueModification& rmod, const String& name);
};
} /* namespace GroupDIA */
#endif /* TOOLSSPECTRUM_H_ */
