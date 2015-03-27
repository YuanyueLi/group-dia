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

#include "Precursor.h"

namespace GroupDIA {
PrecursorIon::PrecursorIon() {
	nearby_intensity.clear();
	nearby_rt.clear();
	peak_rt_start_no = 0, peak_rt_length = 0, peak_rt_start_in_nearby_no = 0;
}

bool PrecursorIon::set_ms1(const SwathExperiment& exp, double mz, double apex_rt,
		int exact_scan_num) {
	int apex_rt_no = exp.get_rt_no(apex_rt);
	if (apex_rt_no >= exp.size() - 1 || apex_rt_no <= 0) {
		return false;
	}
	if (apex_rt_no >= exp.size() - 1)
		apex_rt_no = exp.size() - 1;

	vector<double> peak_intensity;
	bool exist = exp.find_peak(mz, apex_rt_no, peak_rt_start_no, peak_rt_length, peak_intensity);
	if (!exist) {
		return false;
	}

	set_peak(exp, mz, peak_rt_start_no, peak_rt_length, exact_scan_num);
	if (peak_rt_length <= 0
			or peak_rt_start_in_nearby_no + peak_rt_length > nearby_intensity.size())
		return false;

	if (accumulate(nearby_intensity.begin() + peak_rt_start_in_nearby_no,
			nearby_intensity.begin() + peak_rt_start_in_nearby_no + peak_rt_length, 0.0) == 0)
		return false;
	return true;
}

void PrecursorIon::set_peak(const SwathExperiment& exp, double set_mz, int set_peak_rt_start_no,
		int set_peak_rt_length, int exact_scan_num) {
	peak_rt_start_no = set_peak_rt_start_no;
	peak_rt_length = set_peak_rt_length;
//Set peak_rt_start_no peak_rt_length
	if (set_peak_rt_start_no < 0) {
		peak_rt_start_no = 0;
	} else if (set_peak_rt_start_no >= exp.size()) {
		peak_rt_length = 0;
		peak_rt_start_no = exp.size() - 1;
	}

	if (peak_rt_start_no + peak_rt_length > exp.size()) {
		peak_rt_length = exp.size() - peak_rt_start_no;
	}

	exp.get_nearby_intensity(set_mz, peak_rt_start_no, peak_rt_length, exact_scan_num,
			nearby_intensity, peak_rt_start_in_nearby_no);

	int nearby_rt_length = nearby_intensity.size();
	exp.get_rt(get_nearby_rt_start_no(), nearby_rt_length, nearby_rt);
}

void PrecursorIon::add_ms1(const SwathExperiment& exp, double mz, const RTNormalizer& rt,
		int exact_scan_num, PrecursorIon& result) const {
	if (!this->is_exist()) {
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Add error");
	}

//Define this precursor ion as "ref", define wanted result as "sam"
	const auto& ref = *this;
	auto& sam = result;
	double ref_peak_apex_rt = ref.get_peak_apex_rt();
	double sam_peak_apex_rt = rt.get_sample_rt(ref_peak_apex_rt);

//#pragma omp critical
//	cout << ref_peak_apex_rt << "\t" << sam_peak_apex_rt << endl;
	sam.set_peak(exp, mz, exp.get_rt_no(sam_peak_apex_rt), 0, exact_scan_num);

	if (result.peak_rt_start_no + result.peak_rt_length > exp.size()) {
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Add error");
	}
	if (result.peak_rt_start_no - result.peak_rt_start_in_nearby_no + result.nearby_intensity.size()
			> exp.size())
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Add error");
}

void PrecursorIon::add_ms1(const SwathExperiment& exp, double mz, int real_rt_no, int rt_length,
		int exact_scan_num, PrecursorIon& result) const {
	if (!this->is_exist()) {
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Add error");
	}

	auto& sam = result;
	sam.set_peak(exp, mz, real_rt_no, rt_length, exact_scan_num);

	if (result.peak_rt_start_no + result.peak_rt_length > exp.size()) {
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Add error");
	}
	if (result.peak_rt_start_no - result.peak_rt_start_in_nearby_no + result.nearby_intensity.size()
			> exp.size())
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Add error");
}

void PrecursorIon::get_peak_rt(double& rt_start, double& rt_apex, double& rt_end) const {
	int rt_start_in_nearby_no = peak_rt_start_in_nearby_no;
	int length = peak_rt_length;

	rt_start = nearby_rt[rt_start_in_nearby_no];
	rt_end = nearby_rt[rt_start_in_nearby_no + length];
	rt_apex = get_peak_apex_rt();
}

bool PrecursorIon::is_exist() const {
	return peak_rt_length != 0;
}

void PrecursorIon::not_exist() {
	peak_rt_length = 0;
}

double PrecursorIon::get_peak_apex_rt() const {
	if (nearby_rt.empty())
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"The size of peak is empty.");

	return nearby_rt[get_peak_apex_rt_in_nearby_no()];
}

void PrecursorIon::get_peak_intensity(vector<double>& intensity) const {
	if (peak_rt_length == 0)
		intensity = vector<double>(1, 0);
	else
		intensity = vector<double>(nearby_intensity.begin() + peak_rt_start_in_nearby_no,
				nearby_intensity.begin() + peak_rt_start_in_nearby_no + peak_rt_length);
}

void PrecursorIon::get_peak_rt_no(int& rt_start_no, int& length) const {
	rt_start_no = peak_rt_start_no;
	length = peak_rt_length;
}

void PrecursorIon::load_store(fstream& f, int io_type) {
	LoadSave::io(f, peak_rt_start_no, io_type);
	LoadSave::io(f, peak_rt_start_in_nearby_no, io_type);
	LoadSave::io(f, peak_rt_length, io_type);
	LoadSave::vector_io(f, nearby_intensity, io_type);
	LoadSave::vector_io(f, nearby_rt, io_type);
}

void PrecursorIon::set_peak_rt_in_nearby_no(int rt_start_no, int rt_length) {
	if (rt_start_no < 0)
		rt_start_no = 0;
	if (rt_start_no > nearby_intensity.size()) {
		rt_start_no = nearby_intensity.size() - 1;
		rt_length = 0;
	} else if (rt_start_no + rt_length > nearby_intensity.size()) {
		rt_length = nearby_intensity.size() - rt_start_no;
	}
	peak_rt_start_no = peak_rt_start_no + rt_start_no - peak_rt_start_in_nearby_no;
	peak_rt_start_in_nearby_no = rt_start_no;
	peak_rt_length = rt_length;
}

void PrecursorIon::set_peak_rt_no(int rt_start_no, int rt_length) {
	if (rt_start_no < 0)
		rt_start_no = 0;
	nearby_intensity = vector<double>(rt_length, 0);
	nearby_rt = vector<double>(rt_length, 0);

	peak_rt_start_in_nearby_no = 0;
	peak_rt_start_no = rt_start_no;
	peak_rt_length = rt_length;
}

int PrecursorIon::get_peak_apex_rt_in_nearby_no() const {
	if (nearby_intensity.empty())
		return 0;
	double max = nearby_intensity[peak_rt_start_in_nearby_no];
	int max_no = peak_rt_start_in_nearby_no;
	for (int i = 1; i < peak_rt_length; i++) {
		if (nearby_intensity[peak_rt_start_in_nearby_no + i] > max) {
			max_no = peak_rt_start_in_nearby_no + i;
			max = nearby_intensity[max_no];
		}
	}
	return max_no;

//	double sum = 0.0;
//	auto i = nearby_intensity.begin() + peak_rt_start_in_nearby_no;
//	int rt = 0;
//	for (; i != nearby_intensity.begin() + peak_rt_start_in_nearby_no + peak_rt_length; ++i, ++rt) {
//		sum += *i * rt;
//	}
//	double sum_intensity = accumulate(nearby_intensity.begin() + peak_rt_start_in_nearby_no,
//			nearby_intensity.begin() + peak_rt_start_in_nearby_no + peak_rt_length, 0.0);
//	if (sum_intensity > 0)
//		return peak_rt_start_in_nearby_no + sum / sum_intensity;
//	else
//		return peak_rt_start_in_nearby_no;
}

void PrecursorIon::get_nearby_intensity(vector<double>& intensity) const {
	intensity.clear();
	if (nearby_intensity.empty())
		intensity.push_back(0);
	else
		intensity = nearby_intensity;
}

void PrecursorIon::get_peak_rt_in_nearby_no(int& rt_start_in_nearby_no, int& length) const {
	rt_start_in_nearby_no = peak_rt_start_in_nearby_no;
	length = peak_rt_length;
}

void PrecursorIon::get_nearby_rt_no(int& rt_start_no, int& length) const {
	rt_start_no = peak_rt_start_no - peak_rt_start_in_nearby_no;
	length = nearby_intensity.size();
}

int PrecursorIon::get_peak_apex_rt_no() const {
	return peak_rt_start_no - peak_rt_start_in_nearby_no + get_peak_apex_rt_in_nearby_no();
}

void PrecursorIon::get_nearby_rt(vector<double>& rt) const {
	rt = nearby_rt;
}

void PrecursorIon::set_peak(double peak_start, double peak_end) {
	int peak_start_no = ToolsSpectrum::get_rt_no(nearby_rt, peak_start);
	int peak_end_no = ToolsSpectrum::get_rt_no(nearby_rt, peak_end) + 1;
	if (peak_start_no < 0)
		peak_start_no = 0;
	if (peak_end_no > nearby_rt.size())
		peak_end_no = nearby_rt.size();
	set_peak_rt_in_nearby_no(peak_start_no, peak_end_no - peak_start_no);
}

void PrecursorIon::remove_nearby_intensity() {
	vector<double>(nearby_intensity.begin() + peak_rt_start_in_nearby_no,
			nearby_intensity.begin() + peak_rt_start_in_nearby_no + peak_rt_length).swap(
			nearby_intensity);
	vector<double>(nearby_rt.begin() + peak_rt_start_in_nearby_no,
			nearby_rt.begin() + peak_rt_start_in_nearby_no + peak_rt_length).swap(nearby_rt);
	peak_rt_start_in_nearby_no = 0;
}
} /* namespace GroupDIA */
