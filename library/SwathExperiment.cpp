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

#include "SwathExperiment.h"
#include "ToolsSpectrum.h"

template<class T1, class T2> bool small_first(const pair<T1, T2>& i, const pair<T1, T2>& j) {
	return (i.first < j.first);
}

namespace GroupDIA {

int SwathExperiment::get_rt_no(double rt) const {
	return ToolsSpectrum::get_rt_no(retent_time, rt);
}

double SwathExperiment::get_rt(int rt_no) const {
	return ToolsSpectrum::get_rt(retent_time, rt_no);
}

void SwathExperiment::set_exp(string mzml_filename) {
	if (!File::exists(mzml_filename))
		throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, mzml_filename);

	MzMLFile mzml;
	mzml.load(mzml_filename, exp);

	retent_time.resize(exp.size());
	auto iadd = retent_time.begin();
	for (auto it = exp.begin(); it != exp.end(); ++it, ++iadd) {
		*iadd = it->getRT();
	}
}

int SwathExperiment::size() const {
	return exp.size();
}

double SwathExperiment::get_intensity(double mz, int rt_no) const {
	if (rt_no < 0 or rt_no >= exp.size())
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Error in rt_no: " + String(__FILE__) + String(__LINE__)
						+ String(__PRETTY_FUNCTION__));

	//Get the max intensity between [ mz-delta,mz+delta ]
	double delta = get_mz_delta(mz);
	double mz_l = mz - delta;
	double mz_r = mz + delta;

	const SpectrumT& spec = exp[rt_no];
//	vector<double> candidate;
	double intensity = 0;
	for (auto it = spec.begin(); it != spec.end(); ++it) {
		double cmz = it->getMZ();
		if (cmz >= mz_l) {
			if (cmz > mz_r)
				break;
//			candidate.push_back(it->getIntensity());
			intensity += it->getIntensity();
		}
	}
//	if (candidate.empty())
//		return 0;
//	else
	//change here!
//		return *max_element(candidate.begin(), candidate.end());
	return intensity;
}

void SwathExperiment::get_intensity(const vector<double>& mz, int rt_no,
		vector<double>& intensity) const {
//	sort(mz.begin(), mz.end());
	if (is_record_mode) {
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Not in record_mode");
	}
	if (rt_no < 0 or rt_no >= exp.size()) {
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Error in rt_no:\t" + String(rt_no) + "\t" + String(exp.size()) + "\t"
						+ String(__FILE__) + String(__LINE__) + String(__PRETTY_FUNCTION__));
	}
//	if (mz.empty())
//		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
//				"Error in mz: " + String(__FILE__) + String(__LINE__)
//						+ String(__PRETTY_FUNCTION__));

	//Get the all intensity between [ mz-delta,mz+delta ]
	const SpectrumT& spec = exp[rt_no];
	intensity.resize(mz.size(), 0);
	auto imz = mz.begin();
	auto iint = intensity.begin();
	auto ispec = spec.begin();

	while (imz != mz.end() && ispec != spec.end()) {
		double delta = get_mz_delta(*imz);
		double mz_l = *imz - delta;
		double mz_r = *imz + delta;
		if (ispec->getMZ() < mz_l) {
			++ispec;
			continue;
		}
		if (ispec->getMZ() > mz_r) {
			++imz, ++iint;
			continue;
		}
		auto it = ispec;
		while (it != spec.end()) {
			if (it->getMZ() > mz_r)
				break;
//			if (*iint < it->getIntensity())
			*iint += it->getIntensity();
			++it;
		}
		++imz, ++iint;
	}
}

void SwathExperiment::get_record_intensity(int rt_no,
		vector<pair<double, double*> >& mz_intensity) const {
	if (rt_no < 0 or rt_no >= exp.size()) {
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Error in rt_no:\t" + String(rt_no) + "\t" + String(exp.size()) + "\t"
						+ String(__FILE__) + String(__LINE__) + String(__PRETTY_FUNCTION__));
	}
	//Get the all intensity between [ mz-delta,mz+delta ]
	const SpectrumT& spec = exp[rt_no];

	auto imz = mz_intensity.begin();
	auto ispec = spec.begin();

	while (imz != mz_intensity.end() && ispec != spec.end()) {
		double delta = get_mz_delta(imz->first);
		double mz_l = imz->first - delta;
		double mz_r = imz->first + delta;
		if (ispec->getMZ() < mz_l) {
			++ispec;
			continue;
		}
		if (ispec->getMZ() > mz_r) {
			*(imz->second) = 0;
			++imz;
			continue;
		}
		auto it = ispec;
		*(imz->second) = 0;
		while (it != spec.end()) {
			if (it->getMZ() > mz_r)
				break;
			*(imz->second) += it->getIntensity();
			++it;
		}
		++imz;
	}
}

bool SwathExperiment::find_peak(double mz, int apex_rt_no, int& peak_rt_no_from, int& peak_length,
		vector<double>& peak_intensity) const {
	//Get the intensity in apex
	if (apex_rt_no < 0 || apex_rt_no >= exp.size()) {
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Apex rt no error");
	}

	peak_intensity.clear();
	double int_peak;

	bool exist = find_apex(mz, apex_rt_no, apex_rt_no, int_peak);
	if (!exist) {
		peak_rt_no_from = apex_rt_no;
		peak_length = 0;
		return false;
	}

	int rt_c = apex_rt_no - 1;
	double int_c = int_peak;
	double int_old = int_peak;
	list<double> peak_int;
	peak_int.push_back(int_peak);
//From peak to left
	while (rt_c >= 0) {
		int_c = get_intensity(mz, rt_c);
		if (int_c < 0) {
			break;
		} else {
			if (int_c >= int_old) {
				int_c++;
				break;
			}
		}
		peak_int.push_front(int_c);
		int_old = int_c;
		rt_c--;
	}

	peak_rt_no_from = rt_c + 1;
//From peak to right
	int_old = int_peak;
	rt_c = apex_rt_no + 1;
	while (rt_c < exp.size()) {
		int_c = get_intensity(mz, rt_c);
		if (int_c < 0) {
			break;
		} else {
			if (int_c >= int_old) {
				int_c--;
				break;
			}
		}
		peak_int.push_back(int_c);
		int_old = int_c;
		rt_c++;
	}

	peak_length = peak_int.size();
	if (peak_length <= 2) {
		peak_rt_no_from = apex_rt_no;
		peak_length = 0;
		return false;
	}

	peak_intensity.insert(peak_intensity.end(), peak_int.begin(), peak_int.end());
	return true;
}

void SwathExperiment::get_rt(int rt_no_start, int rt_length, vector<double>& rt) const {
	rt.clear();
	rt.reserve(rt_length);
	for (int i = rt_no_start; i < rt_no_start + rt_length; i++) {
		rt.push_back(exp[i].getRT());
	}
}

bool SwathExperiment::find_apex(double mz, int apex_rt_no, int& peak_apex_rt_no,
		double& int_peak) const {
	peak_apex_rt_no = apex_rt_no;

	//determine if apex is exist
	if (apex_rt_no >= this->exp.size() - 1) {
		peak_apex_rt_no = this->exp.size() - 1;
		return false;
	}
	if (apex_rt_no == 0) {
		peak_apex_rt_no = 0;
		return false;
	}

	int_peak = get_intensity(mz, apex_rt_no);
	if (int_peak == 0) {
		return false;
	}
	int left = apex_rt_no - 1;
	int right = apex_rt_no + 1;
	double int_left = get_intensity(mz, left);
	double int_right = get_intensity(mz, right);
	if (int_peak >= int_left and int_peak > int_right) {
		peak_apex_rt_no = apex_rt_no;
		return true;
	} else if (int_left >= int_right) {
		while (left >= 0) {
			int_left = get_intensity(mz, left);
			if (int_left < int_peak) {
				peak_apex_rt_no = apex_rt_no;
				return true;
			} else {
				int_peak = int_left;
				--apex_rt_no, --left;
			}
		}
	} else {
		while (right < exp.size()) {
			int_right = get_intensity(mz, right);
			if (int_right < int_peak) {
				peak_apex_rt_no = apex_rt_no;
				return true;
			} else {
				int_peak = int_right;
				++apex_rt_no, ++right;
			}
		}
	}
	return false;
}

void SwathExperiment::get_nearby_intensity(double mz, int rt_start_no, int rt_length,
		int nearby_outside_length, vector<double>& nearby_intensity,
		int& rt_start_no_in_nearby) const {
	if (rt_start_no < 0 or rt_length < 0) {
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Error in rt_no: " + String(__FILE__) + String(__LINE__)
						+ String(__PRETTY_FUNCTION__));
	}
	if (rt_start_no + rt_length > exp.size()) {
		rt_length = exp.size() - rt_start_no;
	}

	int l = rt_start_no - nearby_outside_length;
	int r = rt_start_no + rt_length + nearby_outside_length;
	if (l < 0) {
		rt_start_no_in_nearby = rt_start_no;
		l = 0;
	} else {
		rt_start_no_in_nearby = rt_start_no - l;
	}

	if (r > exp.size())
		r = exp.size();
	get_intensity(mz, l, r - l, nearby_intensity);
//	smooth_intensity(nearby_intensity);
}

void SwathExperiment::get_intensity(double mz, int rt_start_no, int length,
		vector<double>& intensity) const {
	if (rt_start_no < 0 or rt_start_no + length > exp.size() or length < 0)
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Error in rt_no: " + String(__FILE__) + String(__LINE__)
						+ String(__PRETTY_FUNCTION__));
	intensity.clear();
	intensity.reserve(length);
	for (int rt = rt_start_no; rt < rt_start_no + length; ++rt)
		intensity.push_back(get_intensity(mz, rt));

}

SwathExperiment::SwathExperiment(Param para) {
	_para = para;
	use_ppm = para.getValue("use_ppm").toBool();
	if (use_ppm)
		delta_ppm = para.getValue("delta_ppm");
	else
		delta_mz = para.getValue("delta_mz");
	is_record_mode = false;
}

Param SwathExperiment::get_para() const {
	return _para;
}
double SwathExperiment::get_mz_delta(double mz) const {
	if (use_ppm)
		return mz * delta_ppm / 1000000;
	else
		return delta_mz;
}

void SwathExperiment::prepare_for_record() {
	this->wanted_mz_intensity.clear();
	this->wanted_mz_intensity.resize(exp.size());
	this->map_id_num.clear();
	this->is_record_mode = false;
}

void SwathExperiment::record_mz(int id, int rt_start_no, int rt_length, const vector<double>& mz) {
	for (int i = 0; i < rt_length; i++)
		record_mz(id, rt_start_no + i, mz);
}

void SwathExperiment::record_mz(int id, int rt_no, const vector<double>& mz) {
	if (rt_no >= wanted_mz_intensity.size()) {
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Error in record_mz");
	}

#pragma omp critical
	{
		auto it = map_id_num.find(id);
		int insert_id;
		if (it == map_id_num.end()) {
			insert_id = map_id_num.size();
			map_id_num.insert(pair<int, int>(id, insert_id));
		} else
			insert_id = it->second;

		if (wanted_mz_intensity[rt_no].size() <= insert_id)
			wanted_mz_intensity[rt_no].resize(insert_id + 1);
		wanted_mz_intensity[rt_no][insert_id] = mz;
	}
}

void SwathExperiment::set_exp_for_record() {
	is_record_mode = false;

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < wanted_mz_intensity.size(); i++) {
		vector<pair<double, double*> > wanted_item;
		auto& cur_mz_intensity = wanted_mz_intensity[i];
		for (auto it = cur_mz_intensity.begin(); it != cur_mz_intensity.end(); ++it) {
			for (vector<double>::iterator jt = it->begin(); jt != it->end(); ++jt) {
				wanted_item.push_back(pair<double, double*>(*jt, &(*jt)));
			}
		}

		sort(wanted_item.begin(), wanted_item.end(), small_first<double, double*>);
		get_record_intensity(i, wanted_item);
	}

	is_record_mode = true;
}

void SwathExperiment::get_record_intensity(int id, int rt_no, vector<double>& intensity) const {
	auto it = map_id_num.find(id);
	if (rt_no >= wanted_mz_intensity.size() or rt_no < 0 or it == map_id_num.end()
			or it->second >= wanted_mz_intensity[rt_no].size()) {
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Error in get record intensity");
	}
	intensity = wanted_mz_intensity[rt_no][it->second];
}

const MSSpectrum<Peak1D>& SwathExperiment::get_spec(int rt_no) const {
	return exp[rt_no];
}

} /* namespace GroupDIA */
