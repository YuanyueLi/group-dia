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

#ifndef PRECURSOR_H_
#define PRECURSOR_H_

#include "RTNormalizer.h"
#include "ToolsLoadSave.h"
#include "SwathExperiment.h"
namespace GroupDIA {

class PrecursorIon {
public:
	PrecursorIon();

	bool set_ms1(const SwathExperiment& exp, double mz, double apex_rt, int exact_scan_num);
	void add_ms1(const SwathExperiment& exp, double mz, const RTNormalizer& rt, int exact_scan_num,
			PrecursorIon& result) const;
	void add_ms1(const SwathExperiment& exp, double mz, int real_rt_no, int rt_length,
			int exact_scan_num, PrecursorIon& result) const;

	void set_peak(const SwathExperiment& exp, double mz, int peak_rt_no_start, int peak_rt_length,
			int exact_scan_num);
	void set_peak(double peak_start, double peak_end);
	void set_peak_rt_in_nearby_no(int rt_start_no, int rt_length);
	void set_peak_rt_no(int rt_start_no, int rt_length);
	void remove_nearby_intensity();

	//Information
	bool is_exist() const;
	void not_exist();
	void get_peak_intensity(vector<double>&intensity) const;
	void get_nearby_intensity(vector<double>&intensity) const;
	double get_peak_apex_rt() const;
	int get_peak_apex_rt_in_nearby_no() const;
	int get_peak_apex_rt_no() const;
	void get_peak_rt(double& rt_start, double& rt_apex, double& rt_end) const;

	void get_peak_rt_no(int& rt_start_no, int&length) const;
	void get_nearby_rt_no(int& rt_start_no, int&length) const;
	void get_peak_rt_in_nearby_no(int& rt_start_in_nearby_no, int&length) const;
	void get_nearby_rt(vector<double>& rt) const;

	//load and save
	void load_store(fstream& fo, int type);

private:
	//Peak Info
	int peak_rt_start_no;
	int peak_rt_length;

	//Nearby Info
	int peak_rt_start_in_nearby_no;
	vector<double> nearby_intensity;
	vector<double> nearby_rt;

	int get_nearby_rt_start_no() const {
		return peak_rt_start_no - peak_rt_start_in_nearby_no;
	}

};

} /* namespace GroupDIA */
#endif /* PRECURSOR_H_ */
