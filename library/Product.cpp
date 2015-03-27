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

#include "ToolsLoadSave.h"
#include "Product.h"

namespace GroupDIA {

ProductIons::ProductIons() {
	rt.clear();
	this->clear();
}

int ProductIons::get_rt_size(double rt_start, double rt_end) const {
	int rt_start_no = get_rt_no(rt_start);
	int rt_end_no = get_rt_no(rt_end);
	return rt_end_no - rt_start_no;
}

void ProductIons::load_store(fstream& f, int type) {
	LoadSave::vector_io(f, this->rt, type);
	int size = this->size();
	LoadSave::io(f, size, type);
	if (size > this->size()) {
		this->resize(size);
	}
	for (int i = 0; i < size; i++) {
		ProductIonIntensity& t = this->at(i);
		LoadSave::vector_io(f, t.intensity, type);
	}
}

void ProductIons::add_ms2(const SwathExperiment& exp, const vector<double>& mz, int rt_start_no,
		int length) {
	this->resize(mz.size());
	if (rt_start_no < 0 or length < 0) {
		cout << "Error in add_ms2\n";
		return;
	}
	//add mz
	for (int i = 0; i < mz.size(); i++) {
		this->at(i).intensity.resize(length);
	}

	//add rt
	this->rt.resize(length);
#pragma omp parallel for schedule(dynamic)
	for (int i = rt_start_no; i < rt_start_no + length; i++) {
		vector<double> intensity;
		this->rt[i - rt_start_no] = exp.get_rt(i);
		exp.get_intensity(mz, i, intensity);
		for (int j = 0; j < mz.size(); j++) {
			this->at(j).intensity[i - rt_start_no] = intensity[j];
		}
	}
}

void ProductIons::add_ms2(const SwathExperiment& exp, int scan_num, const vector<double>& mz,
		int rt_start_no, int length) {
	this->resize(mz.size());
	if (rt_start_no < 0 or length < 0) {
		cout << "Error in add_ms2\n";
		return;
	}
	//add mz
	for (int i = 0; i < mz.size(); i++) {
		this->at(i).intensity.resize(length);
	}

	//add rt
	this->rt.resize(length);
#pragma omp parallel for schedule(dynamic)
	for (int i = rt_start_no; i < rt_start_no + length; i++) {
		vector<double> intensity;
		this->rt[i - rt_start_no] = exp.get_rt(i);
//		exp.get_intensity(mz, i, intensity);
		exp.get_record_intensity(scan_num, i, intensity);
		for (int j = 0; j < mz.size(); j++) {
			this->at(j).intensity[i - rt_start_no] = intensity[j];
		}
	}
}

void ProductIons::get_rt(vector<double>& _rt) const {
	_rt = this->rt;
}

void ProductIons::refine_product_ions(const vector<int>& remained_place) {
	ProductIons result;
	result.rt = this->rt;
	result.resize(remained_place.size());
	for (int i = 0; i < result.size(); i++) {
		result[i] = this->at(remained_place[i]);
	}
	result.shrink_to_fit();
	*this = result;
}

void ProductIons::get_intensity(int product_ion_no, int rt_no_from, int length,
		vector<double>& intensity) const {
	intensity.clear();
	if (length == 0)
		intensity.push_back(0);
	else
		intensity.insert(intensity.end(), this->at(product_ion_no).intensity.begin() + rt_no_from,
				this->at(product_ion_no).intensity.begin() + rt_no_from + length);
}

void ProductIons::get_intensity(int product_ion_no, vector<double>& intensity) const {
	OPENMS_PRECONDITION(product_ion_no<this->size(),"product_ion_no too large");
	intensity = this->at(product_ion_no).intensity;
}

void ProductIons::get_intensity(vector<vector<double> >& intensity) const {
	intensity.clear();
	intensity.reserve(this->size());
	for (auto it = this->begin(); it != this->end(); ++it)
		intensity.push_back(it->intensity);
}

void ProductIons::get_intensity(int rt_no_from, int length,
		vector<vector<double> >& intensity) const {
	intensity.clear();
	intensity.reserve(this->size());
	for (auto it = this->begin(); it != this->end(); ++it)
		intensity.push_back(
				vector<double>(it->intensity.begin() + rt_no_from,
						it->intensity.begin() + rt_no_from + length));
}

int ProductIons::get_rt_no(double _rt) const {
	return ToolsSpectrum::get_rt_no(rt, _rt);
}

void ProductIons::remove_nearby_intensity(int peak_start_no, int peak_length) {
	vector<double>(rt.begin() + peak_start_no, rt.begin() + peak_start_no + peak_length).swap(rt);
	for (auto it = this->begin(); it != this->end(); ++it) {
		vector<double>(it->intensity.begin() + peak_start_no,
				it->intensity.begin() + peak_start_no + peak_length).swap(it->intensity);
	}
}
}/* namespace GroupDIA */
