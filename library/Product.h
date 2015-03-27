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

#ifndef PRODUCT_H_
#define PRODUCT_H_
#include "RTNormalizer.h"
#include "ToolsLoadSave.h"
#include "SwathExperiment.h"

using namespace GroupDIA;
namespace GroupDIA {
struct ProductIonIntensity {
	vector<double> intensity;
	vector<double> get_intensity() const {
		return intensity;
	}
	ProductIonIntensity() {
		intensity = vector<double>(1, 0);
	}
};

//该类实现以下功能：
//根据Ms2类，提取对应的子离子的信息
class ProductIons: public vector<ProductIonIntensity> {
public:
	ProductIons();

	void add_ms2(const SwathExperiment& exp, int scan_num, const vector<double>& mz,
			int rt_start_no, int length);
	void add_ms2(const SwathExperiment& exp, const vector<double>& mz, int rt_start_no, int length);
	void remove_nearby_intensity(int peak_start_no, int peak_length);

	//Get the intensity of specific product ion
	void get_intensity(int product_ion_no, int rt_no_from, int length,
			vector<double>& intensity) const;
	void get_intensity(int product_ion_no, vector<double>& intensity) const;
	void get_intensity(vector<vector<double> >& intensity) const;
	void get_intensity(int rt_no_from, int length, vector<vector<double> >& intensity) const;

	void refine_product_ions(const vector<int>& remained_place);
	void get_rt(vector<double>& rt) const;
	int get_rt_no(double rt) const;
	int get_rt_size(double rt_start, double rt_end) const;

	void load_store(fstream& fo, int type);
protected:
	vector<double> rt;
};

} /* namespace GroupDIA */
#endif /* PRODUCT_H_ */
