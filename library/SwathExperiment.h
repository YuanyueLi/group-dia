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

#ifndef SWATHEXPERIMENT_H_
#define SWATHEXPERIMENT_H_
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilterAlgorithm.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <unordered_map>
#include <memory>
#include <vector>
#include <string>
using namespace std;
using namespace OpenMS;

namespace GroupDIA {

class SwathExperiment {

public:
	SwathExperiment(Param para);
	Param get_para() const;
	//set exp
	void set_exp(string mzml_filename);
	void set_exp_for_record();
	const MSSpectrum<Peak1D>& get_spec(int rt_no) const;
	int size() const;

	//Convert rt
	int get_rt_no(double rt) const;
	double get_rt(int rt_no) const;
	void get_rt(int rt_no_start, int rt_length, vector<double>& rt) const;

	//Record the intensity for extra
	void record_mz(int id, int rt_start_no, int rt_length, const vector<double>& mz);
	void record_mz(int id, int rt_no, const vector<double>& mz);

	void prepare_for_record();

	//get the intensity in specific time
	double get_intensity(double mz, int rt_no) const;
	void get_intensity(double mz, int rt_start_no, int length, vector<double>& intensity) const;
	//The mz should be sorted!
	void get_intensity(const vector<double>& mz, int rt_no, vector<double>& intensity) const;
	void get_record_intensity(int rt_no, vector<pair<double, double*> >& intensity) const;
	//The mz should be sorted!
	void get_record_intensity(int id, int rt_no, vector<double>& intensity) const;

	//Get the nearby intensity
	void get_nearby_intensity(double mz, int rt_start_no, int rt_length, int nearby_outside_length,
			vector<double>& nearby_intensity, int& rt_start_no_in_nearby) const;

	//Find peptide peak
	//If the intensity of peak in apex_rt_no is zero, return false, else return true.
	//If return false, peak_intensity is empty.
	bool find_peak(double mz, int apex_rt_no, int& peak_rt_no_start, int& peak_length,
			vector<double>& peak_intensity) const;

	//Get Ms2 intensity
	//The peak_intensity: mz-rt-intensity
	double get_mz_delta(double mz) const;

private:
	OpenMS::MSExperiment<Peak1D> exp;
	vector<double> retent_time;

	bool is_record_mode;
	vector<vector<vector<double> > > wanted_mz_intensity; //rt, id, mz
	map<int, int> map_id_num;

	bool use_ppm;
	double delta_ppm;
	double delta_mz;
	Param _para;

	bool find_apex(double mz, int apex_rt_no, int& peak_apex_rt_no, double& apex_intensity) const;
};

} /* namespace GroupDIA */
#endif /* SWATHEXPERIMENT_H_ */
