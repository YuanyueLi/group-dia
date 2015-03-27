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

#ifndef TOOLSFILENAME_H_
#define TOOLSFILENAME_H_
#include "ToolsSpectrum.h"
namespace GroupDIA {
class ToolsFilename {
public:
	ToolsFilename();
	int read_filelist(string list_filename);
	int get_swath_windows_size() const;

	void get_mz(int windows_num, double& min_mz, double& max_mz) const;
	string get_project_name() const;
	string get_result_path() const;

	//RT filename
	string get_rt_normalizer_filename(int ref_exp_num, int sam_exp_num) const;

	//Ms1 filename
	string get_ms1_mzML_filename(int exp_num) const;
	string get_ms1_feature_filename(int exp_num);

	//Ms2 filename
	string get_ms2_mzML_filename(int exp_num, int window_num) const;
	string get_ms2_splited_chro_filename(int ref_exp_num, int sam_exp_num, int window_num) const;
	string get_feature_chro_filename(int window_num) const;
	string get_ms2_chro_filename(int exp_num, int window_num) const;

	string get_mgf_filename(int exp_num, int window_num) const;
	string get_mzML_filename(int exp_num, int window_num) const;

	//Identification
	string get_identification_result_chro_filename(int window_num);
	string get_identification_chro_filename(int window_num);
	string get_identification_splited_chro_filename(int exp_num, int window_num) const;
	string get_identification_result(int window_num) const;

	string get_openswath_result(int exp_num, int window_num) const;

	string get_mprohet_input_filename() const;
	string get_mprohet_output_filename() const;
	string get_mprohet_input_filename(int exp_num) const;
	string get_mprohet_output_filename(int exp_num) const;
	string get_mprohet_2nd_input_filename() const;
	string get_mprohet_2nd_output_filename() const;

	int read_mzXML_file(string list_filename);
	string get_mzXML_filename(int num) const;
private:
	vector<string> file_pre;
	vector<pair<double, double> > windows;
	String extension;
	String result_file_path;
	String temp_file_path;
	String project_name;

	int load_windows(string win_filename);
};
} /* namespace GroupDIA */
#endif /* TOOLSFILENAME_H_ */
