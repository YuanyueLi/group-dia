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

#include "ToolsFilename.h"

namespace GroupDIA {

int ToolsFilename::read_filelist(string list_filename) {
	fstream fin(list_filename.c_str(), fstream::in);
	getline(fin, project_name);
	project_name.trim();

	//Deal windows
	String line;
	getline(fin, line);
	load_windows(line);

	//Deal result file path
	getline(fin, result_file_path);
	result_file_path.trim();
	result_file_path = result_file_path.ensureLastChar('/');
	if (!File::isDirectory(result_file_path)) {
		throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, result_file_path);
	}

	//Deal temp file path
	getline(fin, temp_file_path);
	temp_file_path.trim();
	temp_file_path = temp_file_path.ensureLastChar('/');
	if (!File::isDirectory(temp_file_path)) {
		throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, temp_file_path);
	}

	while (getline(fin, line)) {
		line.trim();
		//Determine type: mzML or mzML.bz2
		auto cur_name = line;
		extension = "";
		int dot_place = cur_name.rfind('.');
		if (dot_place != string::npos) {
			while (cur_name.substr(dot_place) != ".mzML") {
				cur_name = cur_name.substr(0, dot_place);
				dot_place = cur_name.rfind('.');
			}
			extension = line.substr(dot_place);
			String cur_filename = line.substr(0, dot_place);
			cur_filename = cur_filename.prefix(cur_filename.size() - 8);
			file_pre.push_back(cur_filename);
		} else {
			throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
					"No dot in filename: " + line);
		}
	}
	return file_pre.size();
}

int ToolsFilename::read_mzXML_file(string list_filename) {
	fstream fin(list_filename.c_str(), fstream::in);
	string line;
	while (fin >> line) {
		String cur_filename = File::removeExtension(line);
		file_pre.push_back(cur_filename);
	}
	return file_pre.size();
}

string ToolsFilename::get_mgf_filename(int exp_num, int window_num) const {
	if (exp_num == -1)
		return result_file_path + project_name + "." + String(window_num) + ".mgf";
	else
		return result_file_path + project_name + "." + String(exp_num) + "_" + String(window_num)
				+ ".mgf";
}

string ToolsFilename::get_mzML_filename(int exp_num, int window_num) const {
	if (exp_num == -1)
		return result_file_path + project_name + "." + String(window_num) + ".mzML";
	else
		return result_file_path + project_name + "." + String(exp_num) + "_" + String(window_num)
				+ ".mzML";
}

int ToolsFilename::load_windows(string win_filename) {
	windows.clear();
	fstream fin(win_filename.c_str(), fstream::in);
	double a, b, c;
	while (fin >> a >> b >> c) {
		windows.push_back(pair<double, double>(a, b));
	}
	return windows.size();
}

string ToolsFilename::get_ms1_mzML_filename(int exp_num) const {
	if (exp_num > file_pre.size()) {
		throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Experiment number Error", String(exp_num));
	}

	return file_pre[exp_num] + "_ms1scan" + extension;
}

string ToolsFilename::get_ms1_feature_filename(int exp_num) {
	if (exp_num > file_pre.size()) {
		throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Experiment number Error", String(exp_num));
	}

	return file_pre[exp_num] + "_ms1scan.featureXML";
}

string ToolsFilename::get_ms2_mzML_filename(int exp_num, int window_num) const {
	if (exp_num > file_pre.size()) {
		throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Experiment number Error", String(exp_num));
	}
	string w = "";
	if (window_num < 10) {
		w = "0" + String(window_num);
	} else {
		w = String(window_num);
	}

	return file_pre[exp_num] + "_" + w + extension;
}

//string ToolsFilename::get_ms2_splited_chro_filename(int exp_num, int window_num) const {
//	if (exp_num > file_pre.size()) {
//		throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
//				"Experiment number Error", String(exp_num));
//	}
//
//	string result_file_name = project_name + "." + String(exp_num) + "_" + String(window_num)
//			+ ".part.chro";
//	return temp_file_path + File::basename(result_file_name);
//}

string ToolsFilename::get_ms2_splited_chro_filename(int ref_exp_num, int sam_exp_num,
		int window_num) const {
	string result_file_name;
	if (ref_exp_num == -1)
		result_file_name = project_name + "." + String(sam_exp_num) + "_" + String(window_num)
				+ ".part.chro";
	else
		result_file_name = project_name + "." + String(ref_exp_num) + "_" + String(sam_exp_num)
				+ "_" + String(window_num) + ".part.chro";

	return temp_file_path + result_file_name;
}

string ToolsFilename::get_identification_splited_chro_filename(int exp_num, int window_num) const {
	if (exp_num > file_pre.size()) {
		throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Experiment number Error", String(exp_num));
	}

	string result_file_name = project_name + "." + String(window_num) + "_" + String(exp_num)
			+ ".part.id.chro";
	return temp_file_path + File::basename(result_file_name);
}

string ToolsFilename::get_rt_normalizer_filename(int ref_exp_num, int sam_exp_num) const {
	if (ref_exp_num > file_pre.size() or sam_exp_num > file_pre.size()) {
		throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Experiment number Error", String(ref_exp_num) + "_" + String(sam_exp_num));
	}

	string result_file_name = file_pre[ref_exp_num] + "_" + String(sam_exp_num) + ".rt";
	return result_file_path + File::basename(result_file_name);

}

void ToolsFilename::get_mz(int windows_num, double& min_mz, double& max_mz) const {
	min_mz = windows[windows_num].first;
	max_mz = windows[windows_num].second;
}

int ToolsFilename::get_swath_windows_size() const {
	return windows.size();
}

string ToolsFilename::get_ms2_chro_filename(int exp_num, int window_num) const {
	string result_file_name;
	if (exp_num != -1)
		result_file_name = project_name + "_" + String(exp_num) + "_" + String(window_num)
				+ ".ms2.chro";
	else
		result_file_name = project_name + "." + String(window_num) + ".ms2.chro";

	return result_file_path + result_file_name;
}

string ToolsFilename::get_feature_chro_filename(int window_num) const {
	string result_file_name;
	result_file_name = project_name + "." + String(window_num) + ".fea.chro";
	return temp_file_path + result_file_name;
}

string ToolsFilename::get_identification_result_chro_filename(int window_num) {
	string result_file_name = project_name + "." + String(window_num) + ".id.chro";
	return result_file_path + result_file_name;
}

string ToolsFilename::get_identification_chro_filename(int window_num) {
	string result_file_name = project_name + "." + String(window_num) + ".part.id.chro";
	return temp_file_path + File::basename(result_file_name);
}

string ToolsFilename::get_project_name() const {
	return project_name;
}

string ToolsFilename::get_result_path() const {
	return result_file_path;
}

ToolsFilename::ToolsFilename() {
}

string ToolsFilename::get_identification_result(int window_num) const {
	string result_file_name = project_name + "." + String(window_num) + ".txt";
	return result_file_path + result_file_name;
}

string ToolsFilename::get_mzXML_filename(int num) const {
	return file_pre[num];
}

string ToolsFilename::get_openswath_result(int exp_num, int window_num) const {
	string result_file_name = project_name + "_openswath" + String(exp_num) + "_"
			+ String(window_num) + ".csv";
	return result_file_path + result_file_name;
}

string ToolsFilename::get_mprohet_input_filename() const {
	string result_file_name = project_name + "_mprohet" + ".txt";
	return result_file_path + result_file_name;
}

string ToolsFilename::get_mprohet_output_filename() const {
	string result_file_name = project_name + "_mprohet" + "_with_dscore.csv";
	return result_file_path + result_file_name;
}

string ToolsFilename::get_mprohet_input_filename(int exp_num) const {
	string result_file_name = project_name + "_mprohet_" + String(exp_num) + ".txt";
	return result_file_path + result_file_name;
}

string ToolsFilename::get_mprohet_output_filename(int exp_num) const {
//	string result_file_name = project_name + "_mprohet_" + String(exp_num) + "_with_dscore.csv";
	string result_file_name = project_name + "_mprohet_out_" + String(exp_num) + ".txt";
	return result_file_path + result_file_name;
}

string ToolsFilename::get_mprohet_2nd_input_filename() const {
	string result_file_name = project_name + "_2_mprohet" + ".txt";
	return result_file_path + result_file_name;
}

string ToolsFilename::get_mprohet_2nd_output_filename() const {
	string result_file_name = project_name + "_2_mprohet" + "_with_dscore.csv";
	return result_file_path + result_file_name;
}
}/* namespace GroupDIA */
