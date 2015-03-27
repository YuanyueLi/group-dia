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

#include "AnalysisIdentification.h"

template<class T1, class T2> bool large_first(const pair<T1, T2>& i, const pair<T1, T2>& j) {
	return (i.first > j.first);
}

namespace GroupDIA {
AnalysisIdentification::AnalysisIdentification(const Param& _para, int windows_size,
		int cur_windows_num, int exp_size) {
	map_scan_pep.clear();
	para = _para;
	para_ide = para.copy("identification:", true);
	MIN_PRODUCT_IONS_NUM = para_ide.getValue("min_product_ions_num");
	this->windows_size = windows_size;
	this->cur_windows_num = cur_windows_num;
	this->exp_size = exp_size;
}

void AnalysisIdentification::reexact_intensity(string pepxml_filename, string input_chro_filename,
		string temp_chro_filename, vector<String> ms1_mzML_filename,
		vector<String> ms2_chro_filename, vector<String> ms2_mzML_filename,
		string output_chro_filename) {
////////////////////// Add precursor intensity //////////////////////////////////////////////
	MergeTransitions merge(para, para_ide, windows_size, cur_windows_num, exp_size);
	merge.add_identification(pepxml_filename);
	merge.choose_best_id();
	merge.load_transition(input_chro_filename, temp_chro_filename);

	cout << "Start extract intensity" << endl;
	Param ms1_para = para.copy("spectrum_ms1:", true);
	int exact_scan_num = para.getValue("exact_scan_num");
////////////////////// Add precursor intensity //////////////////////////////////////////////
	for (int cur_ms1_num = 0; cur_ms1_num < ms1_mzML_filename.size(); cur_ms1_num++) {
		string temp_chro_filename_out = temp_chro_filename + ".tmp";
		SwathExperiment exp(ms1_para);
		exp.set_exp(ms1_mzML_filename[cur_ms1_num]);
		fstream fi(temp_chro_filename.c_str(), fstream::in | fstream::binary);
		fstream fo(temp_chro_filename_out.c_str(), fstream::out | fstream::binary);
		omp_lock_t lock1, lock2;
		omp_init_lock(&lock1);
		omp_init_lock(&lock2);
		int batch_size = 1000;
		vector<TransitionIdentified> trans;
		while (!fi.eof()) {
			bool is_exist = true;
			//Load transition
			for (int i = 0; i < batch_size; i++) {
				trans.push_back(TransitionIdentified(para_ide));
				is_exist = trans[i].load_store(fi, LoadSave::IOType::Read,
						LoadSave::Type::Basic_Precursor_Iden);
				if (!is_exist) {
					trans.resize(i, TransitionIdentified(para_ide));
					break;
				}
			}

#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < trans.size(); i++) {
				trans[i].add_ms1_nearby(exp, cur_ms1_num, exact_scan_num);
			}

			//Store transition
			for (int i = 0; i < trans.size(); i++) {
				trans[i].load_store(fo, LoadSave::IOType::Write,
						LoadSave::Type::Basic_Precursor_Iden);
			}
			trans.clear();
		}
		fi.close(), fo.close();
		rename(temp_chro_filename_out.c_str(), temp_chro_filename.c_str());
	}

////////////////////// Add product intensity //////////////////////////////////////////////
	Param swath_para = para.copy("spectrum:", true);
	for (int cur_ms2_num = 0; cur_ms2_num < ms2_mzML_filename.size(); cur_ms2_num++) {
		SwathExperiment exp(swath_para);
		exp.set_exp(ms2_mzML_filename[cur_ms2_num]);
		fstream fi(temp_chro_filename.c_str(), fstream::in | fstream::binary);
		fstream fo(ms2_chro_filename[cur_ms2_num].c_str(), fstream::out | fstream::binary);
		omp_lock_t lock1, lock2;
		omp_init_lock(&lock1);
		omp_init_lock(&lock2);
		int batch_size = 1000;
		vector<TransitionIdentified> trans;
		while (!fi.eof()) {
			bool is_exist = true;
			//Load transition
			for (int i = 0; i < batch_size; i++) {
				trans.push_back(TransitionIdentified(para_ide));
				is_exist = trans[i].load_store(fi, LoadSave::IOType::Read,
						LoadSave::Type::Basic_Precursor_Iden);
				if (!is_exist) {
					trans.resize(i, TransitionIdentified(para_ide));
					break;
				}
			}

#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < trans.size(); i++) {
				trans[i].add_ms2_nearby(exp, cur_ms2_num);
			}

			//Store transition
			for (int i = 0; i < trans.size(); i++) {
				trans[i].load_store(fo, LoadSave::IOType::Write, LoadSave::Type::Product,
						cur_ms2_num);
			}
			trans.clear();
		}
		fi.close(), fo.close();
	}
}

void AnalysisIdentification::score_transition(string input_chro_filename, string temp_chro_filename,
		vector<String> ms2_chro_filename, vector<String> ms2_mzML_filename,
		string output_chro_filename, string output_score_filename) {
	ToolsTabFileWrite trans_score(output_score_filename);
////////////////////// Score transition //////////////////////////////////////////////
	cout << "Start score" << endl;
	int exp_size = ms2_mzML_filename.size();

	fstream fo(output_chro_filename.c_str(), fstream::out | fstream::binary);
	for (int cur_exp_num = 0; cur_exp_num < exp_size; cur_exp_num++) {
		//Score the result
		fstream ft;
		ft.open(temp_chro_filename.c_str(), fstream::in | fstream::binary);
		vector<fstream> fms2(exp_size);
		for (int i = 0; i < exp_size; i++) {
			fms2[i].open(ms2_chro_filename[i].c_str(), fstream::in | fstream::binary);
		}

		MzMLFile swath_file;
		MapType swath_map;
		if (true) {
			swath_file.load(ms2_mzML_filename[cur_exp_num], swath_map);
		}
		OpenSwath::SpectrumAccessPtr swath_ptr =
				SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);

		omp_lock_t lock1, lock2;
		omp_init_lock(&lock1);
		omp_init_lock(&lock2);
#pragma omp parallel
		{
			while (!ft.eof()) {
				TransitionIdentified t(para_ide);
				bool is_exist = true;

				//Load transition
				omp_set_lock(&lock1);
				is_exist = t.load_store(ft, LoadSave::IOType::Read,
						LoadSave::Type::Basic_Precursor_Iden);
				if (is_exist)
					for (int cur_exp_num = 0; cur_exp_num < exp_size; cur_exp_num++)
						if (!t.load_store(fms2[cur_exp_num], LoadSave::IOType::Read,
								LoadSave::Type::Product, cur_exp_num)) {
							throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
									"Error in load when select_ms2");
						}
				omp_unset_lock(&lock1);

				if (is_exist) {
					if (!t.is_exist())
						continue;
					t.reselect_product_ions(para);
					if (t.get_product_num() < MIN_PRODUCT_IONS_NUM or !t.is_exist())
						continue;

					TransitionScore openswath_score(t, para);
					openswath_score.score(cur_exp_num, swath_ptr, trans_score);

					if (cur_exp_num == 0) {
						//Store transition
						omp_set_lock(&lock2);
						t.load_store(fo, LoadSave::IOType::Write, LoadSave::Type::All);
						omp_unset_lock(&lock2);
					}
				}
			}
		}

		for (int i = 0; i < ms2_chro_filename.size(); i++) {
			fms2[i].close();
		}
		ft.close();
	}
	fo.close();
	for (int i = 0; i < ms2_chro_filename.size(); i++) {
		File::remove(ms2_chro_filename[i]);
	}
	File::remove(temp_chro_filename);
}
}/* namespace GroupDIA */
