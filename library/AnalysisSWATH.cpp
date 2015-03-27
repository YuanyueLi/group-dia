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

#include "AnalysisSWATH.h"

namespace GroupDIA {

} /* namespace GroupDIA */

const int BATCH = 100;

GroupDIA::AnalysisSWATH::AnalysisSWATH(const Param& _para) {
	MIN_ALLOWED_PRODUCT_NUM = _para.getValue("similar:min_allowed_product_ion_numbers");
	RT_WINDOWS = _para.getValue("similar:retention_window");
	MIN_ALLOWED_CORRELATION = _para.getValue("similar:min_allowed_correlation_in_crude_refine");

	MAX_ALLOWED_DELAY = _para.getValue("delay:max_allowd_scan_num_delay_in_normalized");
	FINE_NEARBY = _para.getValue("delay:max_allowd_scan_num_delay_in_second_normalized");
	MIN_INTENSITY = _para.getValue("delay:min_intensity_in_calc_delay");

	para = _para.copy("similar:", true);
}

bool GroupDIA::AnalysisSWATH::set_ms1(const Feature& feature, int exp_num, int scan_num,
		const SwathExperiment& exp) {

	TransitionGroup* t = new TransitionGroup();
	bool right = t->set_transition(exp, feature, exp_num, scan_num, RT_WINDOWS);
	if (right) {
#pragma omp critical
		trans.push_back(t);
	} else {
		delete t;
	}

	return right;
}

void GroupDIA::AnalysisSWATH::add_ms1(int ref_exp_num, int add_exp_num, const SwathExperiment& exp,
		string rt_filename) {
	bool not_use = true;
	for (int i = 0; i < size(); i++) {
		if (trans[i]->get_feature_exp_num() == ref_exp_num) {
			not_use = false;
			break;
		}
	}
	if (not_use)
		return;

	RTNormalizer rt;
	rt.load(rt_filename);

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < size(); i++) {
		if (trans[i]->get_feature_exp_num() == ref_exp_num)
			trans[i]->add_ms1(exp, rt, add_exp_num, RT_WINDOWS);
	}
}

void GroupDIA::AnalysisSWATH::set_ms2(int exp_num, SwathExperiment& ms2_exp, double ms1_win_start,
		double ms1_win_end, string output_chro_filename) {
	bool not_use = true;
	for (int i = 0; i < size(); i++) {
		if (trans[i]->get_feature_exp_num() == i) {
			not_use = false;
			break;
		}
	}
	if (not_use)
		return;

	fstream fo(output_chro_filename.c_str(), fstream::out | fstream::binary | fstream::app);

	int istart = 0, iend = 0;
	while (true) {
		istart = iend;
		iend += BATCH;
		if (iend > size())
			iend = size();
		if (istart >= size())
			break;

#pragma omp parallel for schedule(dynamic)
		for (int i = istart; i < iend; i++) {
			if (trans[i]->get_feature_exp_num() != exp_num)
				continue;
			trans[i]->set_product_ions_mz(ms2_exp, ms1_win_start, ms1_win_end);
			trans[i]->add_feature_product(ms2_exp);
			trans[i]->crude_refine(MIN_ALLOWED_CORRELATION);
		}

		for (int i = istart; i < iend; i++) {
			if (trans[i]->get_product_num() >= MIN_ALLOWED_PRODUCT_NUM)
				trans[i]->load_store(fo, LoadSave::IOType::Write, LoadSave::Type::Feature_Product);
			trans[i]->remove_feature_product();
		}
	}
	fo.close();
}

void GroupDIA::AnalysisSWATH::refine() {
	vector<TransitionGroup*> new_trans;
	for (int i = 0; i < trans.size(); i++) {
		if (trans[i]->get_product_num() >= MIN_ALLOWED_PRODUCT_NUM)
			new_trans.push_back(trans[i]);
	}
	new_trans.swap(trans);
}

void GroupDIA::AnalysisSWATH::add_ms2(int exp_num, SwathExperiment& ms2_exp,
		string input_chro_filename, string output_chro_filename) {
	fstream fi(input_chro_filename.c_str(), fstream::in | fstream::binary);
	fstream fo(output_chro_filename.c_str(), fstream::out | fstream::binary);

	int istart = 0, iend = 0;
	while (true) {
		istart = iend;
		iend += BATCH;
		if (iend > size())
			iend = size();
		if (istart >= size())
			break;

		for (int i = istart; i < iend; i++) {
			while (true) {
				bool is_load = trans[i]->load_store(fi, LoadSave::IOType::Read,
						LoadSave::Type::Feature_Product);
				if (is_load)
					break;
				if (fi.eof()) {
					throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
							"Error in load feature product ions when add_ms2, trans_num: "
									+ String(i) + ", size: " + String(size()));
				}
			}
		}

#pragma omp parallel for schedule(dynamic)
		for (int i = istart; i < iend; i++) {
			trans[i]->add_ms2_nearby(ms2_exp, exp_num);
			trans[i]->refined_by_crosscorrelation(exp_num, MAX_ALLOWED_DELAY, FINE_NEARBY,
					MIN_INTENSITY);
			trans[i]->remove_nearby_intensity(exp_num);
		}

		for (int i = istart; i < iend; i++) {
			trans[i]->load_store(fo, LoadSave::IOType::Write, LoadSave::Type::Precursor, exp_num);
			trans[i]->load_store(fo, LoadSave::IOType::Write, LoadSave::Type::Product, exp_num);
			trans[i]->remove_ions();
			trans[i]->remove_feature_product();
		}
	}
	fi.close();
	fo.close();
}

void GroupDIA::AnalysisSWATH::select_ms2(vector<string> input_chro_filename,
		string out_chro_filename, vector<MSSpectrum<Peak1D> >& target_spectrum,
		vector<MSSpectrum<Peak1D> >& decoy_spectrum) {
	int exp_size = input_chro_filename.size();
	vector<fstream> fi(exp_size);
	for (int i = 0; i < exp_size; i++)
		fi[i].open(input_chro_filename[i].c_str(), fstream::binary | fstream::in);

	fstream fo(out_chro_filename.c_str(), fstream::binary | fstream::out | fstream::app);

	int trans_num = 0;
#pragma omp parallel
	{
		while (trans_num < size()) {
			TransitionGroup* t;
			bool is_exist = true;
#pragma omp critical
			{
				if (trans_num < size()) {
					t = trans[trans_num];
					trans_num++;
					for (int cur_exp_num = 0; cur_exp_num < exp_size; cur_exp_num++) {
						if (!t->load_store(fi[cur_exp_num], LoadSave::IOType::Read,
								LoadSave::Type::Precursor, cur_exp_num)) {
							throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
									"Error in load precursor ions when select_ms2, trans_num: "
											+ String(trans_num) + ", size: " + String(size()));
						}
						if (!t->load_store(fi[cur_exp_num], LoadSave::IOType::Read,
								LoadSave::Type::Product, cur_exp_num)) {
							throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
									"Error in load product ions when select_ms2, trans_num: "
											+ String(trans_num) + ", size: " + String(size()));
						}
					}
				} else
					is_exist = false;
			}

			if (is_exist) {
				AnalysisSimilar analysis(para);
				MSSpectrum<Peak1D> target_spec, decoy_spec;
				analysis.cal_similar_by_kmeans(*t, target_spec, decoy_spec);

#pragma omp critical
				{
					if (target_spec.size() >= 5) {
						target_spectrum.push_back(target_spec);

						if (decoy_spec.size() >= 5)
							decoy_spectrum.push_back(decoy_spec);
						//Storage
						t->load_store(fo, LoadSave::IOType::Write, LoadSave::Type::Basic_Precursor);
					}
				}
				delete t;
			}
		}
	}
	for (int i = 0; i < exp_size; i++)
		fi[i].close();
	fo.close();
}
