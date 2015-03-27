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

#ifndef TRANSITIONGROUP_H_
#define TRANSITIONGROUP_H_
#include "RTNormalizer.h"
#include "ToolsLoadSave.h"
#include "SwathExperiment.h"
#include "Precursor.h"
#include "Product.h"
#include "ToolsSpectrumAnalysis.h"
#include "ToolsDeIsotoper.h"
#include <atomic>

namespace GroupDIA {
class TransitionGroup {
	friend class TestTransition;
public:
	bool set_transition(const SwathExperiment& exp, const Feature& feature, int exp_num, int scan_num, int exact_scan_num);

	void add_ms1(const SwathExperiment& exp, const RTNormalizer& rt, int exp_num, int exact_scan_num);

	void set_product_ions_mz(const SwathExperiment& exp, double ms1_win_start, double ms1_win_end);
	void add_ms1_nearby(const SwathExperiment& exp, int exp_num, int exact_scan_num);
	void add_ms2_nearby(const SwathExperiment& exp, int exp_num);

	void add_ms2_nearby_in_record_mode(const SwathExperiment& exp, int exp_num);
	void record_ms2_nearby(SwathExperiment& exp, int exp_num);
	void record_feature_product_peak(SwathExperiment& exp);
	void add_feature_product(const SwathExperiment& exp);

	void crude_refine(double MIN_ALLOWED_CORRELATION);
	void refined_by_crosscorrelation(int cur_exp_num, int MAX_ALLOWED_DELAY, int FINE_NEARBY, double MIN_INTENSITY);

	void remove_nearby_intensity(int exp_num);
	void remove_feature_product();
	void remove_ions();

	void refine_product_ions(const vector<int>& remained_place);

	void set_pep_seq(const string& seq);
	void set_pep_seq(const AASequence& seq);
	void set_uid(UInt64 uid);
	void set_scan_num(int scan_num);

	int size() const;
	int get_scan_num() const;
	UInt64 get_uid() const;
	int get_feature_exp_num() const; //return the feature_precursor exp num

	double get_precursor_rt(int exp_num) const;

	int get_exp_size() const;

	double get_precursor_mz() const;
	int get_precursor_charge() const;
	void get_precursor_nearby_intensity(int exp_num, vector<double>& intensity) const;
	void get_precursor_peak_intensity(vector<double>& intensity) const;
	void get_precursor_rt_in_nearby_no(int exp_num, int& peak_start_no, int& peak_length) const;
	void get_precursor_rt_no(int exp_num, int& peak_start_no, int& peak_length) const;
	void set_precursor_rt_no(int exp_num, int& peak_start_no, int& peak_length);

	int get_product_rt_size(int exp_num, double rt_start, double rt_end) const;
	int get_product_rt_no(int exp_num, double rt) const;
	void get_product_mz(vector<double>& mz) const;
	int get_product_num() const;

	void get_product_intensity_in_each_time(int cur_exp_num, vector<vector<double> >& intensity) const;
	void get_product_intensity_in_each_time_nearby_version(vector<vector<double> >& intensity) const;

	virtual bool load_store(fstream& f, int io_type, int type = LoadSave::Type::Null, int ion_num = 0);

	virtual ~TransitionGroup() {
	}
protected:
	int scan_num;
	UInt64 feature_id;
	int feature_exp_num;
	double feature_score;
	int precursor_charge;
	double precursor_mz;
	int decoy;

//	string sequence;
	AASequence aasequence;
	string accession;
	double identification_score;

	PrecursorIon feature_precursor;
	ProductIons feature_product;

	vector<PrecursorIon> precursor_ions;
	vector<ProductIons> product_ions;
	vector<double> product_mz;

};

} /* namespace GroupDIA */
#endif /* TRANSITIONGROUP_H_ */
