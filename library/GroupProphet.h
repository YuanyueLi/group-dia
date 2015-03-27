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


#ifndef GROUPPROPHET_H_
#define GROUPPROPHET_H_

#include "ToolsSpectrum.h"
#include "ToolsTabFile.h"
#include <linalg.h>
#include <dataanalysis.h>
#include <boost/math/distributions/normal.hpp>
#include <random>
struct Item {
	int group_id;
	int id;
	int type;
	double score;
	double p;
	double ide_p;
	bool is_in_group;
	int feature_no;

	vector<double> values;
	vector<double> intensity;
};

namespace GroupDIA {

class GroupProphet: public DefaultParamHandler {
public:
	enum ClassType {
		Target = 1, Decoy = 0, Unknown = 2
	};

	GroupProphet();

	void updateMembers_();

	int add_values(int group_id, int id, const vector<double>& values,
			const vector<double>& intensity, int class_type, double ide_p, bool is_in_group);

	bool get_final_result(int id, int& no, double& score, double& p, double& ide_score) const;

	bool learn();

	void init();

	vector<int> get_all_id() const;
	vector<int> get_all_gid() const;

private:
	multimap<int, int> map_id_no;
	multimap<int, int> map_gid_id;
	vector<Item> candidates;
	vector<double> final_score;

	double fraction;
	double lambda;
	double it_p;
	double start_p;
	int xval_num;
	int iter_num;

	void learn(int learn_num, vector<double>& para, vector<int>& selected_target_no) const;
	void learn(const vector<int>& target_no, const vector<int>& decoy_no,
			vector<double>& para) const;

	void get_all_gids(vector<int>& target_gid, vector<int>& decoy_gid) const;

	double get_max_score(int id, const vector<double>& all_score, int& no) const;

	vector<double> get_score_by_gid(const vector<int>& gid, const vector<double>& score) const;

	vector<int> get_no_by_id(const vector<int>& id, const vector<double>& score) const;
	vector<int> get_id_by_gid(int gid) const;

	void get_train_data_by_fraction(double fraction, vector<int>& target_gid,
			vector<int>& decoy_gid) const;

	void select_target_first(const vector<int>&target_gid, vector<int>&selected_target_no) const;

	void score_by_para(const vector<double>& para, vector<double>& score) const;

	void select_target_by_p(vector<int> target_gid, vector<int> decoy_gid,
			const vector<int>& decoy_no, const vector<double>& score, vector<double>& ori_p,
			vector<int>&selected_target_no) const;

	vector<double> calculate_p(const vector<int>& decoy_no, const vector<int>& selected_target_no,
			const vector<double>& score, vector<double>& ori_p) const;

	void estimate_p_ts_existspec(const vector<int>& selected_target_gid,
			const vector<double>& p) const;
	double get_p_ts(vector<double> a, vector<double> b) const;

	double min_p_ts;
	double max_p_ts;
};

} /* namespace GroupDIA */

#endif /* GROUPPROPHET_H_ */
