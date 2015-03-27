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


#include "GroupProphet.h"

const double p_ts_in_existspec = 0.7;
const double p_ts_out_existspec = 0.3;

namespace GroupDIA {

} /* namespace GroupDIA */

template<class T> bool large_first_single(const T& i, const T& j) {
	return (i > j);
}

template<class T1, class T2> bool large_first(const pair<T1, T2>& i, const pair<T1, T2>& j) {
	return (i.first > j.first);
}

template<class T1, class T2> bool small_first(const pair<T1, T2>& i, const pair<T1, T2>& j) {
	return (i.first < j.first);
}

GroupDIA::GroupProphet::GroupProphet() :
		DefaultParamHandler("mProhet") {
	defaults_.setValue("fraction", 0.6, "fraction for mProhet");
	defaults_.setValue("lambda", 0.4, "lambda for mProhet");
	defaults_.setValue("it_p", 0.99, "it_fdr for mProhet");
	defaults_.setValue("start_p", 0.8, "start_p for mProhet");
	defaults_.setValue("iter_num", 10, "iter_num for mProhet");
	defaults_.setValue("xval_num", 5, "xval_num for mProhet");
	min_p_ts = 0.2;
	max_p_ts = 0.8;
	// write defaults into Param object param_
	defaultsToParam_();
}

void GroupDIA::GroupProphet::updateMembers_() {
	fraction = (DoubleReal) param_.getValue("fraction");
	lambda = (DoubleReal) param_.getValue("lambda");
	it_p = (DoubleReal) param_.getValue("it_p");
	start_p = (DoubleReal) param_.getValue("start_p");
	iter_num = param_.getValue("iter_num");
	xval_num = param_.getValue("xval_num");
}

void GroupDIA::GroupProphet::learn(int learn_num, vector<double>& para,
		vector<int>& selected_target_no) const {
//	learn_num = 0;
	if (candidates.empty())
		cout << "Error: Not find values" << endl;
//Select for xval
	vector<int> target_gid, decoy_gid;
	get_train_data_by_fraction(fraction, target_gid, decoy_gid);

//Choose target, E-step
//	vector<int> selected_target_no;
	select_target_first(target_gid, selected_target_no);
	vector<int> decoy_no;
	for (auto gt = decoy_gid.begin(); gt != decoy_gid.end(); ++gt) {
		vector<int> ids = get_id_by_gid(*gt);
		for (auto it = ids.begin(); it != ids.end(); ++it) {
			auto it_start = map_id_no.lower_bound(*it);
			auto it_end = map_id_no.upper_bound(*it);
			for (auto jt = it_start; jt != it_end; ++jt) {
				decoy_no.push_back(jt->second);
			}
		}
	}

//	{
//		cout << "--------------------------------------------------" << endl;
//		cout << target_gid.size() << "\t" << selected_target_no.size() << "\t" << decoy_no.size()
//				<< "\t" << endl;
//	}
	if (selected_target_no.size() < 3 or decoy_no.size() < 3) {
		para.clear();
		return;
	}

	learn(selected_target_no, decoy_no, para);

	vector<double> score, ori_p;
	score_by_para(para, score);

//Choose type, M-step
	for (int cur_repeat = 0; cur_repeat < learn_num; cur_repeat++) {
		if (para.empty()) {
			para.clear();
			return;
		}
		select_target_by_p(target_gid, decoy_gid, decoy_no, score, ori_p, selected_target_no);

//		{
//			double min = 1000;
//			for (int i = 0; i < selected_target_no.size(); i++)
//				if (score[selected_target_no[i]] < min)
//					min = score[selected_target_no[i]];
//			cout << selected_target_no.size() << "\t" << decoy_no.size() << "\t" << endl;
//			cout << endl;
//		}

		if (selected_target_no.size() < 3 or decoy_no.size() < 3) {
			para.clear();
			return;
		}

		learn(selected_target_no, decoy_no, para);
		score_by_para(para, score);
	}
	select_target_by_p(target_gid, decoy_gid, decoy_no, score, ori_p, selected_target_no);
}

vector<int> GroupDIA::GroupProphet::get_no_by_id(const vector<int>& id,
		const vector<double>& score) const {
	vector<int> result;
	result.reserve(id.size());
	for (auto it : id) {
		int cur_no;
		get_max_score(it, score, cur_no);
		result.push_back(cur_no);
	}
	return result;
}

vector<double> GroupDIA::GroupProphet::get_score_by_gid(const vector<int>& gid,
		const vector<double>& score) const {
	vector<double> result;
	for (auto g : gid) {
		vector<int> ids = get_id_by_gid(g);
		for (auto it : ids) {
			int cur_no;
			double cur_score = get_max_score(it, score, cur_no);
			result.push_back(cur_score);
		}
	}
	return result;
}

vector<int> GroupDIA::GroupProphet::get_all_id() const {
	if (map_id_no.empty())
		cout << "Error in get_all_id" << endl;
	vector<int> result;
	auto it = map_id_no.begin();
	int cur_id = it->first;
	result.push_back(cur_id);
	for (; it != map_id_no.end(); ++it) {
		if (cur_id != it->first) {
			cur_id = it->first;
			result.push_back(cur_id);
		}
	}
	return result;
}

vector<int> GroupDIA::GroupProphet::get_all_gid() const {
	if (map_gid_id.empty())
		cout << "Error in get_all_id" << endl;
	vector<int> result;
	auto it = map_gid_id.begin();
	int cur_id = it->first;
	result.push_back(cur_id);
	for (; it != map_gid_id.end(); ++it) {
		if (cur_id != it->first) {
			cur_id = it->first;
			result.push_back(cur_id);
		}
	}
	return result;
}

void GroupDIA::GroupProphet::select_target_first(const vector<int>& target_gid,
		vector<int>& selected_target_no) const {
	selected_target_no.clear();
	for (auto gid = target_gid.begin(); gid != target_gid.end(); ++gid) {
		vector<int> all_ids = get_id_by_gid(*gid);
		for (auto it = all_ids.begin(); it != all_ids.end(); ++it) {
			auto it_start = map_id_no.lower_bound(*it);
			auto it_end = map_id_no.upper_bound(*it);
			for (auto jt = it_start; jt != it_end; ++jt) {
				if (candidates[jt->second].type == Target) {
					selected_target_no.push_back(jt->second);
				}
			}
		}
	}
}

void GroupDIA::GroupProphet::select_target_by_p(vector<int> target_gid, vector<int> decoy_gid,
		const vector<int>& decoy_no, const vector<double>& score, vector<double>& ori_p,
		vector<int>& selected_target_no) const {
	vector<double> p = calculate_p(decoy_no, selected_target_no, score, ori_p);
	vector<int> all_gid = target_gid;
	all_gid.insert(all_gid.end(), decoy_gid.begin(), decoy_gid.end());
//	estimate_p_ts_existspec(all_gid, p, parameter);

//Choose target by new p value
	selected_target_no.clear();
	for (auto gid = target_gid.begin(); gid != target_gid.end(); ++gid) {
		vector<int> ids = get_id_by_gid(*gid);
		for (auto id = ids.begin(); id != ids.end(); ++id) {
			int no;
			double cur_p = get_max_score(*id, p, no);
			if (cur_p >= it_p) {
				selected_target_no.push_back(no);
			}
		}
	}
}

void GroupDIA::GroupProphet::learn(const vector<int>& target_no, const vector<int>& decoy_no,
		vector<double>& para) const {
//	vector<int> target_no = get_no_by_id(target_id, score);
//	vector<int> decoy_no = get_no_by_id(decoy_id, score);

	if (target_no.empty() or decoy_no.empty()) {
		cout << "Error in learn: target or decoy no is zero" << endl;
		para.clear();
		return;
	}

	using namespace alglib;
	para.clear();
	ae_int_t npoints = target_no.size() + decoy_no.size();
	ae_int_t nvars = candidates[0].values.size();
	ae_int_t nclasses = 2;
	real_2d_array xy;
	xy.setlength(npoints, nvars + 1);
	int values_var = candidates[0].values.size();
	for (int i = 0; i < target_no.size(); i++) {
		int no = target_no[i];
		for (int j = 0; j < values_var; j++) {
			xy[i][j] = candidates[no].values[j];
		}
		xy[i][nvars] = ClassType::Target;
	}
	int target_size = target_no.size();
	for (int i = 0; i < decoy_no.size(); i++) {
		int no = decoy_no[i];
		for (int j = 0; j < values_var; j++) {
			xy[i + target_size][j] = candidates[no].values[j];
		}
		xy[i + target_size][nvars] = ClassType::Decoy;
	}

	ae_int_t info;
	real_1d_array w;

	fisherlda(xy, npoints, nvars, nclasses, info, w);

	if (info < 1 or w.length() != nvars) {
		cout << "Error in learning, id: " << info << ", size of w is: " << w.length() << endl;
		para.clear();
		return;
	}

//Decide which one is target
	double total_target = 0, total_decoy = 0;
	for (int i = 0; i < target_no.size(); i++) {
		double s = 0;
		for (int j = 0; j < nvars; j++) {
			s += xy[i][j] * w[j];
		}
		total_target += s;
	}
	for (int i = 0; i < decoy_no.size(); i++) {
		double s = 0;
		for (int j = 0; j < nvars; j++) {
			s += xy[i + target_size][j] * w[j];
		}
		total_decoy += s;
	}
	if (total_target / (target_size * 1.0) < total_decoy / (decoy_no.size() * 1.0)) {
		for (int i = 0; i < nvars; i++)
			w[i] = 0 - w[i];
	}

	for (int i = 0; i < nvars; i++) {
		para.push_back(w[i]);
	}
}

double GroupDIA::GroupProphet::get_max_score(int id, const vector<double>& all_score,
		int& no) const {
	auto it_start = map_id_no.lower_bound(id);
	auto it_end = map_id_no.upper_bound(id);

	if (it_start == it_end) {
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Error: not find id: " + String(id));
	}

	no = it_start->second;
	double _score = all_score[no];

	for (auto it = it_start; it != it_end; ++it) {
		int cur_no = it->second;
		double cur_score = all_score[cur_no];
		if (cur_score > _score) {
			_score = cur_score;
			no = cur_no;
		}
	}
	return _score;
}

void GroupDIA::GroupProphet::get_all_gids(vector<int>& target_gid, vector<int>& decoy_gid) const {
	int cur_gid = map_gid_id.begin()->first;
	int cur_id = map_id_no.find(cur_gid)->first;
	int cur_no = map_id_no.find(cur_gid)->second;

	if (candidates[cur_no].type == ClassType::Decoy)
		decoy_gid.push_back(cur_gid);
	else
		target_gid.push_back(cur_gid);

	for (auto it = map_gid_id.begin(); it != map_gid_id.end(); ++it) {
		if (it->first != cur_id) {
			auto id = map_id_no.find(it->second);
			if (candidates[id->second].type == ClassType::Decoy)
				decoy_gid.push_back(it->first);
			else
				target_gid.push_back(it->first);
			cur_id = it->first;
		}
	}

	if (target_gid.empty() or decoy_gid.empty())
		cout << "Target id is empty or decoy id is empty" << endl;
}

void GroupDIA::GroupProphet::get_train_data_by_fraction(double fraction, vector<int>& target_gid,
		vector<int>& decoy_gid) const {
	vector<int> all_target, all_decoy;
	get_all_gids(all_target, all_decoy);

	ToolsSpectrum::random(all_target.begin(), all_target.end());
	ToolsSpectrum::random(all_decoy.begin(), all_decoy.end());

//	target_gid = vector<int>(all_target.begin(), all_target.begin() + all_target.size() * fraction);

	target_gid.clear();
	for (auto gid = all_target.begin(); gid != all_target.begin() + all_target.size() * fraction;
			++gid) {
		auto id = map_gid_id.find(*gid);
		if (id != map_gid_id.end()) {
			target_gid.push_back(*gid);
		}
	}
	decoy_gid = vector<int>(all_decoy.begin(), all_decoy.begin() + all_decoy.size() * fraction);
}

int GroupDIA::GroupProphet::add_values(int group_id, int id, const vector<double>& values,
		const vector<double>& intensity, int class_type, double ide_p, bool is_in_group) {
	int no = candidates.size();
	map_id_no.insert(pair<int, int>(id, no));
	if (map_gid_id.find(group_id) == map_gid_id.end())
		map_gid_id.insert(make_pair(group_id, id));
	else {
		bool is_exist = false;
		for (auto it = map_gid_id.lower_bound(group_id); it != map_gid_id.upper_bound(group_id);
				++it) {
			if (it->second == id) {
				is_exist = true;
				break;
			}
		}
		if (!is_exist)
			map_gid_id.insert(make_pair(group_id, id));
	}

	Item can;
	can.group_id = group_id;
	can.id = id;
	can.values = values;
	can.type = class_type;
	can.ide_p = ide_p;
	can.is_in_group = is_in_group;
	can.intensity = intensity;

	candidates.push_back(can);

	return no;
}

bool GroupDIA::GroupProphet::get_final_result(int _id, int& no, double& score, double& p,
		double& ide_score) const {
	score = get_max_score(_id, final_score, no);
	p = candidates[no].p;
	ide_score = candidates[no].ide_p;

	return true;
}

bool GroupDIA::GroupProphet::learn() {
	init();
	cout << "Start learning" << endl;
	if (candidates.empty())
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Error in learn, value is empty");

	vector<double> para;
	set<int> selected_target_no;

	int num = 0;
	int total_num = 0;
#pragma omp parallel
	{
		while (num < xval_num and total_num < xval_num * 10) {
			vector<double> cur_para;
			vector<int> cur_selected_target_no;
			bool should_learn = false;
#pragma omp critical
			if (num < xval_num) {
				num++;
				should_learn = true;
			}

			if (should_learn) {
				learn(iter_num, cur_para, cur_selected_target_no);
#pragma omp critical
				{
					total_num++;
					if (!cur_para.empty()) {
						if (!para.empty()) {
							for (int i = 0; i < cur_para.size(); i++)
								para[i] += cur_para[i];
						} else {
							para = cur_para;
						}
						for (auto it : cur_selected_target_no)
							selected_target_no.insert(it);
					} else
						num--;
				}
			}
		}
	}

	if (total_num >= xval_num * 10)
		return false;

	if (para.empty()) {
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Error in learn, parameter is empty");
	}

	cout << "Start prepare output" << endl;
	for (auto& it : para)
		it /= (xval_num * 1.0);
	double xval_n = xval_num;
	vector<int> v_selected_target_no;
	for (auto it = selected_target_no.begin(); it != selected_target_no.end(); ++it)
		v_selected_target_no.push_back(*it);

	score_by_para(para, final_score);

	vector<int> target_gid, decoy_gid;
	get_all_gids(target_gid, decoy_gid);

	vector<double> decoy_score;
	decoy_score = get_score_by_gid(decoy_gid, final_score);
	double decoy_mean, decoy_sd;
	ToolsSpectrum::mean_and_std_dev(decoy_score, decoy_mean, decoy_sd);

	for (auto& it : final_score)
		it = (it - decoy_mean) / decoy_sd;

	decoy_score = get_score_by_gid(decoy_gid, final_score);

	vector<int> decoy_no;
	for (auto gid : decoy_gid) {
		vector<int> ids = get_id_by_gid(gid);
		for (auto id : ids) {
			int cur_no;
			get_max_score(id, final_score, cur_no);
			decoy_no.push_back(cur_no);
		}
	}
	vector<double> ori_p;
	vector<double> p = calculate_p(decoy_no, v_selected_target_no, final_score, ori_p);
	p = calculate_p(decoy_no, v_selected_target_no, final_score, p);

	for (int i = 0; i < candidates.size(); i++) {
		candidates[i].score = final_score[i];
		candidates[i].p = p[i];
	}

	cout << "End learning" << endl;
	return true;
}

void GroupDIA::GroupProphet::score_by_para(const vector<double>& para, vector<double>& score) const {
	if (para.empty()) {
		score.clear();
		return;
	}

	int size = para.size();
	score.resize(candidates.size());
	for (int i = 0; i < candidates.size(); i++) {
		score[i] = 0;
		for (int j = 0; j < size; j++) {
			score[i] += para[j] * candidates[i].values[j];
		}
	}

}

vector<double> GroupDIA::GroupProphet::calculate_p(const vector<int>& decoy_no,
		const vector<int>& selected_target_no, const vector<double>& score,
		vector<double>& ori_p) const {
	//Get decoy_score
	vector<double> decoy_score(decoy_no.size());
	for (int i = 0; i < decoy_no.size(); i++) {
		decoy_score[i] = score[decoy_no[i]];
	}

	using namespace boost::math;
	double decoy_mean, decoy_stdev;
	ToolsSpectrum::mean_and_std_dev(decoy_score, decoy_mean, decoy_stdev);
	boost::math::normal_distribution<> decoy_dist(decoy_mean, decoy_stdev);

	//Get target_score
	vector<double> target_score(selected_target_no.size());
	for (int i = 0; i < selected_target_no.size(); i++) {
		target_score[i] = score[selected_target_no[i]];
	}

	using namespace boost::math;
	double target_mean, target_stdev;
	ToolsSpectrum::mean_and_std_dev(target_score, target_mean, target_stdev);
	boost::math::normal_distribution<> target_dist(target_mean, target_stdev);

	vector<double> p = ori_p;
	p.resize(score.size());
	for (int i = 0; i < score.size(); i++) {
		double p_target = 0.5;
		if (!ori_p.empty()) {
			int feature_no = candidates[i].feature_no;
			if (i != feature_no) {
				if (feature_no >= 0) {
					if (ori_p[feature_no] > it_p) {
						p_target = get_p_ts(candidates[i].intensity,
								candidates[feature_no].intensity);
					} else {
						p_target = min_p_ts;
					}
				} else {
					p_target = min_p_ts;
				}
			} else {
				bool is_exist_group = false;
				vector<int> ids = get_id_by_gid(candidates[i].group_id);
				for (auto id = ids.begin(); id != ids.end(); ++id) {
					int no;
					double max_p = get_max_score(*id, ori_p, no);
					if (max_p > it_p and no != i) {
						is_exist_group = true;
						break;
					}
				}
				if (is_exist_group) {
					p_target = get_p_ts(candidates[i].intensity, candidates[feature_no].intensity);
				} else {
					p_target = 0.5;
				}
			}
		}

		double p_decoy = 1 - p_target;
		double p_score_decoy = boost::math::cdf(boost::math::complement(decoy_dist, score[i]));
		double p_score_target = boost::math::cdf(target_dist, score[i]);

		p[i] = (p_score_target * p_target) / (p_score_target * p_target + p_score_decoy * p_decoy);
	}
	return p;
}

vector<int> GroupDIA::GroupProphet::get_id_by_gid(int gid) const {
	auto is = map_gid_id.lower_bound(gid);
	auto ie = map_gid_id.upper_bound(gid);
	vector<int> result;
	for (auto it = is; it != ie; ++it) {
		result.push_back(it->second);
	}
	return result;
}

void GroupDIA::GroupProphet::estimate_p_ts_existspec(const vector<int>& selected_target_gid,
		const vector<double>& p) const {
	int n_ts_in = 0, n_fs_in = 0, n_ts_out = 0, n_fs_out = 0;

	for (auto gid = selected_target_gid.begin(); gid != selected_target_gid.end(); ++gid) {
		auto id = map_id_no.find(*gid);
		if (id != map_id_no.end()) {
			if (candidates[id->second].ide_p <= 0.1)
				continue;
		}

		vector<int> all_ids = get_id_by_gid(*gid);
		int cur_group_num = 0;

		for (auto id = all_ids.begin(); id != all_ids.end(); ++id) {
			int max_no;
			get_max_score(*id, p, max_no);

			auto it_start = map_id_no.lower_bound(*id);
			auto it_end = map_id_no.upper_bound(*id);

			if (p[max_no] >= it_p) {
				for (auto it = it_start; it != it_end; ++it) {
					if (max_no == it->second) {
						//Is selected as max and TrueSpec
						cur_group_num++;
						if (candidates[max_no].is_in_group)
							n_ts_in++;
						else
							n_ts_out++;
					} else {
						if (candidates[it->second].is_in_group)
							n_fs_in++;
						else
							n_fs_out++;
					}
				}
			}
		}
	}
}

void GroupDIA::GroupProphet::init() {
	vector<int> gids = get_all_gid();
	for (auto gid = gids.begin(); gid != gids.end(); ++gid) {
		vector<int> ids = get_id_by_gid(*gid);
		int feature_no = -1;
		for (auto id = ids.begin(); id != ids.end(); ++id) {
			auto ino_start = map_id_no.lower_bound(*id);
			auto ino_end = map_id_no.upper_bound(*id);
			for (auto ino = ino_start; ino != ino_end; ++ino) {
				if (candidates[ino->second].type == Target) {
					feature_no = ino->second;
					break;
				}
			}
			if (feature_no != -1)
				break;
		}
		for (auto id = ids.begin(); id != ids.end(); ++id) {
			auto ino_start = map_id_no.lower_bound(*id);
			auto ino_end = map_id_no.upper_bound(*id);
			for (auto ino = ino_start; ino != ino_end; ++ino) {
				candidates[ino->second].feature_no = feature_no;
			}
		}
	}
}

double GroupDIA::GroupProphet::get_p_ts(vector<double> a, vector<double> b) const {
	int size = a.size();
	double r = 0, ra = 0, rb = 0;
	for (int i = 0; i < size; i++) {
		r += a[i] * b[i];
		ra += a[i] * a[i];
		rb += b[i] * b[i];
	}
	if (r == 0)
		return min_p_ts;
	double sim = r / sqrt(ra * rb);
	return min_p_ts + (sim * (max_p_ts - min_p_ts));
}
