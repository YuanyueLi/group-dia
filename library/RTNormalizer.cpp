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

#include "RTNormalizer.h"

namespace GroupDIA {
void RTNormalizer::store(string rt_filename) const {
	fstream fo(rt_filename, fstream::out);
	fo << anchor.size() << endl;
//	LoadSave::io_write(fo, anchor.size());

	for (auto it = anchor.begin(); it != anchor.end(); ++it) {
		fo << it->first << "\t" << it->second << endl;
//		LoadSave::io_write(fo, it->first);
//		LoadSave::io_write(fo, it->second);
	}

//	fo << dist.size() << "\t" << dist[0].size() << endl;
//	for (int i = 0; i < dist.size(); i++) {
//		for (int j = 0; j < dist.size(); j++) {
//			fo << dist[i][j] << "\t";
//		}
//		fo << endl;
//	}
}

void RTNormalizer::load(string rt_filename) {
	if (!File::exists(rt_filename))
		throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, rt_filename);

	fstream fi(rt_filename, fstream::in);
	int num;
	fi >> num;
//	LoadSave::io(fi, num, LoadSave::IOType::Read);
	for (int i = 0; i < num; i++) {
		double a, b;
//		LoadSave::io(fi, a, LoadSave::IOType::Read);
//		LoadSave::io(fi, b, LoadSave::IOType::Read);
		fi >> a;
		fi >> b;
		anchor.insert(pair<double, double>(a, b));
	}

//	double a, b;
//	fi >> a;
//	fi >> b;
//	dist.resize(a);
//	for (int i = 0; i < a; i++) {
//		dist[i] = vector<double>(b);
//		for (int j = 0; j < b; j++) {
//			fi >> dist[i][j];
//		}
//	}
}

double RTNormalizer::get_sample_rt(double reference_rt) const {
	if (reference_rt < anchor.begin()->first)
		return anchor.begin()->second;

	for (auto it = ++anchor.begin(); it != anchor.end(); ++it) {
		if (reference_rt < it->first) {
			auto ipre = it--;
			if (ipre == anchor.begin())
				return it->second;
			return OpenMS::Math::intervalTransformation(reference_rt, it->first, ipre->first,
					it->second, ipre->second);
		}
	}

	return anchor.rbegin()->second;
}

vector<double> RTNormalizer::align(const Ms1Intensity& exp_ref, const Ms1Intensity& exp_sam) {
	//Init
	int ref_size = exp_ref.size();
	int sam_size = exp_sam.size();

	//Build distance table
	dist = vector<vector<double> >(ref_size);
	for (auto it = dist.begin(); it != dist.end(); ++it)
		*it = vector<double>(sam_size);
#pragma omp parallel for schedule(dynamic) collapse(2)
	for (int i = 0; i < ref_size; i++)
		for (int j = 0; j < sam_size; j++) {
			dist[i][j] = get_dist(exp_ref, i, exp_sam, j);
		}

	vector<vector<double> > score(ref_size);
	vector<vector<int> > path(ref_size);
	for (auto it = path.begin(); it != path.end(); ++it)
		*it = vector<int>(sam_size);
	for (auto it = score.begin(); it != score.end(); ++it)
		*it = vector<double>(sam_size);

//Build score table
//0: i-1	1: j-1	2: i-1,j-1
	score[0][0] = 0;
	for (int i = 1; i < ref_size; i++) {
		for (int j = 1; j < sam_size; j++) {
			double d0 = score[i - 1][j] + dist[i - 1][j];
			double d1 = score[i][j - 1] + dist[i][j - 1];
			double d2 = score[i - 1][j - 1] + dist[i - 1][j - 1] * 2;
			if (d0 >= d1 and d0 > d2 and path[i - 1][j] != 0 and i != 1) {
				score[i][j] = d0;
				path[i][j] = 0;
				continue;
			} else if (d1 > d2 and path[i][j - 1] != 1 and j != 1) {
				score[i][j] = d1;
				path[i][j] = 1;
			} else {
				score[i][j] = d2;
				path[i][j] = 2;
			}
		}
	}

	vector<double> normalized_time_sam;
	normalized_time_sam = get_best_point(exp_ref, exp_sam, path);
	return normalized_time_sam;
}

vector<double> RTNormalizer::find_anchor(const Ms1Intensity& exp_ref, const Ms1Intensity& exp_sam,
		vector<vector<int> >& path) {
	//Find anchor
	vector<double> anchor_ref;
	vector<double> anchor_sam;

	int ref_size = path.size();
	int sam_size = path[0].size();
	int i = ref_size - 1, j = sam_size - 1;
	anchor_ref.push_back(exp_ref.get_origin_time(i));
	anchor_sam.push_back(j);
//	cout << "\tfind_anchor 1\n";
	while (i >= 1 and j >= 1) {
		if (path[i][j] == 0) {
			if (path[i - 1][j] == 1) {
				anchor_ref.push_back(
						(exp_ref.get_origin_time(i) + exp_ref.get_origin_time(i - 1)) / 2);
				anchor_sam.push_back((j + j - 1) / 2);
				i--, j--;
			} else {
				i--;
			}
		} else if (path[i][j] == 1) {
			if (path[i][j - 1] == 0) {
				anchor_ref.push_back(
						(exp_ref.get_origin_time(i) + exp_ref.get_origin_time(i - 1)) / 2);
				anchor_sam.push_back((j + j - 1) / 2);
				i--, j--;
			} else {
				j--;
			}
		} else {
			anchor_ref.push_back((exp_ref.get_origin_time(i) + exp_ref.get_origin_time(i - 1)) / 2);
			anchor_sam.push_back((j + j - 1) / 2);
			i--, j--;
		}
	}
	anchor_ref.push_back(exp_ref.get_origin_time(0));
	anchor_sam.push_back(0);
//	cout << "\tfind_anchor 2\n";
	//Build path
	vector<double> normalized_time = vector<double>(sam_size, -1);
	if (anchor_ref.empty() or anchor_sam.empty()) {
		Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Error in find_anchor, anchor is empty");
	}
	auto ian_i = anchor_ref.begin() + 1;
	auto ian_j = anchor_sam.begin() + 1;
	j = sam_size - 1;
	while (j >= 1 and ian_i != anchor_ref.end() and ian_j != anchor_sam.end()) {
		if (j < *ian_j) {
			++ian_i;
			++ian_j;
			continue;
		} else if (j >= *ian_j) {
			auto ian_i_next = ian_i - 1;
			auto ian_j_next = ian_j - 1;
			normalized_time[j] = (j - *ian_j) * (*ian_i - *ian_i_next) / (*ian_j - *ian_j_next)
					+ *ian_i;
			this->anchor[normalized_time[j]] = exp_sam.get_origin_time(j);
			j--;
			continue;
		}
	}
//	cout << "\tfind_anchor 3\n";
	normalized_time[0] = exp_ref.get_origin_time(0);
	this->anchor[normalized_time[0]] = exp_sam.get_origin_time(0);
	return normalized_time;
}

vector<double> RTNormalizer::get_best_point(const Ms1Intensity& exp_ref,
		const Ms1Intensity& exp_sam, vector<vector<int> >& path) {
	vector<double> point_ref;
	vector<double> point_sam;

	int ref_size = path.size();
	int sam_size = path[0].size();
	int i = ref_size - 1, j = sam_size - 1;
	point_ref.push_back(exp_ref.get_origin_time(i));
	point_sam.push_back(exp_sam.get_origin_time(j));

	while (i >= 1 and j >= 1) {
		point_ref.push_back((exp_ref.get_origin_time(i)));
		point_sam.push_back((exp_sam.get_origin_time(j)));
		if (path[i][j] == 0) {
			if (path[i - 1][j] == 1) {
				i--, j--;
			} else {
				i--;
			}
		} else if (path[i][j] == 1) {
			if (path[i][j - 1] == 0) {
				i--, j--;
			} else {
				j--;
			}
		} else {
			i--, j--;
		}
	}
	point_ref.push_back(exp_ref.get_origin_time(0));
	point_sam.push_back(exp_sam.get_origin_time(0));

	vector<double> normalized_time;
//	cout << point_ref.size() << "\t" << point_sam.size() << "\t" << normalized_time.size() << endl;
	Param lowpar;
	lowpar.setValue("window_size", 30);
	LowessSmoothing smooth;
	smooth.setParameters(lowpar);
	smooth.smoothData(point_sam, point_ref, normalized_time);
//	cout << point_ref.size() << "\t" << point_sam.size() << "\t" << normalized_time.size() << endl;

	for (int i = 0; i < point_ref.size(); i++) {
//		cout << point_ref[i] << "\t" << point_sam[i] << "\t" << normalized_time[i] << endl;
		this->anchor.insert(make_pair(normalized_time[i], point_sam[i]));
	}

	return normalized_time;
}

void Ms1Intensity::load_mzml(const Param& para, string ms1_mzml_filename) {
	//Load mzml file
	MzMLFile mzml;
	MapType exp;
	mzml.load(ms1_mzml_filename, exp);

	//Get information
	double mz_start = (double) para.getValue("mz_start");
	double mz_end = (double) para.getValue("mz_end");
	double mz_space = (double) para.getValue("mz_space");

	if (mz_start >= mz_end or mz_space <= 0 or mz_end < 0) {
		cout << mz_start << "\t" << mz_end << "\t" << mz_space << endl;
		throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Error: parameter error", "");
	}
	int mz_no = ceil((mz_end - mz_start) / mz_space);

	int exp_size = exp.size();

	origin_time.resize(exp_size);
	sum.resize(exp_size, 0);
	intensity.resize(exp_size);
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < exp_size; i++) {
		origin_time[i] = exp[i].getRT();
		auto& cur = intensity[i];
		cur.resize(mz_no, 0);
		auto ic = cur.begin();
		auto ie = exp[i].begin();
		double s = mz_start, e = s + mz_space;
		while (ic != cur.end() and ie != exp[i].end()) {
			if (ie->getMZ() < s) {
				++ie;
			} else if (ie->getMZ() >= e) {
				++ic;
				s = e;
				e += mz_space;
			} else {
				*ic += ie->getIntensity();
				++ie;
			}
		}

		ic = cur.begin();
		double avg = std::accumulate(ic, cur.end(), 0.0) / mz_no;
		for (auto ic = cur.begin(); ic != cur.end(); ++ic) {
			double temp_a = *ic - avg;
			sum[i] += (temp_a * temp_a);
		}
	}
}

double RTNormalizer::get_dist(const Ms1Intensity& exp_ref, int time_ref,
		const Ms1Intensity& exp_sam, int time_sam) {
	if (exp_ref.sum[time_ref] == 0 or exp_sam.sum[time_sam] == 0)
		return 0;

	auto iref = exp_ref.intensity[time_ref].begin(), isam = exp_sam.intensity[time_sam].begin();

	double numerator = 0;
	for (; iref != exp_ref.intensity[time_ref].end(); ++iref, ++isam)
		numerator += (*iref) * (*isam);
	return numerator / sqrt(exp_ref.sum[time_ref] * exp_sam.sum[time_sam]);
}

int Ms1Intensity::size() const {
	return intensity.size();
}

double Ms1Intensity::get_origin_time(int i) const {
	return origin_time[i];
}

}/* namespace GroupDIA */
