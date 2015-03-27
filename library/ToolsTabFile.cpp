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

#include "ToolsTabFile.h"

namespace GroupDIA {

void ToolsTabFile::load(string filename) {
	fstream fi(filename, fstream::in);
	String line;
	getline(fi, line);
	vector<String> items;
	line.split("\t", items);
	for (int i = 0; i < items.size(); i++) {
		item_name.push_back(items[i]);
		item_name_place[items[i]] = i;
	}
	item_values.resize(item_name.size());

	getline(fi, line);
	while (!line.empty()) {
		vector<String> c;
		line.split("\t", c);
		for (int i = 0; i < item_values.size(); i++) {
			if (i < c.size())
				item_values[i].push_back(c[i]);
			else
				item_values[i].push_back("");
		}

		getline(fi, line);
	}
}

void ToolsTabFile::store(string filename) const {
	int name_num = item_name.size();
	int list_num = 0;
	if (name_num == 0) {
		return;
	} else {
		list_num = item_values[0].size();
	}

	for (auto& it : item_values)
		if (it.size() != list_num)
			cout << "Error in store" << endl;

	fstream fo;
	fo.open(filename.c_str(), fstream::out);
	for (int i = 0; i < name_num; i++) {
		fo << item_name[i] << "\t";
	}
	fo << endl;

	for (int i = 0; i < list_num; i++) {
		for (int j = 0; j < name_num; j++) {
			fo << item_values[j][i] << "\t";
		}
		fo << endl;
	}
	fo.close();
}

int ToolsTabFile::get_col_no(const string& name) {
	auto it = item_name_place.find(name);
	if (it != item_name_place.end()) {
		return it->second;
	} else {
		int col_no = item_name.size();
		item_name.push_back(name);
		item_name_place[name] = col_no;

		int max_col_size = 0;
		if (col_no > 0) {
			max_col_size = item_values[col_no - 1].size();
			for (int i = 0; i < col_no - 1; i++)
				if (max_col_size < item_values[i].size())
					max_col_size = item_values[i].size();
		}
		if (max_col_size == 0)
			max_col_size++;

		item_values.push_back(vector<string>(max_col_size, ""));
		return col_no;
	}
}

int ToolsTabFile::get_size() const {
	if (item_name.empty())
		return 0;
	return item_values[0].size();
}

string ToolsTabFile::get_item(int no, string name) const {
	auto in = item_name_place.find(name);
	if (in != item_name_place.end())
		return item_values.at(in->second).at(no);
	else {
		cout << "Error: Not find: " << name << endl;
		return "";
	}
}
} /* namespace GroupDIA */

void GroupDIA::ToolsTabFile::add_new_item() {
	int max_item_num = 0;
	for (auto& it : item_values) {
		if (max_item_num < it.size())
			max_item_num = it.size();
	}

	for (auto &it : item_values) {
		if (max_item_num > it.size()) {
			it.resize(max_item_num, "");
		}
	}

	for (auto &it : item_values)
		it.push_back("");
}

template<class T1, class T2> bool large_first(const pair<T1, T2>& i, const pair<T1, T2>& j) {
	return (i.first > j.first);
}

void GroupDIA::ToolsTabFile::sort_by_column(const string& name) {
	int col_no;
	auto it = item_name_place.find(name);
	if (it != item_name_place.end()) {
		col_no = it->second;
	} else {
		cout << "Error in sort, not find " << name << endl;
		return;
	}

	vector<pair<string, int> > temp;
	int no = 0;
	const auto& item = item_values[col_no];
	for (int i = 0; i < item.size(); i++) {
		temp.push_back(pair<string, int>(item[i], no));
		no++;
	}
	sort(temp.begin(), temp.end(), large_first<string, int>);

	for (int i = 0; i < item_values.size(); ++i) {
		vector<string> new_values = item_values[i];
		for (int j = 0; j < temp.size(); j++) {
			item_values[i][j] = new_values[temp[j].second];
		}
	}
}

void GroupDIA::ToolsTabFile::remove_item(int no) {
	if (no >= get_size())
		return;
	for (int i = 0; i < item_values.size(); i++) {
		item_values[i].erase(item_values[i].begin() + no);
	}
}

void GroupDIA::ToolsTabFile::select_item(vector<int> remain_no) {
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < item_values.size(); i++) {
		vector<string> new_values;
		new_values.reserve(item_values[i].size());
		for (auto j : remain_no) {
			new_values.push_back(item_values[i][j]);
		}
		item_values[i].swap(new_values);
	}
}
