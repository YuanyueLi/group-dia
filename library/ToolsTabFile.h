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

#ifndef TOOLSTABFILE_H_
#define TOOLSTABFILE_H_

#include "ToolsSpectrum.h"
namespace GroupDIA {
class ToolsTabFile {
public:
	void store(string filename) const;
	void load(string filename);

	int get_size() const;

	string get_item(int no, string name) const;

	template<class T> void inline add_item(const string& name, const T& values) {
		int cur_col_no = get_col_no(name);
		*(item_values[cur_col_no].rbegin()) = String(values);
//		item_values[cur_col_no].push_back(String(values));
	}

	template<class T> void inline add_item(const string&place_item, const string& place_value, const string& add_name,
			const T& add_values) {
		int place_col_no = get_col_no(place_item);
		int place_no = 0;
		for (; place_no < item_values[place_col_no].size(); ++place_no) {
			if (item_values[place_col_no][place_no].compare(place_value) == 0)
				break;
		}

		if (place_no == item_values[place_col_no].size()) {
			add_item(add_name, add_values);
			return;
		} else {
			int cur_col_no = get_col_no(add_name);
			item_values[cur_col_no][place_no] = String(add_values);
		}
	}

	template<class T> void inline set_item(int no, const string& name, const T& values) {
		int cur_col_no = get_col_no(name);
		item_values[cur_col_no][no] = String(values);
	}

	void add_new_item();
	void remove_item(int no);
	void select_item(vector<int> remain_no);

	void sort_by_column(const string& name);
private:
	int get_col_no(const string& name);

	vector<string> item_name;
	map<string, int> item_name_place;
	vector<vector<string> > item_values; //col (item_name), row
};
} /* namespace GroupDIA */
#endif /* TOOLSTABFILE_H_ */
