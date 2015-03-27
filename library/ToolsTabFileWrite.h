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

#ifndef TOOLSTABFILEWRITE_H_
#define TOOLSTABFILEWRITE_H_

#include "ToolsSpectrum.h"
namespace GroupDIA {

class ToolsTabFileWrite {
public:
	ToolsTabFileWrite(string filename);
	~ToolsTabFileWrite();

	void end_add_item();

	template<class T> void inline add_item(const string& item_name, const T& item_value) {
		int cur_col_no = get_col_no(item_name);
		values[cur_col_no] = String(item_value);
	}
private:
	int get_col_no(const string& name);

	fstream fo;
	map<string, int> item_name_place;
	vector<string> names;
	vector<string> values;
	bool is_head_out;
};

} /* namespace GroupDIA */

#endif /* TOOLSTABFILEWRITE_H_ */
