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

#include <ToolsTabFileWrite.h>

namespace GroupDIA {

} /* namespace GroupDIA */

GroupDIA::ToolsTabFileWrite::ToolsTabFileWrite(string filename) {
	fo.open(filename.c_str(), fstream::out);
	is_head_out = true;
}

GroupDIA::ToolsTabFileWrite::~ToolsTabFileWrite() {
	fo.close();
}

void GroupDIA::ToolsTabFileWrite::end_add_item() {
	if (names.empty())
		return;
	if (is_head_out) {
		for (int i = 0; i < names.size() - 1; i++) {
			fo << names[i] << "\t";
		}
		fo << *(names.rbegin()) << "\n";
		is_head_out = false;
	}

	for (int i = 0; i < values.size() - 1; i++) {
		fo << values[i] << "\t";
	}
	fo << *(values.rbegin()) << "\n";

	for (auto& it : values)
		it = "";
}

int GroupDIA::ToolsTabFileWrite::get_col_no(const string& name) {
	auto it = item_name_place.find(name);
	if (it != item_name_place.end()) {
		return it->second;
	} else {
		item_name_place.insert(pair<string, int>(name, values.size()));
		this->names.push_back(name);
		values.push_back("");
		return values.size() - 1;
	}
}
