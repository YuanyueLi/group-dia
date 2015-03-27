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

#ifndef MERGETRANSITIONS_H_
#define MERGETRANSITIONS_H_

#include "TransitionIdentified.h"
#include "ScanPepXMLFile.h"

struct PeptideInformation {
	TransitionIdentified tr;
	vector<double> product_mz;
	vector<pair<int, int> > rt_no;
};

namespace GroupDIA {

class MergeTransitions {
public:
	MergeTransitions(Param para, Param para_ide, int windows_size, int cur_windows_num,int exp_size);
	void add_identification(String pepxml_filename);
	void choose_best_id();
	void load_transition(string input_filename, string output_filename);

private:
	Param para;
	Param para_ide;

	int exp_size;
	int windows_size;
	int cur_windows_num;
	int cur_id;

	map<int, PeptideHit> map_scan_pep; //TODO: Should check is unique!

	map<string, vector<pair<double, int> > > max_score;
	map<string, PeptideInformation> result_identification;
	set<int> best_scan_num;

	bool get_next_identification(fstream& fi,
			vector<pair<TransitionIdentified, int> >& result_trans);
	int get_exp_num_from_scan_num(int scan_num) const;
	void merge_transition(fstream& fo);
};

} /* namespace GroupDIA */

#endif /* MERGETRANSITIONS_H_ */
