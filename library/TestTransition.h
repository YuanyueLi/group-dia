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

#ifndef TESTTRANSITION_H_
#define TESTTRANSITION_H_

#include "TransitionIdentified.h"
#include "ToolsFilename.h"
#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
#include "ToolsTabFileWrite.h"

namespace GroupDIA {

class TestTransition {
public:
	void set_transition(TransitionGroup& t);

	bool read(String& id, fstream& fi);
	void write(fstream& fo);
	void write_info(fstream& fo);
private:
	TransitionGroup tr;
};

class OutputToOpenswath {
public:
	void read_transition(TransitionIdentified& t, vector<ReactionMonitoringTransition>& trans);
	void read_transition(TransitionIdentified& t, int exp_num, ToolsTabFileWrite& out);
};

} /* namespace GroupDIA */

#endif /* TESTTRANSITION_H_ */
