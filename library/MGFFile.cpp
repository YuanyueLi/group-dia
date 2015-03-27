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

#include "MGFFile.h"

namespace GroupDIA {

} /* namespace GroupDIA */

void GroupDIA::MGFFile::store_simple(const String& filename, const String& project_name,
		const PeakMap& experiment) {
	fstream os(filename.c_str(), fstream::out);
	//Output spectrum
	for (auto ispec = experiment.begin(); ispec != experiment.end(); ++ispec) {
		auto& spec = *ispec;
		Precursor precursor;
		if (spec.getPrecursors().size() > 0) {
			precursor = spec.getPrecursors()[0];
		}
		if (spec.getPrecursors().size() > 1) {
			cerr
					<< "Warning: The spectrum written to Mascot file has more than one precursor. The first precursor is used!\n";
		}
		DoubleReal mz(precursor.getMZ()), rt(spec.getRT());

		if (mz == 0) {
			continue;
		} else {
			os << "\n";
			os << "BEGIN IONS\n";
			os << "TITLE=" << project_name << "." << spec.getNativeID() << "." << spec.getNativeID()
					<< ".1\n";
			os << "SCANS=" << spec.getNativeID() << "\n";
			os << "PEPMASS=" << precisionWrapper(mz) << "\n";
			os << "RTINSECONDS=" << precisionWrapper(rt) << "\n";

			bool skip_spectrum_charges(param_.getValue("skip_spectrum_charges").toBool());

			int charge(precursor.getCharge());

			if (charge != 0) {
				if (!skip_spectrum_charges) {
					os << "CHARGE=" << charge << "\n";
				}
			}

			os << "\n";

			for (PeakSpectrum::const_iterator it = spec.begin(); it != spec.end(); ++it) {
				os << precisionWrapper(it->getMZ()) << " " << precisionWrapper(it->getIntensity())
						<< "\n";
			}
			os << "END IONS\n";

		}
	}
}
