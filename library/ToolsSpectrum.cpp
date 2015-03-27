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
#include "ToolsSpectrum.h"

namespace GroupDIA {

bool ToolsSpectrum::is_target_identification_result(const vector<String>& acc) {
	for (auto& it : acc) {
		if (it.find("DECOY") == it.npos) {
			return true;
		}
	}
	return false;
}

void ToolsSpectrum::peak_finder(const vector<double>& intensity, int peak_finder_start,
		int peak_finder_length, vector<PeakInfo>& peak_info) {
	peak_info.clear();
	auto intensity_start = intensity.begin();
	auto intensity_end = intensity.end();

	auto range_start = intensity_start + peak_finder_start;
	auto range_end = range_start + peak_finder_length;
	auto ipre = range_start;
	if (ipre == range_end) {
		return;
	}
	auto ic = ipre + 1;
	if (ic == range_end) {
		return;
	}
	auto inext = ic + 1;
	if (inext == range_end)
		return;
	int peak_start = 0;
	int cur_peak_no = peak_finder_start + 1;

	//Move to left
	while (ipre != intensity.begin()) {
		if (*ic >= *inext and *ic > *ipre) {
			break;
		}
		--ipre, --ic, --inext, --cur_peak_no;
	}
	range_start = ipre;

	//Start find
	ipre = range_start;
	ic = ipre + 1;
	inext = ic + 1;

	while (inext != intensity.end()) {
		if (*ic >= *inext and *ic > *ipre) {
			//Find peak start
			int temp_peak_no = cur_peak_no;
			while (ipre != intensity_start) {
				--ic, --ipre, --cur_peak_no;
				if (*ic == 0 or *ipre >= *ic) {
					//Peak start
					peak_start = cur_peak_no;
					break;
				}
			}
			if (peak_start > peak_finder_start + peak_finder_length)
				break;
			cur_peak_no = temp_peak_no + 1;
			ic = inext;
			ipre = ic - 1;
			++inext;
			//Find peak end
			if (inext == intensity_end) {
				int peak_end = cur_peak_no;
				PeakInfo peak(peak_start, peak_end - peak_start);
				peak_info.push_back(peak);
				break;
			}
			while (inext != intensity_end) {
				if (*ic == 0 or *inext >= *ic) {
					int peak_end = cur_peak_no;
					PeakInfo peak(peak_start, peak_end - peak_start);
					peak_info.push_back(peak);
					break;
				}
				++ic, ++inext, ++cur_peak_no;
			}
			cur_peak_no = temp_peak_no + 1;
			ic = ipre + 1;
			inext = ic + 1;
			if (inext == range_end)
				break;
		}
		++ic, ++inext, ++ipre, ++cur_peak_no;
	}
}

ToolsSpectrum::ToolsSpectrum() {
}

ToolsSpectrum::~ToolsSpectrum() {
}

void ToolsSpectrum::mean_and_std_dev(const vector<double>& values, double& mean, double& sd) {
	if (values.empty()) {
		mean = 0;
		sd = 1;
	} else {
		mean = Math::mean(values.begin(), values.end());
		double sq_sum = 0;
		for (auto it = values.begin(); it != values.end(); ++it) {
			sq_sum += (*it - mean) * (*it - mean);
		}
		sd = sqrt(sq_sum / (double) (values.size() - 1));
	}
}

void ToolsSpectrum::median_and_sad(vector<double> v, double& median, double& sad) {
	if (v.empty()) {
		median = 0;
		sad = 0;
	} else {
		nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
		median = v[v.size() / 2];
		vector<double> d = v;
		for (auto &i : d)
			i = abs(i - median);
		nth_element(d.begin(), d.begin() + d.size() / 2, d.end());
		sad = d[d.size() / 2];
	}
}

void ToolsSpectrum::normalize_sum(vector<double>& x) {
	double sumx = 0;
	for (auto i = x.begin(); i != x.end(); ++i) {
		sumx += *i;
	}
	if (sumx == 0.0) {
		return;
	} // do not divide by zero
	for (unsigned int i = 0; i < x.size(); i++) {
		x.at(i) = x.at(i) / sumx;
	}
}

void ToolsSpectrum::normalize_max(vector<double>& x) {
	if (x.empty())
		return;
	double max = *max_element(x.begin(), x.end());
	if (max == 0.0) {
		return;
	} // do not divide by zero
	for (unsigned int i = 0; i < x.size(); i++) {
		x.at(i) = x.at(i) / max;
	}
}

int ToolsSpectrum::get_rt_no(const vector<double>& rt, double r) {
	int num = 0;
	for (auto it_c = rt.begin(); it_c != rt.end(); ++it_c, ++num) {
		if (*it_c > r) {
			if (it_c != rt.begin()) {
				double delta_next = *it_c - r;
				--it_c;
				double delta_pre = r - *it_c;
				if (delta_next > delta_pre)
					return num - 1;
				else
					return num;
			} else
				return num;
		}
	}
	return num - 1;
}

double ToolsSpectrum::get_rt(const vector<double>& rt, int r) {
	if (r < 0)
		return rt.at(0);
	if (r >= rt.size()) {
		return rt.at(rt.size() - 1);
	} else {
		return rt.at(r);
	}
}

void ToolsSpectrum::standardize_data(std::vector<double>& data) {
// subtract the mean and divide by the standard deviation
	double mean = std::accumulate(data.begin(), data.end(), 0.0) / (double) ((data.size()));
	double sqsum = 0;
	for (std::vector<double>::iterator it = data.begin(); it != data.end(); it++) {
		sqsum += (*it - mean) * (*it - mean);
	}
	double std = sqrt(sqsum / data.size()); // standard deviation

	for (std::size_t i = 0; i < data.size(); i++) {
		data.at(i) = (data.at(i) - mean) / std;
	}
}

void ToolsSpectrum::get_best_delay(vector<double> a, vector<double> b, int max_delay,
		vector<pair<int, double> >& dela_score) {
	int a_size = a.size();
	int b_size = b.size();
	if (a_size != b_size or b_size == 0 or a_size == 0)
		throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"Error in calc the cross correlation, vector size error.");

	standardize_data(a);
	standardize_data(b);

	dela_score.clear();
	for (int delay = -max_delay; delay <= max_delay; delay++) {
		double sxy = 0;
		for (int i = 0; i < a_size; i++) {
			int j = i + delay;
			if (j < 0 || j >= b_size) {
				continue;
			}
			sxy += a[i] * b[j];
		}

		dela_score.push_back(pair<int, double>(delay, sxy));
	}
}

bool ToolsSpectrum::is_same_in_ppm(const double a, const double b, const double ppm) {
	return abs(a - b) * 1000 * 1000 / a < ppm;
}

void ToolsSpectrum::get_pearson_score(const vector<double>& precursor_intensity,
		const vector<vector<double> >& product_intensity, vector<double>& score) {
	score.resize(product_intensity.size());
	auto ipro = product_intensity.begin();
	auto is = score.begin();
	for (; is != score.end(); ++is, ++ipro) {
		if (precursor_intensity.size() != ipro->size()) {
			cout << "Error in pearson_score" << endl;
			cout << precursor_intensity.size() << "\t" << ipro->size() << endl;
			throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
					"Error in pearson_score");
		}
		*is = Math::pearsonCorrelationCoefficient(precursor_intensity.begin(),
				precursor_intensity.end(), ipro->begin(), ipro->end());
	}
}

void ToolsSpectrum::convert_aasequence_to_peptide(const AASequence& aa_sequence,
		TargetedExperiment::Peptide& result_peptide) {
	ModificationsDB* mod_db = ModificationsDB::getInstance();
	std::vector<TargetedExperiment::Peptide::Modification> mods;
	result_peptide.sequence = aa_sequence.toUnmodifiedString();

	// in TraML, the modification the AA starts with residue 1 but the
	// OpenMS objects start with zero -> we start counting with zero here
	// and the TraML handler will add 1 when storing the file.
	if (aa_sequence.isValid() && aa_sequence.toString() != aa_sequence.toUnmodifiedString()) {
		if (!aa_sequence.getNTerminalModification().empty()) {
			ResidueModification rmod = mod_db->getTerminalModification(
					aa_sequence.getNTerminalModification(), ResidueModification::N_TERM);
			add_modification_(mods, -1, rmod, aa_sequence.getNTerminalModification());
		}
		if (!aa_sequence.getCTerminalModification().empty()) {
			ResidueModification rmod = mod_db->getTerminalModification(
					aa_sequence.getCTerminalModification(), ResidueModification::C_TERM);
			add_modification_(mods, aa_sequence.size(), rmod,
					aa_sequence.getCTerminalModification());
		}
		for (Size i = 0; i != aa_sequence.size(); i++) {
			if (aa_sequence[i].isModified()) {
				// search the residue in the modification database (if the sequence is valid, we should find it)
				TargetedExperiment::Peptide::Modification mod;
				ResidueModification rmod = mod_db->getModification(
						aa_sequence.getResidue(i).getOneLetterCode(),
						aa_sequence.getResidue(i).getModification(), ResidueModification::ANYWHERE);
				add_modification_(mods, i, rmod, aa_sequence.getResidue(i).getModification());
			}
		}
	}

	result_peptide.mods = mods;
}

void ToolsSpectrum::add_modification_(std::vector<TargetedExperiment::Peptide::Modification>& mods,
		int location, ResidueModification& rmod, const String& name) {
	TargetedExperiment::Peptide::Modification mod;
	String unimod_str = rmod.getUniModAccession();
	mod.location = location;
	mod.mono_mass_delta = rmod.getDiffMonoMass();
	mod.avg_mass_delta = rmod.getDiffAverageMass();
	// CV term with the full unimod accession number and name
	CVTerm unimod_name;
	unimod_name.setCVIdentifierRef("UNIMOD");
	unimod_name.setAccession(unimod_str.toUpper());
	unimod_name.setName(name);
	mod.addCVTerm(unimod_name);
	mods.push_back(mod);
}
}/* namespace GroupDIA */
