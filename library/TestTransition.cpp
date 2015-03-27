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

#include "TestTransition.h"

namespace GroupDIA {

} /* namespace GroupDIA */

void GroupDIA::TestTransition::set_transition(TransitionGroup& t) {
	tr = t;
}

bool GroupDIA::TestTransition::read(String& id, fstream& fi) {
	while (!fi.eof()) {
		TransitionGroup t;
		t.load_store(fi, LoadSave::IOType::Read, LoadSave::Type::All);
		if (String(t.scan_num) == id) {
			tr = t;
			return true;
		}
	}
	return false;
}

void GroupDIA::TestTransition::write_info(fstream& fo) {
	for (auto it = tr.precursor_ions.begin(); it != tr.precursor_ions.end(); ++it) {
		vector<double> rt;
		it->get_nearby_rt(rt);
		fo << rt.size() << "\t";
	}
	fo << endl;

	vector<double> mz;
	tr.get_product_mz(mz);
	for (auto it = mz.begin(); it != mz.end(); ++it) {
		fo << *it << "\t";
	}
	fo << endl;

	for (auto it = tr.precursor_ions.begin(); it != tr.precursor_ions.end(); ++it) {
		int start, length;
		it->get_peak_rt_in_nearby_no(start, length);
		fo << start << "\t" << length << "\t";
	}
	fo << endl;

	fo << tr.feature_exp_num << endl;
}

void GroupDIA::TestTransition::write(fstream& fo) {
	//Output ms1

	//output rt
	for (auto it = tr.precursor_ions.begin(); it != tr.precursor_ions.end(); ++it) {
		vector<double> rt;
		it->get_nearby_rt(rt);
		for (auto jt = rt.begin(); jt != rt.end(); ++jt) {
			fo << *jt << "\t";
		}
	}
	fo << endl;

	//output intensity
	for (auto it = tr.precursor_ions.begin(); it != tr.precursor_ions.end(); ++it) {
		vector<double> intensity;
		it->get_nearby_intensity(intensity);
		for (auto jt = intensity.begin(); jt != intensity.end(); ++jt) {
			fo << *jt << "\t";
		}
	}
	fo << endl;

	//Output ms2
	//output rt
	for (auto it = tr.product_ions.begin(); it != tr.product_ions.end(); ++it) {
		vector<double> rt;
		it->get_rt(rt);
		for (auto jt = rt.begin(); jt != rt.end(); ++jt) {
			fo << *jt << "\t";
		}
	}
	fo << endl;

	//output intensity
	for (int mz = 0; mz < tr.get_product_num(); mz++) {
		for (auto it = tr.product_ions.begin(); it != tr.product_ions.end(); ++it) {
			vector<double> intensity;
			it->get_intensity(mz, intensity);
			for (auto jt = intensity.begin(); jt != intensity.end(); ++jt) {
				fo << *jt << "\t";
			}
		}
		fo << endl;
	}
	fo << endl;

}

void GroupDIA::OutputToOpenswath::read_transition(TransitionIdentified& t, vector<ReactionMonitoringTransition>& trans) {
	vector<double> product_mz;
	vector<double> product_intensity;
	vector<double> product_charge;
	vector<string> product_type;

	t.get_product_mz(product_mz);
	double precursor_intensity;
	t.get_library_intensity(precursor_intensity, product_intensity);
	vector<string> type;
	t.get_product_type(type);
	for (int i = 0; i < type.size(); i++) {
		String s = type[i];
		vector<String> ss;
		s.split("_", ss);
		product_charge.push_back(ss[0].toDouble());
		product_type.push_back(ss[1]);
	}

	trans.clear();
	for (int i = 0; i < product_mz.size(); i++) {
		ReactionMonitoringTransition rm_trans;
		rm_trans.setNativeID(String(t.get_uid()) + String(i));
		rm_trans.setPrecursorMZ(t.get_precursor_mz());
		rm_trans.setProductMZ(product_mz[i]);
		rm_trans.setPeptideRef(String(t.get_uid()));
		rm_trans.setLibraryIntensity(product_intensity[i]);
		{
			OpenMS::ReactionMonitoringTransition::Product p = rm_trans.getProduct();
			p.setChargeState(product_charge[i]);
			rm_trans.setProduct(p);
		}
		// add interpretation
		{
			OpenMS::ReactionMonitoringTransition::Product p = rm_trans.getProduct();
			CVTermList interpretation;
			CVTerm rank;
			rank.setCVIdentifierRef("MS");
			rank.setAccession("MS:1000926");
			rank.setName("product interpretation rank");
			rank.setValue(1); // we only store the best interpretation
			CVTerm frag_nr;
			frag_nr.setCVIdentifierRef("MS");
			frag_nr.setAccession("MS:1000903");
			frag_nr.setName("product ion series ordinal");
			frag_nr.setValue(String(i));

			// figure out which fragment it is
			if (product_type[i] == "y") {
				CVTerm ion;
				ion.setCVIdentifierRef("MS");
				ion.setAccession("MS:1001220");
				ion.setName("frag: y ion");
				interpretation.addCVTerm(ion);
			} else if (product_type[i] == "b") {
				CVTerm ion;
				ion.setCVIdentifierRef("MS");
				ion.setAccession("MS:1001224");
				ion.setName("frag: b ion");
				interpretation.addCVTerm(ion);
			}

			interpretation.addCVTerm(rank);
			interpretation.addCVTerm(frag_nr);
			p.addInterpretation(interpretation);
			rm_trans.setProduct(p);
		}

		if (t.get_decoy_type() == 0) {
			rm_trans.setDecoyTransitionType(ReactionMonitoringTransition::TARGET);
		} else {
			rm_trans.setDecoyTransitionType(ReactionMonitoringTransition::DECOY);
		}

		rm_trans.setMetaValue("annotation", "");
		trans.push_back(rm_trans);
	}
}

void GroupDIA::OutputToOpenswath::read_transition(TransitionIdentified& t, int exp_num, ToolsTabFileWrite& out) {
	if (t.get_decoy_type() == TransitionIdentified::Decoy_Rand)
		return;

	vector<double> product_mz;
	vector<double> product_intensity;
	vector<double> product_charge;
	vector<string> product_type;

	t.get_product_mz(product_mz);
	double precursor_intensity;
	t.get_library_intensity(precursor_intensity, product_intensity);
	vector<string> type;
	t.get_product_type(type);
	for (int i = 0; i < type.size(); i++) {
		String s = type[i];
		vector<String> ss;
		s.split("_", ss);
		product_charge.push_back(ss[0].toDouble());
		product_type.push_back(ss[1] + ss[2]);
	}

	for (int i = 0; i < product_mz.size(); i++) {
		out.add_item("PrecursorMz", t.get_precursor_mz());
		out.add_item("ProductMz", product_mz[i]);
		out.add_item("Tr_recalibrated", t.get_precursor_rt(exp_num));
		out.add_item("transition_name",
				t.get_accession() + "_" + t.get_pep_seq() + "/" + String(t.get_uid()) + "_" + product_type[i] + "_"
						+ String(product_charge[i]));
		out.add_item("CE", 0);
		out.add_item("LibraryIntensity", product_intensity[i]);
		out.add_item("transition_group_id", t.get_accession() + "_" + t.get_pep_seq() + "/" + String(t.get_uid()));
		out.add_item("decoy", 0);
		AASequence aas;
		t.get_pep_aaseq(aas);
		out.add_item("PeptideSequence", aas.toUnmodifiedString());
		out.add_item("ProteinName", t.get_accession());
		out.add_item("Annotation", product_type[i]);
		out.add_item("FullUniModPeptideName", t.get_pep_seq());
		out.add_item("MissedCleavages", 0);
		out.add_item("Replicates", 0);
		out.add_item("NrModifications", 0);
		out.add_item("PrecursorCharge", t.get_precursor_charge());
		out.add_item("GroupLabel", "light");

		out.end_add_item();
	}
}
