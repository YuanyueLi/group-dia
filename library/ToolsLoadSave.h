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

#ifndef TOOLSLOADSAVE_H_
#define TOOLSLOADSAVE_H_
#include "ToolsSpectrum.h"
namespace GroupDIA {
namespace LoadSave {
namespace IOType {
enum {
	Read, Write, Null
};
}
namespace Type {
enum {
	Basic, Feature_Product, Precursor, Product, Null, All, Basic_Precursor, Basic_Precursor_Iden
};
}
template<class T> void inline io_write(fstream& f, const T& value) {
	f.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

template<class T> void inline io(fstream& f, T& value, int type) {
	if (type == IOType::Write) {
		f.write(reinterpret_cast<const char*>(&value), sizeof(value));
	} else if (type == IOType::Read) {
		f.read(reinterpret_cast<char*>(&value), sizeof(value));
	}
}

template<class T> void inline ion_io(fstream& f, T& ion, int type) {
	ion.load_store(f, type);
}

template<class T> void inline vector_io(fstream& f, vector<T>& vec, int io_type) {
	if (io_type == IOType::Read) {
		int size;
		io(f, size, io_type);
		vec.resize(size);
		if (size != 0)
			f.read(reinterpret_cast<char*>(&vec[0]), size * sizeof(T));
	} else if (io_type == IOType::Write) {
		int size = vec.size();
		io(f, size, io_type);
		if (!vec.empty())
			f.write(reinterpret_cast<const char*>(&vec[0]), vec.size() * sizeof(T));
	}
}

template<class T1, class T2> void inline map_io(fstream& f, map<T1, T2>& vec, int io_type) {
	if (io_type == IOType::Read) {
		int size;
		io(f, size, io_type);
		for (int i = 0; i < size; i++) {
			T1 a;
			T2 b;
			io(f, a, io_type);
			io(f, b, io_type);
			vec.insert(pair<T1, T2>(a, b));
		}
	} else if (io_type == IOType::Write) {
		int size = vec.size();
		io(f, size, io_type);
		for (auto it = vec.begin(); it != vec.end(); ++it) {
			T1 a = it->first;
			T2 b = it->second;
			io(f, a, io_type);
			io(f, b, io_type);
		}
	}
}

template<> void inline io(fstream& f, string& str, int type) {
	if (type == IOType::Write) {
		vector<char> data(str.begin(), str.end());
		vector_io(f, data, type);
	} else if (type == IOType::Read) {
		vector<char> data;
		vector_io(f, data, type);
		str = string(data.begin(), data.end());
	}
}

template<class T> void inline vector_ion_io(fstream& f, vector<T>& vec, int io_type) {
	int size = vec.size();
	io(f, size, io_type);
	if (io_type == IOType::Read) {
		vec = vector<T>(size);
	}
	for (int i = 0; i < size; i++) {
		ion_io(f, vec[i], io_type);
	}
}
}
}

#endif /* TOOLSLOADSAVE_H_ */
