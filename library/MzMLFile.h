// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef MZML_H
#define MZML_H

#include <OpenMS/FORMAT/XMLFile.h>
#include "MzMLHandler.h"
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <iostream>
namespace GroupDIA {
/**
 @brief File adapter for MzML files

 This implementation does currently not support the whole functionality of MzML.
 Some minor features are still missing:
 - chromatograms

 @ Implement chromatograms (Andreas)

 @ingroup FileIO
 */
class MzMLFile: public OpenMS::Internal::XMLFile, public OpenMS::ProgressLogger {
public:
	///Default constructor
	MzMLFile();
	///Destructor
	~MzMLFile();

	/// Mutable access to the options for loading/storing
	OpenMS::PeakFileOptions& getOptions();

	/// Non-mutable access to the options for loading/storing
	const OpenMS::PeakFileOptions& getOptions() const;

	/// set options for loading/storing
	void setOptions(const OpenMS::PeakFileOptions &);

	/**
	 @brief Loads a map from a MzML file.

	 @p map has to be a MSExperiment or have the same interface.

	 @exception Exception::FileNotFound is thrown if the file could not be opened
	 @exception Exception::ParseError is thrown if an error occurs during parsing
	 */
	template<typename MapType>
	void load(const OpenMS::String& filename, MapType& map) {
		map.reset();

		//set DocumentIdentifier
		map.setLoadedFileType(filename);
		map.setLoadedFilePath(filename);

		MzMLHandler<MapType> handler(map, filename, getVersion(), *this);
		handler.setOptions(options_);
		//handler can throw parse error and other errors - catch those here - they are the cause for a parse error - report accordingly
		try {

			parse_(filename, &handler);
		} catch (OpenMS::Exception::BaseException& e) {
			std::string expr;
			expr.append(e.getFile());
			expr.append("@");
			std::stringstream ss;
			ss << e.getLine(); // we need c++11!! maybe in 2012?
			expr.append(ss.str());
			expr.append("-");
			expr.append(e.getFunction());
			std::string mess = "- due to that error of type ";
			mess.append(e.getName());
			throw OpenMS::Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, expr,
					mess);
		}

	}

	/**
	 @brief Stores a map in a MzML file.

	 @p map has to be a MSExperiment or have the same interface.

	 @exception Exception::UnableToCreateFile is thrown if the file could not be created
	 */
	template<typename MapType>
	void store(const OpenMS::String& filename, const MapType& map) const {
		MzMLHandler<MapType> handler(map, filename, getVersion(), *this);
		handler.setOptions(options_);
		save_(filename, &handler);
	}

	/**
	 @brief Checks if a file validates against the XML schema.

	 @exception Exception::FileNotFound is thrown if the file cannot be found.
	 */
	bool isValid(const OpenMS::String& filename, std::ostream& os = std::cerr);

	/**
	 @brief Checks if a file is valid with respect to the mapping file and the controlled vocabulary.

	 @param filename File name of the file to be checked.
	 @param errors Errors during the validation are returned in this output parameter.
	 @param warnings Warnings during the validation are returned in this output parameter.

	 @exception Exception::FileNotFound is thrown if the file could not be opened
	 */
	bool isSemanticallyValid(const OpenMS::String& filename, OpenMS::StringList& errors,
			OpenMS::StringList& warnings);

private:

	/// Options for loading / storing
	OpenMS::PeakFileOptions options_;

	/// Location of indexed mzML schema
	OpenMS::String indexed_schema_location_;
};

}
// namespace LISwath

#endif // MZMLFILE_H
