// $Id$
//
//  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <Geometry/point.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/BadFileException.h>
#include <boost/foreach.hpp>
#include "ConformerParser.h"

#include <fstream>


namespace RDKit{
  namespace ConformerParser {

    INT_VECT addConformersFromList(ROMol &mol, const std::vector<std::vector<double> > &coords,
                               int numConf){
      unsigned int numAtomsPerConf = mol.getNumAtoms();
      unsigned int numCoordPerConf = 3 * numAtomsPerConf;
      PRECONDITION(numConf <= int(coords.size()), "numConf greater than number of conformations");
      if (numConf < 0) {
        numConf = coords.size();
      }
      // loop over the conformers
      INT_VECT confIds;
      for (unsigned int i = 0; i < numConf; ++i) {
        if (coords[i].size() != numCoordPerConf) {
          throw ValueErrorException("Wrong number of coordinates");
        }
        RDKit::Conformer *conf = new RDKit::Conformer(numAtomsPerConf);
        // loop over atoms
        for (unsigned int atom = 0; atom < numAtomsPerConf; ++atom) {
          //RDGeom::Point3D p(coords[i][3*atom], coords[i][3*atom+1], coords[i][3*atom+2]);
          (*conf).setAtomPos(atom, RDGeom::Point3D(coords[i][3*atom], coords[i][3*atom+1], coords[i][3*atom+2]));
        }
        int confId = mol.addConformer(conf, true);
        confIds.push_back(confId);
      }
      return confIds;
    }

    void readAmberTrajectory(std::string fName, std::vector<std::vector<double> > &coords,
                             unsigned int numAtoms) {
      std::ifstream inStream(fName.c_str());
      if (!inStream || (inStream.bad()) ) {
        std::ostringstream errout;
        errout << "Bad input file " << fName;
        throw BadFileException(errout.str());
      }

      // clear coords
      coords.resize(0);

      std::string tempStr;
      // title
      std::getline(inStream, tempStr);
      // read coordinates
      std::vector<double> tmpCoords;
      while (true) {
        double c;
        if (!(inStream >> c)) {
          if (!inStream.eof()) {
            throw ValueErrorException("Error while reading file");
          }
          break;
        }
        tmpCoords.push_back(c);
      }
      // convert to conformers
      unsigned int numCoordsPerConf = 3 * numAtoms;
      if (tmpCoords.size() % numCoordsPerConf != 0) {
        throw ValueErrorException("Wrong number of coordinates");
      }
      unsigned int numConfs = tmpCoords.size() / numCoordsPerConf;
      unsigned int c = 0;
      for (unsigned int i = 0; i < numConfs; ++i) {
        std::vector<double> coordConf(numCoordsPerConf);
        for (unsigned int atom = 0; atom < numCoordsPerConf; ++atom, ++c) {
          coordConf[atom] = tmpCoords[c];
        }
        coords.push_back(coordConf);
      }
    }

    void readGromosTrajectory(std::string fName, std::vector<std::vector<double> > &coords,
							 unsigned int numAtoms) {
	  std::ifstream inStream(fName.c_str());
	  if (!inStream || (inStream.bad()) ) {
		std::ostringstream errout;
		errout << "Bad input file " << fName;
		throw BadFileException(errout.str());
	  }

	  // clear coords
	  coords.resize(0);

	  unsigned int numCoordsPerConf = 3 * numAtoms; // number of coordinates per conformer

	  std::string tempStr;
	  while (!inStream.eof()) {
		  std::getline(inStream, tempStr);
		  //std::cerr << tempStr << std::endl;
		  if (tempStr == "TITLE") { // title block - will be ignored
			  while(tempStr != "END") {
				  std::getline(inStream, tempStr);
			  }
		  } else if (tempStr == "TIMESTEP") { // timestep block - will be ignored
			  while(tempStr != "END") {
				  std::getline(inStream, tempStr);
			  }
		  } else if (tempStr == "POSITIONRED") { // these are the positions
			  std::vector<double> coordConf;
			  for (unsigned int i = 0; i < numAtoms; ++i) {
				  std::getline(inStream, tempStr);
				  if (inStream.eof() || tempStr == "END") {
					  throw ValueErrorException("Wrong number of coordinates");
				  }
				  if (tempStr.find("#") != std::string::npos) { // ignore comments
					  --i; // reset the atom counter
					  continue;
				  }
				  std::stringstream ls(tempStr);
				  double x, y, z;
				  if (!(ls >> x >> y >> z)) {
					  throw ValueErrorException("Error while reading file");
				  }
				  // store the coordinates
				  coordConf.push_back(x);
				  coordConf.push_back(y);
				  coordConf.push_back(z);
			  }
			  if (coordConf.size() != numCoordsPerConf) {
				  throw ValueErrorException("Wrong number of coordinates");
			  }
			  std::getline(inStream, tempStr); // the END line
			  if (tempStr != "END") {
				  throw ValueErrorException("Wrong number of coordinates");
			  }
			  coords.push_back(coordConf);
		  } else if (tempStr == "GENBOX") { // box information block - will be ignored
			  while(tempStr != "END") {
				  std::getline(inStream, tempStr);
			  }
		  } else {
			  if (!inStream.eof()) {
			      throw ValueErrorException("Unsupported block in file: "+tempStr+". Supported blocks are TITLE, TIMESTEP, POSITIONRED, GENBOX.");
			  }
		  }
	  } // read file
    }

  } // end namespace ConformerParser
} // end namespace RDKit
