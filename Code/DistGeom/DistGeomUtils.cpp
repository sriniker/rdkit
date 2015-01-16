// $Id$
//
//  Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "BoundsMatrix.h"
#include "DistGeomUtils.h"
#include "DistViolationContrib.h"
#include "ChiralViolationContrib.h"
#include "FourthDimContrib.h"
#include <Numerics/Matrix.h>
#include <Numerics/SymmMatrix.h>
#include <Numerics/Vector.h>
#include <RDGeneral/Invariant.h>
#include <Numerics/EigenSolvers/PowerEigenSolver.h>
#include <RDGeneral/utils.h>
#include <ForceField/ForceField.h>
#include <ForceField/UFF/DistanceConstraint.h>
#include <ForceField/UFF/TorsionConstraint.h>

namespace DistGeom {
  const double EIGVAL_TOL=0.001;


  double pickRandomDistMat(const BoundsMatrix &mmat, 
                           RDNumeric::SymmMatrix<double> &distMat,
                           int seed) {
    if(seed>0){
      RDKit::getRandomGenerator(seed);
    }
    return pickRandomDistMat(mmat,distMat,RDKit::getDoubleRandomSource());
  }

  double pickRandomDistMat(const BoundsMatrix &mmat, 
                           RDNumeric::SymmMatrix<double> &distMat,
                           RDKit::double_source_type &rng) {
    // make sure the sizes match up
    unsigned int npt = mmat.numRows();
    CHECK_INVARIANT(npt == distMat.numRows(), "Size mismatch");

    double largestVal=-1.0;
    double *ddata = distMat.getData();
    for (unsigned int i = 1; i < npt; i++) {
      unsigned int id = i*(i+1)/2;
      for (unsigned int j = 0; j < i; j++) {
        double ub = mmat.getUpperBound(i,j);
        double lb = mmat.getLowerBound(i,j);
        CHECK_INVARIANT(ub >= lb, "");
        double rval = rng();
        //std::cerr<<i<<"-"<<j<<": "<<rval<<std::endl;
        double d = lb + (rval)*(ub - lb);
        ddata[id+j] = d;
        if(d>largestVal){
          largestVal=d;
        }
      }
    }
    return largestVal;
  }

  bool computeInitialCoords(const RDNumeric::SymmMatrix<double> &distMat,  
                            RDGeom::PointPtrVect &positions, bool randNegEig, 
                            unsigned int numZeroFail,
                            int seed) {
    if(seed>0){
      RDKit::getRandomGenerator(seed);
    }
    return computeInitialCoords(distMat,positions,RDKit::getDoubleRandomSource(),
                                randNegEig,numZeroFail);

  }
  bool computeInitialCoords(const RDNumeric::SymmMatrix<double> &distMat,  
                            RDGeom::PointPtrVect &positions,
                            RDKit::double_source_type &rng,
                            bool randNegEig, 
                            unsigned int numZeroFail){
    unsigned int N = distMat.numRows();
    unsigned int nPt = positions.size();
    CHECK_INVARIANT(nPt == N, "Size mismatch");
    
    unsigned int dim = positions.front()->dimension();
    
    const double *data = distMat.getData();
    RDNumeric::SymmMatrix<double> sqMat(N), T(N, 0.0); 
    RDNumeric::DoubleMatrix eigVecs(dim,N);
    RDNumeric::DoubleVector eigVals(dim);
    
    double *sqDat = sqMat.getData();
    
    unsigned int dSize = distMat.getDataSize();
    double sumSqD2  = 0.0;      
    for (unsigned int i = 0; i < dSize; i++) {
      sqDat[i] = data[i]*data[i];
      sumSqD2 += sqDat[i];
    }
    sumSqD2 /= (N*N);

    RDNumeric::DoubleVector sqD0i(N, 0.0);
    double *sqD0iData = sqD0i.getData();
    for (unsigned int i = 0; i < N; i++) {
      for (unsigned int j = 0; j < N; j++) {
        sqD0iData[i] += sqMat.getVal(i,j);
      }
      sqD0iData[i] /= N;
      sqD0iData[i] -= sumSqD2;
      
      if ((sqD0iData[i] < EIGVAL_TOL) && (N > 3)){
        return false;
      }
    }

    for (unsigned int i = 0; i < N; i++) {
      for (unsigned int j = 0; j <= i; j++) {
        double val = 0.5*(sqD0iData[i] + sqD0iData[j] - sqMat.getVal(i,j));
        T.setVal(i,j, val);
      }
    }
    unsigned int nEigs = (dim < N) ? dim : N;
    RDNumeric::EigenSolvers::powerEigenSolver(nEigs, T, eigVals, eigVecs,
                                              (int)(sumSqD2*N));
    
    double *eigData = eigVals.getData();
    bool foundNeg = false;
    unsigned int zeroEigs = 0;
    for (unsigned int i = 0; i < dim; i++) {
      if (eigData[i] > EIGVAL_TOL) {
        eigData[i] = sqrt(eigData[i]);
      } else if (fabs(eigData[i]) < EIGVAL_TOL) {
        eigData[i] = 0.0;
        zeroEigs++;
      } else {
        foundNeg = true;
      }
    }
    if ((foundNeg) && (!randNegEig) ) {
      return false;
    }

    if ((zeroEigs >= numZeroFail) && (N > 3)) {
      return false;
    }

    for (unsigned int i = 0; i < N; i++) {
      RDGeom::Point *pt = positions[i];
      for (unsigned int j = 0; j < dim; ++j) {
        if (eigData[j] >= 0.0) {
          (*pt)[j] = eigData[j]*eigVecs.getVal(j,i);
        } else {
          //std::cerr<<"!!! "<<i<<"-"<<j<<": "<<eigData[j]<<std::endl;
          (*pt)[j] = 1.0 - 2.0*rng();
        }
      }
    }
    return true;
  }

  bool computeRandomCoords(RDGeom::PointPtrVect &positions, double boxSize,
                           int seed){
    if(seed>0){
      RDKit::getRandomGenerator(seed);
    }
    return computeRandomCoords(positions,boxSize,
                               RDKit::getDoubleRandomSource());

  }
  bool computeRandomCoords(RDGeom::PointPtrVect &positions, double boxSize,
                           RDKit::double_source_type &rng){
    CHECK_INVARIANT(boxSize>0.0, "bad boxSize");
    
    for(RDGeom::PointPtrVect::iterator ptIt=positions.begin();
        ptIt!=positions.end();++ptIt){
      RDGeom::Point *pt = *ptIt;
      for (unsigned int i = 0; i<pt->dimension(); ++i) {
        (*pt)[i]=boxSize*(rng()-0.5);
      }
    }
    return true;
  }

  
  ForceFields::ForceField *constructForceField(const BoundsMatrix &mmat,
                                               RDGeom::PointPtrVect &positions,
                                               const VECT_CHIRALSET & csets,
                                               double weightChiral,
                                               double weightFourthDim,
                                               std::map< std::pair<int,int>,double> *extraWeights,
                                               double basinSizeTol) {
    unsigned int N = mmat.numRows();
    CHECK_INVARIANT(N == positions.size(), "");
    ForceFields::ForceField *field=new ForceFields::ForceField(positions[0]->dimension());
    for(unsigned int i=0; i < N; i++){
      field->positions().push_back(positions[i]);
    }
    
    for (unsigned int i = 1; i < N; i++) {
      for (unsigned int j = 0; j < i; j++) {
        double w = 1.0;
        double l = mmat.getLowerBound(i,j);
        double u = mmat.getUpperBound(i,j);
        bool includeIt=false;
        if(extraWeights){
          std::map< std::pair<int,int>,double>::const_iterator mapIt;
          mapIt = extraWeights->find(std::make_pair(i,j));
          if(mapIt != extraWeights->end()){
            w = mapIt->second;
            includeIt=true;
          }
        }
        if(u-l <= basinSizeTol) {
          includeIt=true;
        }
        if(includeIt){
          DistViolationContrib *contrib = new DistViolationContrib(field, i, j, u, l, w);
          field->contribs().push_back(ForceFields::ContribPtr(contrib));
        }
      }
    }

    // now add chiral constraints
    if (weightChiral > 1.e-8) {
      for (VECT_CHIRALSET::const_iterator csi = csets.begin();
           csi != csets.end(); csi++) {
        ChiralViolationContrib *contrib = new ChiralViolationContrib(field, csi->get(),
                                                                     weightChiral);
        field->contribs().push_back(ForceFields::ContribPtr(contrib));
      }
    }

    // finally the contribution from the fourth dimension if we need to
    if ((field->dimension() == 4) && (weightFourthDim > 1.e-8)) {
      for (unsigned int i = 1; i < N; i++) {
        FourthDimContrib *contrib = new FourthDimContrib(field,i,weightFourthDim);
        field->contribs().push_back(ForceFields::ContribPtr(contrib));
      }
    }
    return field;
  }

  ForceFields::ForceField *construct3DForceField(const BoundsMatrix &mmat,
                                               RDGeom::Point3DPtrVect &positions,
                                               std::vector<std::pair<int, int> > &bonds,
                                               std::vector<std::pair<int, int> > &angles,
                                               std::vector<std::vector<int> > &expTorsionAtoms,
                                               std::vector<std::pair<double, double> > &expTorsionAngles,
                                               double basinSizeTol) {
    unsigned int N = mmat.numRows();
    CHECK_INVARIANT(N == positions.size(), "");
    CHECK_INVARIANT(expTorsionAtoms.size() == expTorsionAngles.size(), "");
    ForceFields::ForceField *field = new ForceFields::ForceField(positions[0]->dimension());
    for (unsigned int i = 0; i < N; ++i){
      field->positions().push_back(positions[i]);
    }

    // current angles
    std::vector<double> currentAngles;
    for (unsigned int t = 0; t < expTorsionAngles.size(); ++t) {
      RDGeom::Point3D p1 = (*positions[expTorsionAtoms[t][0]]);
      RDGeom::Point3D p2 = (*positions[expTorsionAtoms[t][1]]);
      RDGeom::Point3D p3 = (*positions[expTorsionAtoms[t][2]]);
      RDGeom::Point3D p4 = (*positions[expTorsionAtoms[t][3]]);
      double ang = RDGeom::computeSignedDihedralAngle(p1, p2, p3, p4);
      //std::cout << ang/M_PI * 180.0 << std::endl;
      currentAngles.push_back(ang);
    }

    // torsion constraints
    double ftorsion = 1000.0; // force constant
    for (unsigned int t = 0; t < expTorsionAtoms.size(); ++t) {
    	int ang1 = expTorsionAngles[t].first/M_PI * 180.0;
    	int ang2 = expTorsionAngles[t].second/M_PI * 180.0;
    	int i = expTorsionAtoms[t][0];
    	int j = expTorsionAtoms[t][1];
    	int k = expTorsionAtoms[t][2];
    	int l = expTorsionAtoms[t][3];
    	if (currentAngles[t] < 0) {
    	  double tmp = -ang2;
    	  ang2 = -ang1; ang1 = tmp;
    	}
    	//std::cout << ang1 << " " << ang2 << std::endl;
    	ForceFields::UFF::TorsionConstraintContrib *contrib = new ForceFields::UFF::TorsionConstraintContrib(field, i, j, k, l, ang1, ang2, ftorsion);
    	field->contribs().push_back(ForceFields::ContribPtr(contrib));
    }

    // 1,2 distance constraints
	  double fdist = 100.0; // force constant
	  std::vector<std::pair<int, int> >::iterator bi;
	  for (bi = bonds.begin(); bi != bonds.end(); bi++) {
		  unsigned int i = bi->first;
		  unsigned int j = bi->second;
		  double d = ((*positions[i]) - (*positions[j])).length();
		  double l = d-0.01;
		  double u = d+0.01;
		  ForceFields::UFF::DistanceConstraintContrib *contrib = new ForceFields::UFF::DistanceConstraintContrib(field, i, j, l, u, fdist);
		  field->contribs().push_back(ForceFields::ContribPtr(contrib));
	  }

	  // 1,3 distance constraints
		for (bi = angles.begin(); bi != angles.end(); bi++) {
		  unsigned int i = bi->first;
		  unsigned int j = bi->second;
		  double d = ((*positions[i]) - (*positions[j])).length();
		  double l = d-0.01;
		  double u = d+0.01;
		  ForceFields::UFF::DistanceConstraintContrib *contrib = new ForceFields::UFF::DistanceConstraintContrib(field, i, j, l, u, fdist);
		  field->contribs().push_back(ForceFields::ContribPtr(contrib));
	  }

    return field;
  }

  ForceFields::ForceField *construct3DForceField2(const BoundsMatrix &mmat,
                                                 RDGeom::Point3DPtrVect &positions,
                                                 std::vector<std::pair<int, int> > &bonds,
                                                 std::vector<std::pair<int, int> > &angles,
                                                 std::vector<std::vector<int> > &expTorsionAtoms,
                                                 std::vector<std::pair<double, double> > &expTorsionAngles,
                                                 double basinSizeTol) {
      unsigned int N = mmat.numRows();
      CHECK_INVARIANT(N == positions.size(), "");
      CHECK_INVARIANT(expTorsionAtoms.size() == expTorsionAngles.size(), "");
      ForceFields::ForceField *field = new ForceFields::ForceField(positions[0]->dimension());
      for (unsigned int i = 0; i < N; ++i){
        field->positions().push_back(positions[i]);
      }

      // torsion constraints
      double ftorsion = 1000.0; // force constant
      for (unsigned int t = 0; t < expTorsionAtoms.size(); ++t) {
      	int ang1 = expTorsionAngles[t].first;
      	int ang2 = expTorsionAngles[t].second;
      	int i = expTorsionAtoms[t][0];
      	int j = expTorsionAtoms[t][1];
      	int k = expTorsionAtoms[t][2];
      	int l = expTorsionAtoms[t][3];
      	ForceFields::UFF::TorsionConstraintContrib *contrib = new ForceFields::UFF::TorsionConstraintContrib(field, i, j, k, l, ang1, ang2, ftorsion);
      	field->contribs().push_back(ForceFields::ContribPtr(contrib));
      } // torsion constraints


      double fdist = 100.0; // force constant
      // 1,2 distance constraints
      std::vector<std::pair<int, int> >::iterator bi;
      for (bi = bonds.begin(); bi != bonds.end(); bi++) {
    	  unsigned int i = bi->first;
    	  unsigned int j = bi->second;
    	  double d = ((*positions[i]) - (*positions[j])).length();
		  double l = d-0.01;
		  double u = d+0.01;
		  ForceFields::UFF::DistanceConstraintContrib *contrib = new ForceFields::UFF::DistanceConstraintContrib(field, i, j, l, u, fdist);
		  field->contribs().push_back(ForceFields::ContribPtr(contrib));
      }

      // 1,3 distance constraints
		for (bi = angles.begin(); bi != angles.end(); bi++) {
		  unsigned int i = bi->first;
		  unsigned int j = bi->second;
		  double d = ((*positions[i]) - (*positions[j])).length();
		  double l = d-0.01;
		  double u = d+0.01;
		  ForceFields::UFF::DistanceConstraintContrib *contrib = new ForceFields::UFF::DistanceConstraintContrib(field, i, j, l, u, fdist);
		  field->contribs().push_back(ForceFields::ContribPtr(contrib));
	  }

  	  /*for (unsigned int i = 0; i < N-2; ++i) {
	  	  for (unsigned int j = i+1; j <= i+2; ++j) { // only add 1-2 and 1-3 distances
			double d = ((*positions[i]) - (*positions[j])).length();
			double l = d-0.01;
			double u = d+0.01;
			bool includeIt = false;
			if (u-l <= basinSizeTol) {
			  includeIt = true;
			}
			if (includeIt) {
			  ForceFields::UFF::DistanceConstraintContrib *contrib = new ForceFields::UFF::DistanceConstraintContrib(field, i, j, l, u, fdist);
			  field->contribs().push_back(ForceFields::ContribPtr(contrib));
			}
		  }
  	  } // distance constraints*/

      return field;
    }

  ForceFields::ForceField *construct3DForceField3(const BoundsMatrix &mmat,
                                               RDGeom::Point3DPtrVect &positions,
                                               double basinSizeTol) {
    unsigned int N = mmat.numRows();
    CHECK_INVARIANT(N == positions.size(), "");
    ForceFields::ForceField *field = new ForceFields::ForceField(positions[0]->dimension());
    for (unsigned int i = 0; i < N; ++i){
      field->positions().push_back(positions[i]);
    }

    double fdist = 100.0; // force constant
	for (unsigned int i = 1; i < N; ++i) {
	  	for (unsigned int j = 0; j < i; ++j) {
	  		double l = mmat.getLowerBound(i,j);
	  		double u = mmat.getUpperBound(i,j);
			bool includeIt = false;
			if (u-l <= basinSizeTol) {
			  includeIt = true;
			}
			if (includeIt) {
			  ForceFields::UFF::DistanceConstraintContrib *contrib = new ForceFields::UFF::DistanceConstraintContrib(field, i, j, l, u, fdist);
			  field->contribs().push_back(ForceFields::ContribPtr(contrib));
			}
		  }
	} // distance constraints

    return field;
  }
}
