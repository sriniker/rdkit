//
//  Copyright (c) 2014, Novartis Institutes for BioMedical Research Inc.
// All rights reserved.
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
#ifndef __RD_MULTI_RANGE_BOUNDS_MATRIX_H__
#define __RD_MULTI_RANGE_BOUNDS_MATRIX_H__

#include <RDGeneral/Invariant.h>
#include <boost/smart_ptr.hpp>
#include <iostream>
#include <iomanip>
#include <vector>
#include <stdlib.h>
#include <Numerics/MultiRangeSquareMatrix.h>
#include "BoundsMatrix.h"

namespace DistGeom {

  class vect4mat {
    private:
      std::vector<double> ranges;
    public:
      vect4mat() {}; // Default constructor
      vect4mat(double val) {
        ranges.push_back(val);
      }
      vect4mat(std::vector<double> v) {
        ranges = v;
      }
      std::vector<double> getVector() {
        std::vector<double> newV(ranges.size());
        for (unsigned int i = 0; i < ranges.size(); ++i) {
          newV[i] = ranges[i];
        }
        return newV;
      }
      void push_back(double val) {
        ranges.push_back(val);
      }
      unsigned int size() {
        return ranges.size();
      }
      bool empty() {
        return ranges.empty();
      }
      vect4mat & operator = (double val) {
        ranges.clear();
        ranges.push_back(val);
        return *this;
      }

      vect4mat & operator = (vect4mat &vm) {
        ranges = vm.getVector();
        return *this;
      }
      vect4mat & operator *= (double val) {
        // operation not supported
      }
      vect4mat & operator *= (const vect4mat &vm) {
        // operation not supported
      }
      vect4mat & operator * (double val) {
        // operation not supported
      }
      vect4mat & operator * (const vect4mat &vm) {
        // not supported
      }
      vect4mat & operator /= (const vect4mat &vm) {
        // operation not supported
      }
      vect4mat & operator -= (const vect4mat &vm) {
        // operation not supported
      }
      vect4mat & operator += (const vect4mat &vm) {
        // operation not supported
      }
  };

  //typedef std::vector<double> vect4mat;

  //! Class to store the multi-range distance bound
  /*! 
    Basically a N by N matrix of std::vectors
    with lower distance bounds on the lower triangle and upper bounds in the upper
    triangle
  */
  class MultiRangeBoundsMatrix : public RDNumeric::MultiRangeSquareMatrix<vect4mat> {
  public:
    typedef boost::shared_array<vect4mat> DATA_SPTR;

    explicit MultiRangeBoundsMatrix(unsigned int N) : RDNumeric::MultiRangeSquareMatrix<vect4mat>(N) {};
    MultiRangeBoundsMatrix(unsigned int N, DATA_SPTR data) :
      RDNumeric::MultiRangeSquareMatrix<vect4mat>(N,data) {};

    //! Get the upper bound between points i and j
    inline std::vector<double> getUpperBound(unsigned int i, unsigned int j) const {
      RANGE_CHECK(0, i, d_nRows-1);
      RANGE_CHECK(0, j, d_nCols-1);

      if (i < j) {
        return getVal(i,j).getVector();
      } else {
        return getVal(j,i).getVector();
      }
    }
    
    //! Set the lower bound between points i and j
    inline void setUpperBound(unsigned int i, unsigned int j, double val) {
      RANGE_CHECK(0, i, d_nRows-1);
      RANGE_CHECK(0, j, d_nCols-1);
      CHECK_INVARIANT(val >= 0.0, "Negative upper bound");
      if (i < j) {
        setVal(i,j,val);
      } else {
        setVal(j,i,val);
      }
    }

    //! Set the lower bound between points i and j 
    inline void setLowerBound(unsigned int i, unsigned int j, double val) {
      RANGE_CHECK(0, i, d_nRows-1);
      RANGE_CHECK(0, j, d_nCols-1);
      CHECK_INVARIANT(val >= 0.0, "Negative lower bound");
      if (i < j) {
        setVal(j,i,val);
      } else {
        setVal(i,j,val);
      }
    }
    
    
    //! Get the lower bound between points i and j
    inline std::vector<double> getLowerBound(unsigned int i, unsigned int j) const {
      RANGE_CHECK(0, i, d_nRows-1);
      RANGE_CHECK(0, j, d_nCols-1);

      if (i < j) {
        return getVal(j,i).getVector();
      } else {
        return getVal(i,j).getVector();
      }
    }

    //! Get a random "normal" BoundsMatrix
    void getRandomBoundsMatrix(DistGeom::BoundsMatPtr mmat, std::map<std::pair<int, int>, int > &ranges) {
      unsigned int i, j;
      for (i = 1; i < d_nRows; i++) {
        for (j = 0; j < i; j++) {
          std::vector<double> ub = getUpperBound(i,j);
          if (ub.size() > 1) {
            // choose random range
            int n = rand() % ub.size();
            mmat->setUpperBound(i,j,ub[n]);
            mmat->setLowerBound(i,j,getLowerBound(i,j)[n]);
            if (i < j) { // insert with order
              ranges.insert(std::make_pair(std::make_pair(i, j), n));
            } else {
              ranges.insert(std::make_pair(std::make_pair(j, i), n));
            }
          } else { // only 1 range
            mmat->setUpperBound(i,j,ub[0]);
            mmat->setLowerBound(i,j,getLowerBound(i,j)[0]);
            if (i < j) { // insert with order
			  ranges.insert(std::make_pair(std::make_pair(i, j), 0));
			} else {
			  ranges.insert(std::make_pair(std::make_pair(j, i), 0));
			}
          }
        }
      }
    }

    //! Copy bounds from a "normal" BoundsMatrix
    void copyFromBoundsMatrix(DistGeom::BoundsMatPtr mmat) {
      unsigned int i, j;
      for (i = 1; i < d_nRows; i++) {
        for (j = 0; j < i; j++) {
          if (getUpperBound(i,j).empty()) {
            setUpperBound(i, j, mmat->getUpperBound(i,j));
            setLowerBound(i, j, mmat->getLowerBound(i,j));
          }
        }
      }
    }

    //! sets a particular element of the matrix
    inline void setVal(unsigned int i, unsigned int j, double val) {
      PRECONDITION(i<d_nRows,"bad index");
      PRECONDITION(j<d_nCols,"bad index");
      unsigned int id = i*d_nCols + j;

      d_data[id].push_back(val);
    }

  }; 

  typedef boost::shared_ptr<MultiRangeBoundsMatrix> MultiRangeBoundsMatPtr;
}

#endif
          
