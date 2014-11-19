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
#ifndef __RD_MULTI_RANGE_SQUARE_MATRIX_H__
#define __RD_MULTI_RANGE_SQUARE_MATRIX_H__

#include "Matrix.h"

namespace RDNumeric {
  template <typename TYPE> class MultiRangeSquareMatrix : public Matrix<TYPE> {
  public:
    //! brief Square matrix of size N
    MultiRangeSquareMatrix() {};

    explicit MultiRangeSquareMatrix(unsigned int N) : Matrix<TYPE>(N, N) {};
          
    //MultiRangeSquareMatrix(unsigned int N, TYPE val) : Matrix<TYPE>(N, N, val) {};

    MultiRangeSquareMatrix(unsigned int N, typename Matrix<TYPE>::DATA_SPTR data) : Matrix<TYPE>(N, N, data) {};

  };

}

#endif
    
