/*------------------------------------------------------------------------
Some functions required for kernel processing

Copyright 2020 Grigorii Trofimiuk, Peter Trifonov

Licensed under the Apache License,
Version 2.0(the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
--------------------------------------------------------------------------*/
#pragma once

#include "misc.hpp"
#include <cstdint>

void setPairEW(LLRType y, LLRType &a0, LLRType &a1);


/// compute Q(a, b) = sign(a)*sign(b)*min(a,b)
void SoftXOR2Elements(LLRType &pDest, const LLRType &a, const LLRType &b);
/// compute Q(a, b) = sign(a)*sign(b)*min(a,b) of an array
void SoftXOR(LLRType *pDest, const LLRType *pS, const int &N);

/// compute P(a,b,c) = (-1)^u*a+b
void SoftCombine2Elements(LLRType &pDest, const LLRType &a, const LLRType &b, const tBit &u);
/// compute P(a,b,c) = (-1)^u*a+b of an array
void SoftCombine(LLRType *pDest, const LLRType *pS, const tBit *pC, const int &N);

/// compute FHT of length 4 vector
void hadamardTransform4(const LLRType *in, LLRType *out);

/// inplace calculation of xF^{\otimes m}, where F is the Arikan kernel
template <class T>
void Arikan(unsigned LogLength, T *pData)
{
	unsigned N = 1u << LogLength;
	// implement Arikan's butterfly
	unsigned L = 1;
	while (N > 1) {
		N >>= 1;
		for (unsigned k = 0; k < L; k++) {
			for (unsigned j = 2 * k * N; j < (2 * k + 1) * N; j++) {
				pData[j] ^= pData[j + N];
			};
		};
		L <<= 1;
	};
};
