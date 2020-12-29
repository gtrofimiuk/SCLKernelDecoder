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
#include "Utils.hpp"
#include <cstring>
#include <math.h>

LLRType sign(LLRType a) { return (a <= 0) ? -1 : 1; }

void setPairEW(LLRType y, LLRType &a0, LLRType &a1)
{
	if (y >= 0)
		a1 -= y;
	else
		a0 += y;
	return;
}

// Calculation of P function:
// P(a,b,c) = (-1)^c*a+b
void SoftCombine(LLRType *pDest, const LLRType *pS, const tBit *pC, const int &N)
{
	for (uint32_t beta = 0; beta < N; ++beta) SoftCombine2Elements(pDest[beta], pS[beta], pS[beta + N], pC[beta]);
}

// Calculation of Q function:
// Q(a, b) = sign(a)*sign(b)*min(a,b)
void SoftXOR(LLRType *pDest, const LLRType *pS, const int &N)
{
	for (uint32_t beta = 0; beta < N; beta++) SoftXOR2Elements(pDest[beta], pS[beta], pS[beta + N]);
}

// Calculation of P function:
// P(a,b,c) = (-1)^c*a+b
void SoftCombine2Elements(LLRType &pDest, const LLRType &a, const LLRType &b, const tBit &u)
{
	pDest = (u != 0) ? (b - a) : (b + a);
	SUM_COUNT;
}

// Calculation of Q function:
// Q(a, b) = sign(a)*sign(b)*min(a,b)
void SoftXOR2Elements(LLRType &pDest, const LLRType &a, const LLRType &b)
{
	pDest = sign(a) * sign(b) * std::min(abs(a), abs(b));
	CMP_COUNT;
}

void hadamardTransform4(const LLRType *in, LLRType *out)
{
	LLRType D01 = in[0] - in[1];
	LLRType S01 = in[0] + in[1];
	LLRType D23 = in[2] - in[3];
	LLRType S23 = in[2] + in[3];
	out[0] = D01 - D23;
	out[1] = S01 - S23;
	out[2] = D01 + D23;
	out[3] = S01 + S23;
}
