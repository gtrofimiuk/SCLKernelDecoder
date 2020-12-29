/*------------------------------------------------------------------------
Implementation of K_32 kernel from
"Window Processing of Binary Polarization Kernels" by Grigorii Trofimiuk and Peter Trifonov
Now is avaliable at https://arxiv.org/abs/2010.07884

Copyright 2020 Grigorii Trofimiuk

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

#include "KernelProcessor.hpp"

class KernelProcessorTrofimiuk32_342 : public KernelProcessor
{
private:
	void ProcessU0();
	void ProcessU1();
	void ProcessU2();
	void ProcessU4();
	void ProcessU5();
	void ProcessU6();
	void ProcessU7();
	void ProcessU8();
	void ProcessU9();
	void ProcessU10();
	void ProcessU11();
	void ProcessU12();
	void ProcessU13();
	void ProcessU16();
	void ProcessU17();
	void ProcessU20();
	void ProcessU27();

	// Described in Section V-B
	void RunRecursiveMaximumComputation(unsigned depth);

	// Expression (20)
	void UpdatePathScores(unsigned numOfPathScores);

	// Expression (20), but indexing depends on i
	void UpdatePathScoresWithSign(unsigned numOfPathScores, unsigned signMask);

	// See Section V-C
	void ComputeReusedMaximum(unsigned numOfPathScores, unsigned sign, unsigned prevMaxBit = 0);

	/// Array for Arikan SC intermediate LLRs
	LLRType *pS[6];

	//temporary arrays, used for convenience
	LLRType* tmpLLR;
	tBit c[32];

	/// CSE arrays, see Section IV
	/// One dimensional arrays can be used, but two-dimensional representation is more convenient
	LLRType* X1;     // 32 values
	LLRType* X2[8];  // 8 * 2 values
	LLRType* X3[4];  // 4 * 2 values
	LLRType* X4[2];  // 2 * 8 values 
	LLRType* X5;

	//path scores of Arikan SC decoding
	LLRType* EW; //32 values

	/// We need these values to reuse maximums from previous phases, see Section V-C
	LLRType* maxEW;
	unsigned maxCW[2];

	/// This array is used for the recursive path score maximization, see Section V-B
	LLRType *m_pMaxStorage[5];
	unsigned shift;

	static const unsigned StateSize = 120;
	LLRType pState[StateSize];

	/// Input vector of Arikan matrix F_4, see expression (7)
	tBit v[32];

	LLRType pInputLLRs[32];

public:
	virtual tBit *GetCodeword();

	virtual void ProcessNextSymbol();

	virtual void CopyState(KernelProcessor &dst);

	KernelProcessorTrofimiuk32_342();
	~KernelProcessorTrofimiuk32_342(){};
};
