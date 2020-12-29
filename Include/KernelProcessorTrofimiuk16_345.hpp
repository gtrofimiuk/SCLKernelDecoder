/*------------------------------------------------------------------------
Implementation of K_16 kernel from 
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

class KernelProcessorTrofimiuk16_345 : public KernelProcessor
{
private:
	void ProcessU0();
	void ProcessU1();
	void ProcessU2();
	void ProcessU3();
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
	void ProcessU14();
	void ProcessU15();

	/// Array for Arikan SC intermediate LLRs
	LLRType *pS[5];

	/// Partial sums
	tBit pC[16];

	/// Input vector of Arikan matrix F_4, see expression (7)
	tBit v[16];

	/// Temporary array
	LLRType *pTmp;  // 16 values

	/// CSE arrays, see Section IV
	/// One dimensional arrays can be used, but two-dimensional representation is more convenient
	LLRType *X1;   // 16 values
	LLRType *X2[4];  // 4 * 2 values
	LLRType *X3[2];  // 2 * 8 values
	LLRType *X4;   // 8 values

	/// Path scores R_y, see expression (20)
	LLRType *pPathScores;  // 16 values
	
	/// This array is used for the recursive path score maximization, see Section V-B
	LLRType *pMaxStorage[4];
	unsigned shift;

	/// We need these values to reuse maximums from previous phases, see Section V-C
	LLRType *maxEW;  // 2 values
	unsigned maxCW[2];

	static const unsigned StateSize = 78;
	LLRType pState[StateSize];

	LLRType InputLLRs[16];
public:
	virtual tBit *GetCodeword();

	virtual void ProcessNextSymbol();

	virtual void CopyState(KernelProcessor &dst);

	KernelProcessorTrofimiuk16_345();
	~KernelProcessorTrofimiuk16_345(){};
};