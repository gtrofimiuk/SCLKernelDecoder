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

#include "KernelProcessorTrofimiuk32_342.hpp"
#include "KernelProcessorTrofimiuk32_342_CSE.hpp"
#include "Utils.hpp"
#include <immintrin.h>
#include <iostream>
#include <math.h>

//Described in Section V-B
void KernelProcessorTrofimiuk32_342::RunRecursiveMaximumComputation(unsigned depth)
{
	LLRType *pPreviousStorage = EW;
	// we save intermediate comparisons to reuse it
	for (int i = depth; i >= 0; --i) {
		LLRType *pCurrentStorage = m_pMaxStorage[i];
		unsigned range = 1u << (i + 1);
		for (unsigned j = 0; j < range; ++j) {
			pCurrentStorage[j] = std::max(pPreviousStorage[2 * j], pPreviousStorage[2 * j + 1]);
			CMP_COUNT;
		}
		pPreviousStorage = pCurrentStorage;
	}
	m_CurrentLLR = m_pMaxStorage[0][0] - m_pMaxStorage[0][1];
	SUM_COUNT;
}

//Expression (20)
void KernelProcessorTrofimiuk32_342::UpdatePathScores(unsigned numOfPathScores)
{
	for (unsigned i = 0; i < numOfPathScores; ++i) {
		EW[i + numOfPathScores] = EW[i];
		setPairEW(X5[i], EW[i], EW[i + numOfPathScores]);
	}
	SUM_BLOCK_COUNT(numOfPathScores);
}

// Expression (20), but indexing depends on i
void KernelProcessorTrofimiuk32_342::UpdatePathScoresWithSign(unsigned numOfPathScores, unsigned signMask)
{
	for (unsigned i = 0; i < numOfPathScores; ++i) {
		unsigned sign = (_mm_popcnt_u32(i & signMask) % 2);
		unsigned Index0 = i + numOfPathScores * sign, Index1 = i + numOfPathScores * (sign ^ 1);
		EW[(m_CurrentC ^ 1) * numOfPathScores + i] = EW[m_CurrentC * numOfPathScores + i];
		setPairEW(X5[i], EW[Index0], EW[Index1]);
	}
	SUM_BLOCK_COUNT(numOfPathScores);
}

//See Section V-C
void KernelProcessorTrofimiuk32_342::ComputeReusedMaximum(unsigned numOfPathScores, unsigned a, unsigned prevMaxBit)
{
	maxCW[a] = maxCW[prevMaxBit];
	maxEW[a] = EW[maxCW[a] + numOfPathScores * a];

	unsigned offset = numOfPathScores * (a ^ 1);
	maxEW[a ^ 1] = EW[offset], maxCW[a ^ 1] = 0;
	for (unsigned i = 1; i < numOfPathScores; ++i) {
		if (EW[i + offset] > maxEW[a ^ 1]) {
			maxEW[a ^ 1] = EW[i + offset], maxCW[a ^ 1] = i;
		}
	}
	CMP_BLOCK_COUNT(numOfPathScores - 1);

	m_CurrentLLR = maxEW[0] - maxEW[1];
	SUM_COUNT;
}

KernelProcessorTrofimiuk32_342::KernelProcessorTrofimiuk32_342()
{
	m_LastCalledPhase = -1;
	m_KernelSize = 32;

	m_pLLR = pInputLLRs;

	// We pack state of the processor in a single array pState
	// Here we allocate memory for CSE
	// We do not minimized the size of pState, we suppose that it can be less
	X1 = pState; //32 values
	tmpLLR = pState + 32;
	maxEW = pState + 72;
    for (unsigned i = 0; i < 8; ++i) X2[i] = pState + 32 + 2 * i; //16 values
	for (unsigned i = 0; i < 4; ++i) X3[i] = pState + 48 + 2 * i; //8 values
	X4[0] = pState + 56, X4[1] = pState + 64;
	X5 = pState + 72;
	EW = pState + 88; //32 values
	pS[1] = pState + 56;
	pS[2] = pState + 96; //do not ask me why
	pS[3] = pState + 104;
	pS[4] = pState + 108;
	pS[5] = pState + 110;
	m_pMaxStorage[0] = pState + 72; //2 values
	m_pMaxStorage[1] = pState + 74; //4 values
	m_pMaxStorage[2] = pState + 78; //8 values
	m_pMaxStorage[3] = pState + 56; //16 values
}

tBit *KernelProcessorTrofimiuk32_342::GetCodeword()
{
	v[31] = m_CurrentC;

	memcpy(c, v, sizeof(tBit) * 32);
	Arikan(5, c);

	return c;
};
/*--- These phases are equivalent to min-sum Arikan SC decoding ---*/
void KernelProcessorTrofimiuk32_342::ProcessU0()
{
	memset(v, 0, sizeof(tBit) * m_KernelSize);

	SoftXOR(pS[1], m_pLLR, 16);
	SoftXOR(pS[2], pS[1], 8);
	SoftXOR(pS[3], pS[2], 4);
	SoftXOR(pS[4], pS[3], 2);
	SoftXOR(pS[5], pS[4], 1);

	m_CurrentLLR = pS[5][0];
}
void KernelProcessorTrofimiuk32_342::ProcessU1()
{
	SoftCombine(pS[5], pS[4], &m_CurrentC, 1);

	m_CurrentLLR = pS[5][0];
}
void KernelProcessorTrofimiuk32_342::ProcessU2()
{
	SoftCombine(pS[4], pS[3], c, 2);
	SoftXOR(pS[5], pS[4], 1);

	m_CurrentLLR = pS[5][0];
}
void KernelProcessorTrofimiuk32_342::ProcessU4()
{
	SoftCombine(pS[3], pS[2], c, 4);
	SoftXOR(pS[4], pS[3], 2);
	SoftXOR(pS[5], pS[4], 1);

	m_CurrentLLR = pS[5][0];
}
/*-----------------------------------------------------------------*/

void KernelProcessorTrofimiuk32_342::ProcessU5()
{
	// we had y0,y1,...,y7 LLRs. if v4 = 1, then signs of y0 and y4 should be inverted
	// moreover, we can skip the summation of LLRs, which has already been done in
	// SoftCombine(m_pS[3], m_pS[2], c, 4) of ProcessU4()

	pS[3][0] *= 1 - 2 * m_CurrentC;

	LLRType *tR = tmpLLR;
	// tR[0] -> 0110 tR[1] -> 0011 tR[2] -> 0101 tR[3] -> 0000
	hadamardTransform4(pS[3], tR);
	SUM_BLOCK_COUNT(8);

	// We want to order path scores of codewords in indexing v5|v6|v7, where 001 corresponds to v7 = 1
	EW[0] = tR[3], EW[1] = -tR[3], EW[2] = -tR[2], EW[3] = tR[2], EW[4] = -tR[1], EW[5] = tR[1], EW[6] = tR[0], EW[7] = -tR[0];
	for (unsigned i = 0; i < 8; ++i) EW[i] /= 2;

	// here we find a maximal path on the v7 phase
	// the EW of this path will be also amximal among EWs after v8 phase
	LLRType Hmax = abs(EW[0]);
	maxCW[0] = unsigned(EW[0] < 0);
	for (unsigned i = 1; i < 4; ++i) {
		if (abs(EW[2 * i]) > Hmax) Hmax = abs(EW[2 * i]), maxCW[0] = 2 * i ^ unsigned(EW[2 * i] < 0);
	}
	CMP_BLOCK_COUNT(3);
	/*--- Now we move to 8 phase of Arikan SC decoding ---*/

	// In fact, we need to calculate how to use signs of previously known values of v
	// Instead of LLRs, we use values of m_pS[1]

	// bit-reversal permutation
	for (unsigned i = 0; i < 16; ++i) {
		unsigned reversed = 0, mult = 1;
		for (unsigned j = 0; j < 4; ++j) reversed += ((i >> (3 - j) & 1)) << j;
		tmpLLR[i] = pS[1][reversed];
		if (c[reversed] != 0) tmpLLR[i] *= -1;
	}

	// L[0] = y0+y1, L[1] = -y0+y1, L[2] = y2+y3, L[3] = -y2+y3...
	// Here we compute X2 CSE, but it is convenient to do it in X1 array
	LLRType *cX2 = X1;
	for (unsigned i = 0; i < 8; ++i) {
		cX2[2 * i] = tmpLLR[2 * i] + tmpLLR[2 * i + 1];
		cX2[2 * i + 1] = -tmpLLR[2 * i] + tmpLLR[2 * i + 1];
	}
	SUM_BLOCK_COUNT(16);

	for (unsigned i = 0; i < 4; ++i) {
		SoftXOR2Elements(X3[i][0], cX2[4 * i + 0], cX2[4 * i + 2]);  // 00
		SoftXOR2Elements(X3[i][1], cX2[4 * i + 1], cX2[4 * i + 3]);  // 11
	}

	// The structure of CSE allow as to obtain its indices directly from index j
	for (unsigned i = 0; i < 2; ++i) {
		for (unsigned j = 0; j < 4; ++j) {
			unsigned ind0 = ((j >> 1) & 1), ind1 = j & 1;
			SoftXOR2Elements(X4[i][j], X3[2 * i][ind0], X3[2 * i + 1][ind1]);
		}
	}

	for (unsigned i = 0; i < 8; ++i) SoftXOR2Elements(X5[i], X4[0][U8_X5[i][0]], X4[1][U8_X5[i][1]]);

	UpdatePathScores(8);

	unsigned v8 = (X5[maxCW[0]] < 0) ? 1 : 0;
	ComputeReusedMaximum(8, v8);
}

void KernelProcessorTrofimiuk32_342::ProcessU6()
{
	for (unsigned i = 0; i < 8; ++i) SoftCombine2Elements(X5[i], X4[0][U8_X5[i][0]], X4[1][U8_X5[i][1]], m_CurrentC);

	UpdatePathScoresWithSign(8, 6);  // the mask corresponds to v5 ^ v6

	// Again, we already know one max
	// we get v5^v6 of maxCW[m_CurrentC]
	unsigned v56m = ((maxCW[m_CurrentC] >> 2) & 1) ^ ((maxCW[m_CurrentC] >> 1) & 1);
	unsigned signZ = (X5[maxCW[m_CurrentC]] < 0) ? 1 : 0;
	unsigned signM = v56m ^ signZ;

	ComputeReusedMaximum(8, signM, m_CurrentC);
}
void KernelProcessorTrofimiuk32_342::ProcessU7()
{
	// we constructed CSE in assumtion that v8 = v9 = 0
	// in actual processing, in SoftCombine we need to adjust values
	for (unsigned i = 0; i < 2; ++i) {
		for (unsigned j = 0; j < 8; ++j) {
			tBit add = c[i] ^ ((j >= 4) * 1);  // this due our construction of CSE
			SoftCombine2Elements(X4[i][j], X3[2 * i][U10_X4[j][0]], X3[2 * i + 1][U10_X4[j][1]], add);
		}
	}

	for (unsigned i = 0; i < 8; ++i) SoftXOR2Elements(X5[i], X4[0][U10_X5[i][0]], X4[1][U10_X5[i][1]]);

	UpdatePathScoresWithSign(8, 4);  // v5
	RunRecursiveMaximumComputation(2);
}

/*--- These phases are responsible for recursive maximum computation, see Section V-B---*/
void KernelProcessorTrofimiuk32_342::ProcessU8()
{
	shift = m_CurrentC << 1;
	m_CurrentLLR = m_pMaxStorage[1][shift] - m_pMaxStorage[1][shift + 1];
	SUM_COUNT;
}
void KernelProcessorTrofimiuk32_342::ProcessU9()
{
	shift = (shift << 1) + (m_CurrentC << 1);
	m_CurrentLLR = m_pMaxStorage[2][shift] - m_pMaxStorage[2][shift + 1];
	SUM_COUNT;
}
void KernelProcessorTrofimiuk32_342::ProcessU10()
{
	shift = (shift << 1) + (m_CurrentC << 1);
	m_CurrentLLR = EW[shift] - EW[shift + 1];
	SUM_COUNT;
}
/*-------------------------------------------------------------------------------------*/
void KernelProcessorTrofimiuk32_342::ProcessU11()
{
	v[7] = m_CurrentC;
	unsigned pathIndex = v[5] * 4 + v[6] * 2 + v[7];

	// We can backtrack indices of CSE and obtain LLRs for a single path
	SoftCombine2Elements(X5[pathIndex], X4[0][U10_X5[pathIndex][0]], X4[1][U10_X5[pathIndex][1]], v[10]);
	m_CurrentLLR = X5[pathIndex];
}
void KernelProcessorTrofimiuk32_342::ProcessU12()
{
	v[11] = m_CurrentC;

	// Here we have a decoding window given by 12,13,14,15
	// at first, we need to compute path scores
	memset(c, 0, sizeof(tBit) * 8);
	c[5] = v[5], c[6] = v[6], c[7] = v[7];  // we already considered v0,...,v4 in L array

	Arikan(3, c);
	for (unsigned i = 0; i < 8; ++i) pS[2][BR3[i]] = (c[BR3[i]] == 0) ? X1[2 * i] : X1[2 * i + 1];

	memcpy(c, v + 8, sizeof(tBit) * 4);
	Arikan(2, c);
	LLRType y[4];

	// here we use Corollary 1 to obtain path scores
	SoftCombine(pS[3], pS[2], c, 4);

	for (unsigned i = 0; i < 4; ++i) {
		y[i] = -abs(pS[3][i]);
		c[i] = 1 * (pS[3][i] < 0);
	}
	Arikan(2, c);

	// here we compute sum_{i=0}^3 \tau(y_i,c_i), where c \in F_2^4

	unsigned sh = 0;
	for (unsigned i = 0; i < 4; ++i) sh ^= (1u << i) * (c[3 - i] != 0);
	maxCW[0] = sh;

	EW[0 ^ sh] = 0, EW[8 ^ sh] = y[0], EW[12 ^ sh] = y[1], EW[10 ^ sh] = y[2], EW[15 ^ sh] = y[3];
	EW[2 ^ sh] = y[0] + y[2], EW[3 ^ sh] = y[1] + y[3], EW[4 ^ sh] = y[0] + y[1];
	EW[5 ^ sh] = y[2] + y[3], EW[6 ^ sh] = y[1] + y[2], EW[7 ^ sh] = y[0] + y[3];
	EW[9 ^ sh] = EW[5 ^ sh] + y[1], EW[11 ^ sh] = EW[3 ^ sh] + y[0];
	EW[13 ^ sh] = EW[5 ^ sh] + y[0], EW[14 ^ sh] = EW[6 ^ sh] + y[0];
	EW[1 ^ sh] = EW[9 ^ sh] + y[0];
	SUM_BLOCK_COUNT(11);

	// the decoding window is given by 12,13,14,15
	memset(c, 0, sizeof(tBit) * 32);
	memcpy(c, v, sizeof(tBit) * 12);
	Arikan(4, c);
	// bit-reversal permutation, can be omitted, but we use it for convenience
	for (unsigned i = 0; i < 32; ++i) {
		unsigned reversed = 0, mult = 1;
		for (unsigned j = 0; j < 5; ++j) reversed += ((i >> (4 - j) & 1)) << j;
		tmpLLR[i] = m_pLLR[reversed];
		if (c[reversed] != 0) tmpLLR[i] *= -1;
	}

	for (unsigned i = 0; i < 16; ++i) {
		X1[2 * i] = tmpLLR[2 * i] + tmpLLR[2 * i + 1];
		X1[2 * i + 1] = -tmpLLR[2 * i] + tmpLLR[2 * i + 1];
	}
	SUM_BLOCK_COUNT(32);

	for (unsigned i = 0; i < 8; ++i) {
		SoftXOR2Elements(X2[i][0], X1[4 * i + 0], X1[4 * i + 2]);  // 00
		SoftXOR2Elements(X2[i][1], X1[4 * i + 1], X1[4 * i + 3]);  // 11
	}

	for (unsigned i = 0; i < 4; ++i) {
		SoftXOR2Elements(X3[i][0], X2[2 * i][0], X2[2 * i + 1][0]);  // 0000
		SoftXOR2Elements(X3[i][1], X2[2 * i][1], X2[2 * i + 1][1]);  // 1111
	}

	for (unsigned i = 0; i < 2; ++i) {
		for (unsigned j = 0; j < 4; ++j) {
			unsigned ind0 = ((j >> 1) & 1), ind1 = j & 1;
			SoftXOR2Elements(X4[i][j], X3[2 * i][ind0], X3[2 * i + 1][ind1]);
		}
	}

	for (unsigned i = 0; i < 16; ++i) SoftXOR2Elements(X5[i], X4[0][U16_X5[i][0]], X4[1][U16_X5[i][1]]);

	UpdatePathScores(16);

	unsigned v16 = (X5[maxCW[0]] < 0) ? 1 : 0;
	ComputeReusedMaximum(16, v16);
}
void KernelProcessorTrofimiuk32_342::ProcessU13()
{
	v[16] = m_CurrentC;

	for (unsigned i = 0; i < 16; ++i) SoftCombine2Elements(X5[i], X4[0][U16_X5[i][0]], X4[1][U16_X5[i][1]], v[16]);

	UpdatePathScoresWithSign(16, 8);  // v12
	RunRecursiveMaximumComputation(3);
}
void KernelProcessorTrofimiuk32_342::ProcessU16()
{
	v[13] = m_CurrentC;

	/*--- obtain path scores from previous phase ---*/
	// here we have the decoding window reduced from 4 to 2
	// {12, 13, 14, 15} -> {14, 15}
	unsigned prevShift = shift;
	shift = (shift << 1) + (v[13] << 1);
	shift <<= 1;

	// This copying is done just for convenience
	if (shift > 0) memcpy(EW, EW + shift, sizeof(LLRType) * 4);

	maxEW[0] = EW[0];
	maxCW[0] = 0;
	// we also know maximum for this array, but we need the corresponding index
	for (unsigned i = 1; i < 4; ++i)
		if (EW[i] == m_pMaxStorage[2][prevShift + v[13]]) maxCW[0] = i;

	/*-------------------------------------------------*/

	// we need to chose those ones in cse, which corresponds to coset given by v12 and v13
	shift = (10 * v[13]) ^ (8 * v[12]);  // bit magick, induced by CSE structure
	c[0] = v[16] ^ v[17], c[1] = v[17];

	for (unsigned i = 0; i < 2; ++i) {
		unsigned sh = (shift >> 2 * (1 - i)) & 3;
		unsigned ind0 = sh >> 1, ind1 = sh & 1;
		SoftCombine2Elements(X4[i][0], X3[2 * i][ind0], X3[2 * i + 1][ind1], c[i]);
		SoftCombine2Elements(X4[i][1], X3[2 * i][ind0 ^ 1], X3[2 * i + 1][ind1 ^ 1], c[i]);
	}

	for (unsigned i = 0; i < 4; ++i) SoftXOR2Elements(X5[i], X4[0][U18_X5[i][0]], X4[1][U18_X5[i][1]]);

	UpdatePathScores(4);

	unsigned v18 = (X5[maxCW[0]] < 0) ? 1 : 0;
	ComputeReusedMaximum(4, v18);
}
void KernelProcessorTrofimiuk32_342::ProcessU17()
{
	v[18] = m_CurrentC;
	unsigned len = 4;  // The number of path scores to be computed = 2^|D_{17}|

	for (unsigned i = 0; i < len; ++i) SoftCombine2Elements(X5[i], X4[0][U18_X5[i][0]], X4[1][U18_X5[i][1]], v[18]);

	UpdatePathScoresWithSign(4, 2);  // v14
	RunRecursiveMaximumComputation(1);
}
void KernelProcessorTrofimiuk32_342::ProcessU20()
{
	v[15] = m_CurrentC;

	unsigned pathIndex = v[12] * 8 + v[13] * 4 + v[14] * 2 + v[15];

	// here we retrieve the intermediate LLRs for a single path in Arikan SC decoding
	unsigned ZInd[2], YInd[4], XInd[8];
	ZInd[0] = U16_X5[pathIndex][0], ZInd[1] = U16_X5[pathIndex][1];
	for (unsigned i = 0; i < 2; ++i) YInd[2 * i] = ((ZInd[i] >> 1) & 1), YInd[2 * i + 1] = ZInd[i] & 1;
	for (unsigned i = 0; i < 4; ++i) XInd[2 * i] = XInd[2 * i + 1] = YInd[i];
	for (unsigned i = 0; i < 8; ++i) pS[2][BR3[i]] = X2[i][XInd[i]];

	// for the next phase
	for (unsigned i = 0; i < 8; ++i)
		pS[1][BR4[2 * i]] = X1[4 * i + XInd[i]], pS[1][BR4[2 * i + 1]] = X1[4 * i + 2 + XInd[i]];

	c[0] = v[16], c[1] = v[17], c[2] = v[18], c[3] = v[19];
	Arikan(2, c);
	SoftCombine(pS[3], pS[2], c, 4);
	SoftXOR(pS[4], pS[3], 2);
	SoftXOR(pS[5], pS[4], 1);

	m_CurrentLLR = pS[5][0];
}
void KernelProcessorTrofimiuk32_342::ProcessU27()
{
	v[23] = m_CurrentC;
	unsigned pathIndex = v[21] * 4 + v[22] * 2 + v[23];

	// for the next phase
	unsigned ZInd[2], YInd[4];
	ZInd[0] = U10_X5[pathIndex][0], ZInd[1] = U10_X5[pathIndex][1];
	for (unsigned i = 0; i < 2; ++i) YInd[2 * i] = ((ZInd[i] >> 1) & 1), YInd[2 * i + 1] = ZInd[i] & 1;
	for (unsigned i = 0; i < 4; ++i)
		pS[2][BR3[2 * i]] = X1[4 * i + YInd[i]], pS[2][BR3[2 * i + 1]] = X1[4 * i + 2 + YInd[i]];

	// We can backtrack indices of CSE
	SoftCombine2Elements(X5[pathIndex], X4[0][ZInd[0]], X4[1][ZInd[1]], v[26]);

	m_CurrentLLR = X5[pathIndex];
}

// process next required symbol
void KernelProcessorTrofimiuk32_342::ProcessNextSymbol()
{
	unsigned phase = m_LastCalledPhase + 1;

	switch (phase) {
	case 0: ProcessU0(); break;
	case 1:
		v[0] = m_CurrentC;
		ProcessU1();
		break;
	case 2:
		v[1] = m_CurrentC;
		c[0] = v[0] ^ v[1], c[1] = v[1];
		ProcessU2();
		break;
	case 3:
		// phase 3 is equivalent to phase 1 due to structure of Arikan SC decoding
		v[2] = m_CurrentC;
		ProcessU1();
		break;
	case 4:
		v[3] = m_CurrentC;
		memcpy(c, v, 4 * sizeof(tBit));

		c[0] = v[0] ^ v[1] ^ v[2] ^ v[3], c[1] = v[1] ^ v[3], c[2] = v[2] ^ v[3], c[3] = v[3];
		ProcessU4();
		break;
	case 5:
		v[4] = m_CurrentC;
		memset(c, 0, sizeof(tBit) * 16);
		c[0] = v[0] ^ v[1] ^ v[2] ^ v[3] ^ v[4], c[1] = v[1] ^ v[3], c[2] = v[2] ^ v[3], c[3] = v[3], c[4] = v[4];
		ProcessU5();
		break;
	case 6:
		v[8] = m_CurrentC;
		ProcessU6();
		break;
	case 7:
		v[9] = m_CurrentC;
		c[0] = v[8] ^ v[9], c[1] = v[9];
		ProcessU7();
		break;
	case 8:
		v[10] = m_CurrentC;  // v5 + v10
		ProcessU8();
		break;
	case 9:
		v[5] = m_CurrentC;
		v[10] ^= v[5];
		ProcessU9();
		break;
	case 10:
		v[6] = m_CurrentC;
		v[9] ^= v[6] ^ v[5];
		ProcessU10();
		break;
	case 11: ProcessU11(); break;
	case 12: ProcessU12(); break;
	case 13: ProcessU13(); break;
	case 14:
		v[17] = m_CurrentC;  // v12 + v17
		ProcessU8();
		break;
	case 15:
		v[12] = m_CurrentC;
		v[17] ^= v[12];
		ProcessU9();
		break;
	case 16: ProcessU16(); break;
	case 17: ProcessU17(); break;
	case 18:
		v[19] = m_CurrentC;  // v14 + v19
		ProcessU8();
		break;
	case 19:
		v[14] = (m_CurrentC);
		v[19] ^= v[14];
		ProcessU10();
		break;
	case 20: ProcessU20(); break;
	case 21:
		v[20] = m_CurrentC;
		memset(c, 0, sizeof(tBit) * 16);
		c[0] = v[16] ^ v[17] ^ v[18] ^ v[19] ^ v[20], c[1] = v[17] ^ v[19], c[2] = v[18] ^ v[19], c[3] = v[19], c[4] = v[20];
		/// In the implementation of this kernel we show how processing algorithms of different phases can be reused.
		/// This is possible due to the structure of the kernel. Namely, (u_5 = v_8, u_{21} = v_{24}), (u_6 = v_5 + v_6
		/// + v_9, u_{22} = v_{21}+v_{22}+v_{25}), and (u_7 = v_5 + v_{10}, u_{23} = v_{21} + v_{26}) has absolutely the
		/// same structure of CSE, which is demonstrated in this code
		ProcessU5();
		break;
	case 22:
		v[24] = m_CurrentC;
		ProcessU6();
		break;
	case 23:
		v[25] = m_CurrentC;
		c[0] = v[24] ^ v[25], c[1] = v[25];
		ProcessU7();
		break;
	case 24:
		v[26] = m_CurrentC;
		ProcessU8();
		break;
	case 25:
		v[21] = m_CurrentC;
		v[26] ^= v[21];
		ProcessU9();
		break;
	case 26:
		v[22] = m_CurrentC;
		v[25] ^= v[22] ^ v[21];
		ProcessU10();
		break;
	case 27: ProcessU27(); break;
	case 28:
		v[27] = m_CurrentC;
		c[0] = v[24], c[1] = v[25], c[2] = v[26], c[3] = v[27];
		Arikan(2, c);
		ProcessU4();
		break;
	case 29:
		v[28] = m_CurrentC;

		ProcessU1();
		break;
	case 30:
		v[29] = m_CurrentC;

		c[0] = v[28] ^ v[29], c[1] = v[29];
		ProcessU2();
		break;
	case 31:
		v[30] = m_CurrentC;
		ProcessU1();
		break;
	default: break;
	}
	m_LastCalledPhase = phase;
};

void KernelProcessorTrofimiuk32_342::CopyState(KernelProcessor &srca)
{
	KernelProcessorTrofimiuk32_342 &src = (KernelProcessorTrofimiuk32_342 &)srca;

	m_CurrentC = src.m_CurrentC;
	m_LastCalledPhase = src.m_LastCalledPhase;

	// copy input LLR
	// Can be done without copying by sharing the pointers
	std::copy_n(src.m_pLLR, m_KernelSize, m_pLLR);

	// copy state variable, contains LLRType values in kernel processor
	// Actually, obe can share the pointers between kernel processors.
	// Moreover, we do not need all state variables for each phase.
	// However, we did not implemented this yet
	std::copy_n(src.pState, StateSize, pState);

	shift = src.shift;

	maxCW[0] = src.maxCW[0], maxCW[1] = src.maxCW[1];

	std::copy_n(src.v, 32, v);
}
