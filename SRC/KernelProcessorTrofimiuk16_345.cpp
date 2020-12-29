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
#include "KernelProcessorTrofimiuk16_345.hpp"
#include "Utils.hpp"

#define UPPER_BOUND 1e+25

KernelProcessorTrofimiuk16_345::KernelProcessorTrofimiuk16_345()
{
	m_LastCalledPhase = -1;
	m_KernelSize = 16;
	m_pLLR = InputLLRs;

	// We pack state of the processor in a single array pState
	// Here we allocate memory for CSE
	X1 = pState;
	pS[1] = pState;
	pS[2] = pState + 16;
	X2[0] = pState + 16, X2[1] = pState + 16 + 2, X2[2] = pState + 16 + 4, X2[3] = pState + 16 + 6;
	pS[3] = pState + 24;
	X3[0] = pState + 24, X3[1] = pState + 32;
	pS[4] = pState + 40;
	X4 = pState + 40;
	pPathScores = pState + 48;

	pMaxStorage[0] = pState + 64;
	maxEW = pState + 64;
	pMaxStorage[1] = pState + 66;
	pMaxStorage[2] = pState + 70;
	// pState + 80

	pTmp = pState + 16;
}

// Phase 0, just Arikan SC decoding
void KernelProcessorTrofimiuk16_345::ProcessU0()
{
	// Expressions (8)-(11)
	SoftXOR(pS[1], m_pLLR, 8);
	SoftXOR(pS[2], pS[1], 4);
	SoftXOR(pS[3], pS[2], 2);
	SoftXOR(pS[4], pS[3], 1);
	m_CurrentLLR = pS[4][0];
}

// Phase 1, just Arikan SC decoding
void KernelProcessorTrofimiuk16_345::ProcessU1()
{
	v[0] = m_CurrentC;

	SoftCombine(pS[4], pS[3], &m_CurrentC, 1);
	m_CurrentLLR = pS[4][0];
}

// Phase 2, just Arikan SC decoding
void KernelProcessorTrofimiuk16_345::ProcessU2()
{
	v[1] = m_CurrentC;

	// Compute partial sums
	pC[0] = v[0] ^ v[1], pC[1] = v[1];
	SoftCombine(pS[3], pS[2], pC, 2);
	SoftXOR(pS[4], pS[3], 1);
	m_CurrentLLR = pS[4][0];
}

// Phase 3, just Arikan SC decoding
void KernelProcessorTrofimiuk16_345::ProcessU3()
{
	v[2] = m_CurrentC;

	SoftCombine(pS[4], pS[3], &m_CurrentC, 1);
	m_CurrentLLR = pS[4][0];
}

// Phase 4, just Arikan SC decoding
void KernelProcessorTrofimiuk16_345::ProcessU4()
{
	v[3] = m_CurrentC;

	// partial sums
	std::copy_n(v, 4, pC);
	Arikan(2, pC);

	SoftCombine(pS[2], pS[1], pC, 4);
	SoftXOR(pS[3], pS[2], 2);
	SoftXOR(pS[4], pS[3], 1);
	m_CurrentLLR = pS[4][0];
}

// Phase 5, u_5 = v_8, the decoding window = {5, 6, 7}
void KernelProcessorTrofimiuk16_345::ProcessU5()
{
	v[4] = m_CurrentC;

	/*--- Algorithm 3, Step 1. Here we compute 8 path scores, see Section V-A and VI-B for more details ---*/

	// we use values form pS[2] array as the input for FHT, we only need to change LLRs if v[4] == 1
	pS[2][0] *= 1 - 2 * v[4];

	// pTmp[0] -> 0110 pTmp[1] -> 0011 pTmp[2] -> 0101 pTmp[3] -> 0000
	hadamardTransform4(pS[2], pTmp);
	SUM_BLOCK_COUNT(8);

	// Use results of FHT to obtain all correlations
	pPathScores[0] = pTmp[3], pPathScores[1] = -pTmp[3], pPathScores[2] = -pTmp[2], pPathScores[3] = pTmp[2];
	pPathScores[4] = -pTmp[1], pPathScores[5] = pTmp[1], pPathScores[6] = pTmp[0], pPathScores[7] = -pTmp[0];
	for (unsigned i = 0; i < 8; ++i) pPathScores[i] /= 2;

	// here we find a maximal path on the internal phase 7
	// the ellisodal weight (EW) of this path will be also amximal among EWs after internal phase 8
	LLRType Hmax = abs(pPathScores[0]);
	maxCW[0] = unsigned(pPathScores[0] < 0);
	for (unsigned i = 1; i < 4; ++i) {
		if (abs(pPathScores[2 * i]) > Hmax) {
			Hmax = abs(pPathScores[2 * i]);
			maxCW[0] = 2 * i ^ unsigned(pPathScores[2 * i] < 0);
		}
	}
	CMP_BLOCK_COUNT(3);
	/*--- The end of Algorithm 3, Step 1 ---*/

	/*--- Alogrithm 3, Step 2. Use CSE to compute LLRs, see Algorithm 2 --- */
	// All CSE were identified offline, here we use just integer indices

	// One can switch signs of kernel LLRs or change c in SoftCombine function, here we choose the first
	std::copy_n(m_pLLR, 16, pTmp);

	pC[0] = v[0] ^ v[1] ^ v[2] ^ v[3] ^ v[4], pC[1] = v[1] ^ v[3], pC[2] = v[2] ^ v[3], pC[3] = v[3], pC[4] = v[4];

	for (unsigned i = 0; i < 5; ++i) pTmp[i] *= 1 - 2 * pC[i];

	for (unsigned i = 0; i < 8; ++i) {
		X1[i] = pTmp[8 + i] + pTmp[i];
		X1[8 + i] = pTmp[8 + i] - pTmp[i];
	}
	SUM_BLOCK_COUNT(16);

	// fill the X2 array
	for (unsigned i = 0; i < 4; ++i) {
		SoftXOR2Elements(X2[i][0], X1[i], X1[i + 4]);
		SoftXOR2Elements(X2[i][1], X1[i + 8], X1[i + 12]);
	}

	// fill the X3 array
	for (unsigned i = 0; i < 2; ++i)
		for (unsigned j = 0; j < 4; ++j) SoftXOR2Elements(X3[i][j], X2[i][((j >> 1) & 1)], X2[i + 2][(j & 1)]);


	// fill the X4 array
	SoftXOR2Elements(X4[0], X3[0][0], X3[1][0]), SoftXOR2Elements(X4[1], X3[0][3], X3[1][3]);
	SoftXOR2Elements(X4[2], X3[0][3], X3[1][0]), SoftXOR2Elements(X4[3], X3[0][0], X3[1][3]);
	SoftXOR2Elements(X4[4], X3[0][2], X3[1][2]), SoftXOR2Elements(X4[5], X3[0][1], X3[1][1]);
	SoftXOR2Elements(X4[6], X3[0][1], X3[1][2]), SoftXOR2Elements(X4[7], X3[0][2], X3[1][1]);
	/*--- The end of Algorithm 3, Step 2 ---*/

	/*--- Algorithm 3, Step 3. Here we use X4 array to compute R_y*/
	for (unsigned i = 0; i < 8; ++i) {
		pPathScores[i + 8] = pPathScores[i];
		setPairEW(X4[i], pPathScores[i], pPathScores[i + 8]);
	}
	SUM_BLOCK_COUNT(8);
	/*--- The end of Algorithm 3, Step 3 ---*/

	/*--- Algorithm 3, Step 4. We have one maximum from FHT, compute the remainig one ---*/
	unsigned v8 = X4[maxCW[0]] < 0;

	maxCW[v8] = maxCW[0];
	maxEW[v8] = pPathScores[maxCW[v8] + 8 * v8];

	// we store the values of maximal path score for each v8 to use it at the next phase, see Section V-C
	maxEW[v8 ^ 1] = pPathScores[8 * (v8 ^ 1)], maxCW[v8 ^ 1] = 0;
	for (unsigned i = 1; i < 8; ++i) {
		if (pPathScores[i + 8 * (v8 ^ 1)] > maxEW[v8 ^ 1]) {
			maxEW[v8 ^ 1] = pPathScores[i + 8 * (v8 ^ 1)];
			maxCW[v8 ^ 1] = i;
		}
	}
	CMP_BLOCK_COUNT(7);
	/*--- The end of Algorithm 3, Step 4 ---*/

	/*--- Algorithm 3, Step 5 ---*/
	m_CurrentLLR = maxEW[0] - maxEW[1];
	SUM_COUNT;
}

// Phase 6, u_6 = v_6 + v_9, the decoding window = {5, 6, 7}
void KernelProcessorTrofimiuk16_345::ProcessU6()
{
	v[8] = m_CurrentC;
	unsigned startIndex = 8 * v[8];

	/*--- Algorithm 3, Step 2. Here we reuse the intermediate LLRs from the previous phase---*/
	SoftCombine2Elements(X4[0], X3[0][0], X3[1][0], v[8]), SoftCombine2Elements(X4[4], X3[0][2], X3[1][2], v[8]);
	SoftCombine2Elements(X4[1], X3[0][3], X3[1][3], v[8]), SoftCombine2Elements(X4[5], X3[0][1], X3[1][1], v[8]);
	SoftCombine2Elements(X4[2], X3[0][3], X3[1][0], v[8]), SoftCombine2Elements(X4[6], X3[0][1], X3[1][2], v[8]);
	SoftCombine2Elements(X4[3], X3[0][0], X3[1][3], v[8]), SoftCombine2Elements(X4[7], X3[0][2], X3[1][1], v[8]);
	/*--- The end of Algorithm 3, Step 2 ---*/

	/*--- Algorithm 3, Step 3.---*/
	for (unsigned i = 0; i < 8; ++i) {
		unsigned v5 = (i >> 2) & 1, v6 = (i >> 1) & 1;

		pPathScores[8 * (v[8] ^ 1) + i] = pPathScores[8 * v[8] + i];
		unsigned Index0 = i + 8 * v6, Index1 = i + 8 * (v6 ^ 1);

		setPairEW(X4[i], pPathScores[Index0], pPathScores[Index1]);
	}
	SUM_BLOCK_COUNT(8);
	/*--- The end of Algorithm 3, Step 3 ---*/

	/*--- Algorithm 3, Step 4a.--*/
	unsigned v6m = (maxCW[v[8]] >> 1) & 1;
	unsigned signX4 = (X4[maxCW[v[8]]] < 0) ? 1 : 0;
	unsigned signM = v6m ^ signX4;

	maxEW[signM] = maxEW[v[8]];
	maxCW[signM] = maxCW[v[8]];
	maxEW[signM ^ 1] = pPathScores[8 * (signM ^ 1)];
	maxCW[signM ^ 1] = 0;

	for (unsigned i = 1; i < 8; ++i) {
		if (pPathScores[i + 8 * (signM ^ 1)] > maxEW[signM ^ 1]) {
			maxEW[signM ^ 1] = pPathScores[i + 8 * (signM ^ 1)];
			maxCW[signM ^ 1] = i;
		}
	}
	CMP_BLOCK_COUNT(7);
	/*--- The end of Algorithm 3, Step 4a ---*/

	/*--- Algorithm 3, Step 5 ---*/
	m_CurrentLLR = maxEW[0] - maxEW[1];
	SUM_COUNT;
}

// Phase 7, u_7 = v_5 + v_6 + v_{10}, the decoding window = {5, 6, 7}
void KernelProcessorTrofimiuk16_345::ProcessU7()
{
	v[9] = m_CurrentC;  // v6 + v9
	unsigned v69 = m_CurrentC;

	/*--- Algorithm 3, Step 2.---*/
	// We obtained these array indices according to Algorithm 1. Here is a nontivial example of lines12-13 of Algorithm
	// 2 Observe that we compute only X3 and X4
	SoftCombine2Elements(X3[0][0], X2[0][0], X2[2][0], v[8] ^ v[9]), SoftCombine2Elements(X3[1][0], X2[1][0], X2[3][0], v[9]);
	SoftCombine2Elements(X3[0][1], X2[0][1], X2[2][1], v[8] ^ v[9]), SoftCombine2Elements(X3[1][1], X2[1][1], X2[3][1], v[9]);
	SoftCombine2Elements(X3[0][2], X2[0][1], X2[2][1], v[8] ^ v[9] ^ 1),
	SoftCombine2Elements(X3[1][2], X2[1][0], X2[3][0], v[9] ^ 1);
	SoftCombine2Elements(X3[0][3], X2[0][0], X2[2][0], v[8] ^ v[9] ^ 1),
	SoftCombine2Elements(X3[1][3], X2[1][1], X2[3][1], v[9] ^ 1);
	SoftCombine2Elements(X3[0][4], X2[0][1], X2[2][0], v[8] ^ v[9]), SoftCombine2Elements(X3[1][4], X2[1][1], X2[3][0], v[9]);
	SoftCombine2Elements(X3[0][5], X2[0][0], X2[2][1], v[8] ^ v[9]), SoftCombine2Elements(X3[1][5], X2[1][0], X2[3][1], v[9]);
	SoftCombine2Elements(X3[0][6], X2[0][0], X2[2][1], v[8] ^ v[9] ^ 1),
	SoftCombine2Elements(X3[1][6], X2[1][1], X2[3][0], v[9] ^ 1);
	SoftCombine2Elements(X3[0][7], X2[0][1], X2[2][0], v[8] ^ v[9] ^ 1),
	SoftCombine2Elements(X3[1][7], X2[1][0], X2[3][1], v[9] ^ 1);

	for (unsigned i = 0; i < 8; ++i) SoftXOR2Elements(X4[i], X3[0][i], X3[1][i]);
	/*--- The end of Algorithm 3, Step 2.---*/

	/*--- Algorithm 3, Step 3.---*/
	for (unsigned i = 0; i < 8; ++i) {
		unsigned v5 = (i >> 2) & 1, v6 = (i >> 1) & 1, v7 = i & 1;
		// This indexing is convenient for the recursive maximum computation, see Section V-B
		unsigned indexOfWord0 = 8 * (v6 ^ v5) + 4 * v5 + 2 * v6 + v7;
		unsigned indexOfWord1 = 8 * (v6 ^ v5 ^ 1) + 4 * v5 + 2 * v6 + v7;
		pPathScores[8 * (v69 ^ 1) + i] = pPathScores[8 * v69 + i];
		setPairEW(X4[i], pPathScores[indexOfWord0], pPathScores[indexOfWord1]);
	}
	SUM_BLOCK_COUNT(8);
	/*--- The end of Algorithm 3, Step 3.---*/

	/*--- Algorithm 3, Step 4b. Here we use recursive maximum computation, see Section V-B---*/
	LLRType *pPreviousStorage = pPathScores;

	for (int i = 2; i >= 0; --i) {
		LLRType *pCurrentStorage = pMaxStorage[i];
		unsigned range = 1u << (i + 1);
		for (unsigned j = 0; j < range; ++j) {
			pCurrentStorage[j] = std::max(pPreviousStorage[2 * j], pPreviousStorage[2 * j + 1]);
			CMP_COUNT;
		}
		pPreviousStorage = pCurrentStorage;
	}
	/*--- The end of Algorithm 3, Step 4b.---*/

	/*--- Algorithm 3, Step 5 ---*/
	m_CurrentLLR = pMaxStorage[0][0] - pMaxStorage[0][1];
	SUM_COUNT;
}

// Phase 8, u_8 = v_5, the decoding window = {6, 7}
void KernelProcessorTrofimiuk16_345::ProcessU8()
{
	v[10] = m_CurrentC;  // v5 + v6 + v10

	// We use recursive maximums obtained at phase 7
	shift = v[10] << 1;
	m_CurrentLLR = pMaxStorage[1][shift] - pMaxStorage[1][shift + 1];
	SUM_COUNT;
}

// Phase 9, u_9 = v_6, the decoding window = {7}
void KernelProcessorTrofimiuk16_345::ProcessU9()
{
	v[5] = m_CurrentC;

	shift = (shift << 1) + (v[5] << 1);
	m_CurrentLLR = pMaxStorage[2][shift] - pMaxStorage[2][shift + 1];
	SUM_COUNT;
}

// Phase 10, u_{10} = v_7
void KernelProcessorTrofimiuk16_345::ProcessU10()
{
	v[6] = m_CurrentC;
	v[9] ^= v[6], v[10] ^= v[5] ^ v[6];

	shift = (shift << 1) + (v[6] << 1);
	m_CurrentLLR = pPathScores[shift] - pPathScores[shift + 1];
	SUM_COUNT;
}

// Phase 11, just Arikan SC decoding
void KernelProcessorTrofimiuk16_345::ProcessU11()
{
	v[7] = m_CurrentC;
	unsigned pathIndex = v[5] * 4 + v[6] * 2 + v[7];

	// we use CSE LLRs from the previous phases
	SoftCombine2Elements(X4[pathIndex], X3[0][pathIndex], X3[1][pathIndex], v[10]);

	m_CurrentLLR = X4[pathIndex];
}

// Phase 12, just Arikan SC decoding
void KernelProcessorTrofimiuk16_345::ProcessU12()
{
	v[11] = m_CurrentC;

	std::fill_n(pC, 8, 0);
	pC[5] = v[5], pC[6] = v[6], pC[7] = v[7];
	Arikan(3, pC);

	// Again, we USE CSE from prewious phase
	for (unsigned j = 0; j < 8; ++j) pS[1][j] = (pC[j] != 0) ? X1[j + 8] : X1[j];

	std::copy_n(v + 8, 4, pC);
	Arikan(2, pC);

	SoftCombine(pS[2], pS[1], pC, 4);
	SoftXOR(pS[3], pS[2], 2);
	SoftXOR(pS[4], pS[3], 1);

	m_CurrentLLR = pS[4][0];
}

// Phase 13, just Arikan SC decoding
void KernelProcessorTrofimiuk16_345::ProcessU13()
{
	v[12] = m_CurrentC;
	SoftCombine(pS[4], pS[3], v + 12, 1);

	m_CurrentLLR = pS[4][0];
}

// Phase 14, just Arikan SC decoding
void KernelProcessorTrofimiuk16_345::ProcessU14()
{
	v[13] = m_CurrentC;
	pC[0] = v[12] ^ v[13], pC[1] = v[13];
	SoftCombine(pS[3], pS[2], pC, 2);
	SoftXOR(pS[4], pS[3], 1);

	m_CurrentLLR = pS[4][0];
}

// Phase 15, just Arikan SC decoding
void KernelProcessorTrofimiuk16_345::ProcessU15()
{
	v[14] = m_CurrentC;

	SoftCombine(pS[4], pS[3], &m_CurrentC, 1);
	m_CurrentLLR = pS[4][0];
}

tBit *KernelProcessorTrofimiuk16_345::GetCodeword()
{
	v[15] = m_CurrentC;

	std::copy_n(v, 16, pC);
	Arikan(4, pC);

	return pC;
};

// process next required symbol
void KernelProcessorTrofimiuk16_345::ProcessNextSymbol()
{
	unsigned phase = m_LastCalledPhase + 1;

	switch (phase) {
	case 0: ProcessU0(); break;
	case 1: ProcessU1(); break;
	case 2: ProcessU2(); break;
	case 3: ProcessU3(); break;
	case 4: ProcessU4(); break;
	case 5: ProcessU5(); break;
	case 6: ProcessU6(); break;
	case 7: ProcessU7(); break;
	case 8: ProcessU8(); break;
	case 9: ProcessU9(); break;
	case 10: ProcessU10(); break;
	case 11: ProcessU11(); break;
	case 12: ProcessU12(); break;
	case 13: ProcessU13(); break;
	case 14: ProcessU14(); break;
	case 15: ProcessU15(); break;
	default: break;
	}

	m_LastCalledPhase = phase;
};

void KernelProcessorTrofimiuk16_345::CopyState(KernelProcessor &srca)
{
	KernelProcessorTrofimiuk16_345 &src = (KernelProcessorTrofimiuk16_345 &)srca;

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
	std::copy_n(src.v, 16, v);

	maxCW[0] = src.maxCW[0], maxCW[1] = src.maxCW[1];
}
