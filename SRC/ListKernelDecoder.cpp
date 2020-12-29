/*------------------------------------------------------------------------
SCL decoder of polar codes, supports polar subcodes

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

#include "ListKernelDecoder.hpp"
#include "Utils.hpp"
#include <fstream>
#include <iostream>
#include <set>
#include <string>

#include "KernelProcessorFactoryTrofimiuk16_345.hpp"
#include "KernelProcessorFactoryTrofimiuk32_342.hpp"

CListKernelPolarDecoder::CListKernelPolarDecoder(std::ifstream &SpecFile  // code specification file
                                                 ,unsigned ListSize)
: CTVEngineKernel(ListSize)
{
	unsigned NumOfShortened, NumOfPunctured;
	SpecFile >> m_Length >> m_Dimension >> m_MinDistance >> m_NumOfLayers >> NumOfShortened >> NumOfPunctured;

	//now shortening and 
	if (!SpecFile) throw Exception("Error reading file header");
	if (m_Dimension > m_Length) throw Exception("Code dimension cannot exceed code length");
	if (NumOfShortened + NumOfPunctured) throw Exception("Puncturing and shortening are not supported yet");

	set<string> kernelNames;
	for (unsigned i = 0; i < m_NumOfLayers; ++i) {
		string kernelName;
		SpecFile >> kernelName;
		kernelNames.insert(kernelName);
	}

	if (kernelNames.size() > 1) throw Exception("Mixed kernels are not supported yet");

	KernelSelector(*kernelNames.begin());

	if (m_Length != pow(m_pKernelProcessorFactory->GetKernelSize(), m_NumOfLayers))
		throw Exception("Code length does not match size of the kernel and number of layers");

	CTVEngineKernel::Init(m_NumOfLayers * m_Length / m_KernelSize, m_NumOfLayers - 1, ~0);

	// arrays for path management in list decoder
	m_pR = new LLRType[m_ListSize];
	m_pActivePaths = new bool[m_ListSize];
	m_pSortingBuffer = new tPairOfLLRTypeUnsigned[m_ListSize * 2];

	m_pEncodingBuffer = new tBit[m_Length];
	m_pEncodingBuffer2 = new tBit[m_Length];

	// arrays and function for processing of dynamic frozen symbols
	m_pSymbolFreezingMasks = new tDFMask[m_Length];
	std::fill_n(m_pSymbolFreezingMasks, m_Length, 0);
	m_pDFValueBits = new unsigned[m_Length];
	std::fill_n(m_pDFValueBits, m_Length, ~0);

	m_NumOfChecks = m_Length - m_Dimension;
	m_ppFreezingConstraints = new unsigned *[m_NumOfChecks];
	// load information symbol constraints
	m_pDecisionPoints = new unsigned[m_Length];
	std::fill_n(m_pDecisionPoints, m_Length, ~0);

	fillFreezingConstraintsFromSpec(SpecFile);

	// DFMask initialization
	set<unsigned> AvailableBits;
	for (unsigned i = 0; i < 8 * sizeof(tDFMask); i++) AvailableBits.insert(i);
	for (unsigned i = 0; i < m_Length; i++) {
		if (m_pDFValueBits[i] < 8 * sizeof(tDFMask))
			// this bit is no longer needed
			AvailableBits.insert(m_pDFValueBits[i]);
		// check if this symbol is a starting one in any of the following dynamic freezing constraints
		for (unsigned j = i + 1; j < m_Length; j++) {
			if ((m_pDecisionPoints[j] != ~0) && m_ppFreezingConstraints[m_pDecisionPoints[j]] &&
			    (*(m_ppFreezingConstraints[m_pDecisionPoints[j]]) == i)) {
				if (AvailableBits.size() == 0)
					throw Exception("Too many dynamic freezing constraints are simultaneously active");
				unsigned B = *AvailableBits.begin();
				m_pDFValueBits[j] = B;
				AvailableBits.erase(B);
				unsigned *pFC = m_ppFreezingConstraints[m_pDecisionPoints[j]];
				do {
					m_pSymbolFreezingMasks[*pFC] ^= 1ull << B;
				} while (*(++pFC));
				// cleanup the bit in order to allow its reuse
				m_pSymbolFreezingMasks[j] ^= 1ull << B;
			};
		};
	}
}

// select kernel processor
// naive implementation
void CListKernelPolarDecoder::KernelSelector(string kernelName)
{
	if (kernelName == "Trofimiuk16_345")
		m_pKernelProcessorFactory = new KernelProcessorFactoryTrofimiuk16_345();
	else if (kernelName == "Trofimiuk32_342")
		m_pKernelProcessorFactory = new KernelProcessorFactoryTrofimiuk32_342();
	else
		throw Exception("Wrong kernel name");

	m_KernelSize = m_pKernelProcessorFactory->GetKernelSize();
}

void CListKernelPolarDecoder::fillFreezingConstraintsFromSpec(std::ifstream &SpecFile)
{
	for (unsigned i = 0; i < m_NumOfChecks; i++) {
		unsigned CurWeight;
		SpecFile >> CurWeight;
		if (CurWeight == 1) {
			m_ppFreezingConstraints[i] = 0;
			unsigned S;
			SpecFile >> S;
			if (!SpecFile) throw Exception("Failed to read specification file");
			if (m_pDecisionPoints[S] != ~0) {
				throw Exception("Duplicated dynamic freezing constraint on symbol %d at line %d", S, i);
			}
			m_pDecisionPoints[S] = i;
		} else {
			m_ppFreezingConstraints[i] = new unsigned[CurWeight];
			for (unsigned j = 0; j < CurWeight; j++) {
				unsigned S;
				SpecFile >> S;
				m_ppFreezingConstraints[i][j] = S;
				if (m_ppFreezingConstraints[i][j] >= m_Length) throw Exception("Invalid dynamic freezing constraint");
			};
			sort(m_ppFreezingConstraints[i], m_ppFreezingConstraints[i] + CurWeight);
			unsigned S = m_ppFreezingConstraints[i][CurWeight - 1];
			if (!SpecFile) throw Exception("Failed to read specification file");
			if (m_pDecisionPoints[S] != ~0) {
				throw Exception("Duplicated dynamic freezing constraint on symbol %d at line %d", S, i);
			}
			m_ppFreezingConstraints[i][CurWeight - 1] = 0;
			m_pDecisionPoints[S] = i;
		};
	};
}

void CListKernelPolarDecoder::Encode(const tBit *pSrc, tBit *pEncoded)
{
	for (unsigned i = 0; i < m_Length; i++) {
		if (m_pDecisionPoints[i] == ~0)
			m_pEncodingBuffer[i] = *(pSrc++);
		else {
			const unsigned *pFC = m_ppFreezingConstraints[m_pDecisionPoints[i]];
			tBit C = 0;
			if (pFC) {
				while (*pFC) {
					C ^= m_pEncodingBuffer[*pFC];
					pFC++;
				};
			};
			assert(C < 2);
			m_pEncodingBuffer[i] = C;
		};
	};

	const tBit *pKernel = m_pKernelProcessorFactory->GetKernel();

	unsigned Stride = 1;
	for (int L = m_NumOfLayers - 1; L >= 0; L--) {
		unsigned NextStride = Stride * m_KernelSize;
		unsigned NumOfBlocks = m_Length / NextStride;
		for (unsigned i = 0; i < NumOfBlocks; i++) {
			tBit *pSrc = m_pEncodingBuffer + i * NextStride;
			tBit *pDest = m_pEncodingBuffer2 + i * NextStride;
			MatrixMultiply(m_KernelSize, m_KernelSize, Stride, pSrc, pKernel, pDest);
		}
		swap(m_pEncodingBuffer, m_pEncodingBuffer2);
		Stride = NextStride;
	}
	std::copy_n(m_pEncodingBuffer, m_Length, pEncoded);
}

int CListKernelPolarDecoder::Decode(const LLRType *pLLR,  // log(P{c_i=0|y_i}/P{c_i=1|y_i})
                                    tBit *pCodeword       // output codeword
)
{
	CTVEngineKernel::Cleanup();

	unsigned initialPath = CTVMemoryEngine::AssignInitialPath();
	CTVEngineKernel::LoadLLRs(pLLR, initialPath);

	std::fill_n(m_pR, m_ListSize, 0);
	std::fill_n(m_pActivePaths, m_ListSize, 0);
	m_pActivePaths[initialPath] = true;

	for (unsigned phi = 0; phi < m_Length; ++phi) {
		unsigned l = 0;
		bool isSymbolUnfrozen = m_pDecisionPoints[phi] == ~0;

		for (unsigned i = 0; i < m_ListSize; ++i) {
			if (m_pActivePaths[i]) {
				unsigned *pCurrentIndexArray = GetIndexArrayPointer(i);

				if (phi > 0) IterativelyUpdateC(phi, pCurrentIndexArray);
				IterativelyCalcS(phi, pCurrentIndexArray);

				KernelProcessor **proc = GetProcessorPointerToRead(m_NumOfLayers - 1, pCurrentIndexArray);
				KernelProcessor *bottomProcessor = proc[0];
				LLRType currentLLR = bottomProcessor->GetCurrentLLR();

				if (!isSymbolUnfrozen) {

					tBit CurValue;
					if (m_ppFreezingConstraints[m_pDecisionPoints[phi]]) {
						// update dynamic frozen symbol
						CurValue = (m_pDFMask[i] >> m_pDFValueBits[phi]) & 1;
						if (CurValue) m_pDFMask[i] ^= m_pSymbolFreezingMasks[phi];
					} else
						CurValue = 0;

					if (CurValue) {
						if (currentLLR > 0) m_pR[i] -= currentLLR;
					} else {
						if (currentLLR < 0) m_pR[i] += currentLLR;
					};

					if (m_ListSize > 1) {
						SUM_COUNT;
						CMP_COUNT;
					}
					
					bottomProcessor->LoadNextSymbol(CurValue);
				} else {
					bool X = currentLLR < 0;

					m_pSortingBuffer[l].first = m_pR[i];
					m_pSortingBuffer[l++].second = 2 * i + X;

					m_pSortingBuffer[l].first = m_pR[i] - fabs(currentLLR);
					m_pSortingBuffer[l++].second = 2 * i + !X;

					if (m_ListSize > 1) SUM_COUNT;
				}
			}
		}

		if (isSymbolUnfrozen) {
			if (l > m_ListSize) {
				if (m_ListSize > 1) sort(m_pSortingBuffer, m_pSortingBuffer + l, CPairComparator());
				l = m_ListSize;
			}
			// compute the number of clones for each path
			unsigned *pNumOfClones = new unsigned[m_ListSize];
			std::fill_n(pNumOfClones, m_ListSize, 0);

			for (unsigned i = 0; i < l; i++) {
				pNumOfClones[m_pSortingBuffer[i].second >> 1]++;
			}

			for (unsigned i = 0; i < m_ListSize; i++)
				if (m_pActivePaths[i] && !pNumOfClones[i]) {
					KillPath(i);
					m_pActivePaths[i] = false;
				}

			for (unsigned i = 0; i < m_ListSize; i++) {
				if (pNumOfClones[i]) {
					KernelProcessor **proc = GetProcessorPointerToRead(m_NumOfLayers - 1, i);
					KernelProcessor *bottomProcessor = proc[0];
					LLRType S = bottomProcessor->GetCurrentLLR();
					tBit X = (S < 0);
					bottomProcessor->LoadNextSymbol(X);
					if (pNumOfClones[i] == 1) {
						if (X) m_pDFMask[i] ^= m_pSymbolFreezingMasks[phi];
					} else {
						unsigned PID2 = ClonePath(i);
						m_pActivePaths[PID2] = true;
						KernelProcessor **clonedProcessors = GetProcessorPointerToWrite(m_NumOfLayers - 1, PID2);
						KernelProcessor *clonedProc = clonedProcessors[0];
						clonedProc->LoadNextSymbol(X ^ 1);
						m_pR[PID2] = m_pR[i] - fabs(S);

						if(m_ListSize > 1) SUM_COUNT;

						m_pDFMask[PID2] = m_pDFMask[i];
						if (X)
							m_pDFMask[i] ^= m_pSymbolFreezingMasks[phi];
						else
							m_pDFMask[PID2] ^= m_pSymbolFreezingMasks[phi];
					}
				}
			}
		}
	}

	// find the best path and write the codeword
	LLRType BestR = -HUGE_VAL;
	unsigned BestPID = ~0;
	for (unsigned i = 0; i < m_ListSize; i++) {
		if (m_pActivePaths[i]) {
			if (m_ListSize > 1) CMP_COUNT;

			if (m_pR[i] > BestR) {
				BestR = m_pR[i];
				BestPID = i;
			}
		}
	}

	// extract codeword and information bits
	IterativelyUpdateC(m_Length, BestPID);
	WriteCodeword(pCodeword, BestPID);

	return 1;
}

/// compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l, 0\leq i<d
/// where l is kernel dimension (m_Size), and d is stride
void MatrixMultiply(unsigned NumOfRows, unsigned NumOfColumns, unsigned Stride, const tBit *pInputVector, const tBit *pMatrix, tBit *pOutputVector)
{
	std::fill_n(pOutputVector, Stride * NumOfColumns, 0);
	for (unsigned j = 0; j < NumOfRows; j++) {
		const tBit *pS = pInputVector + j * Stride;
		for (unsigned i = 0; i < NumOfColumns; i++) {
			if (pMatrix[j * NumOfColumns + i]) {
				tBit *pD = pOutputVector + i * Stride;
				// XOR
				for (unsigned k = 0; k < Stride; ++k) pD[k] ^= pS[k];
			}
		}
	}
}
