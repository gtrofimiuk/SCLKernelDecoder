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
#ifndef LIST_KERNEL_DECODER
#define LIST_KERNEL_DECODER

#include "TVEngineKernel.hpp"
#include "KernelProcessorFactory.hpp"

using namespace std;

class CListKernelPolarDecoder: public CTVEngineKernel
{
	unsigned m_Dimension;
	unsigned m_MinDistance;
	unsigned m_NumOfChecks;

	tBit *m_pEncodingBuffer;
	tBit *m_pEncodingBuffer2;

	/// path scores in SCL
	LLRType *m_pR;
	bool *m_pActivePaths;
	tPairOfLLRTypeUnsigned *m_pSortingBuffer;

	/// arrays for synaaic frozen symbols processing
	tDFMask *m_pSymbolFreezingMasks;
	unsigned *m_pDFValueBits;
	unsigned *m_pDecisionPoints;
	unsigned **m_ppFreezingConstraints;

	/// select kernel processor from specification
	void KernelSelector(string kernelName);
public:
	CListKernelPolarDecoder(std::ifstream &SpecFile  // code specification file
	                  ,unsigned MaxNumOfPaths);
	~CListKernelPolarDecoder(){};

	unsigned GetLength() { return m_Length; }
	unsigned GetDimension() { return m_Dimension; }
	string GetKernelName() { return m_pKernelProcessorFactory->GetName(); }

	void Encode(const tBit *pSrc, tBit *pEncoded);
	void fillFreezingConstraintsFromSpec(std::ifstream &SpecFile);

	// soft-input decoder. Produces  a list of information vectors
	// returns the number of vectors obtained, or -1 in case of error
	int Decode(const LLRType *pLLR,  // log(P{c_i=1|y_i}/P{c_i=0|y_i})
	           tBit *pCodeword  // output codeword
	);
};

/** \brief Compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l, 0\leq i<d
    where l is kernel dimension (m_Size), and d is stride
  \param NumOfRows number of rows in a matrix
  \param NumOfColumns number of columns in a matrix
  \param Stride size of an input block in a matrix
  \param pInputVector the input block vector
  \param pMatrix the matrix
  \param pOutputVector the output vector */
void MatrixMultiply(unsigned NumOfRows, unsigned NumOfColumns, unsigned Stride, const tBit *pInputVector, const tBit *pMatrix, tBit *pOutputVector);
#endif
