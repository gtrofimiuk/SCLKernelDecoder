/*------------------------------------------------------------------------
Extension of SCL decoding for large kernels

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
#ifndef TV_ENGINE_KERNEL
#define TV_ENGINE_KERNEL

#include "KernelProcessorFactory.hpp"
#include "TVMemoryEngine.hpp"

#include <iostream>

using namespace std;

class CTVEngineKernel : public CTVMemoryEngine
{
protected:
	/// Length of code
	unsigned m_Length;
	unsigned m_KernelSize;
	/// Length = m_KernelSize^{m_NumOfLayers}
	unsigned m_NumOfLayers;

	KernelProcessorFactory *m_pKernelProcessorFactory;
	KernelProcessor ***m_ppProcessorPointer;

public:
	CTVEngineKernel(unsigned ListSize)
	: CTVMemoryEngine(ListSize){};
	virtual ~CTVEngineKernel(){};

	void Init(unsigned LogLength, unsigned UpperLayerBoundary, unsigned MemoryLimit);

	/// calculate depth of SC decoding at each phase
	unsigned CalculateDepth(unsigned phi)
	{
		unsigned depth = 0;
		unsigned p = phi;
		for (unsigned i = 0; i < m_NumOfLayers - 1; ++i) {
			if (p % m_KernelSize == 0) {
				depth++;
				p /= m_KernelSize;
			} else
				break;
		}
		return depth;
	}

	KernelProcessor **GetProcessorPointerToWrite(unsigned lambda, unsigned *pIndexArray);

	KernelProcessor **GetProcessorPointerToWrite(unsigned lambda, unsigned l)
	{
		unsigned *pIndexArray = GetIndexArrayPointer(l);
		return GetProcessorPointerToWrite(lambda, pIndexArray);
	}

	KernelProcessor **GetProcessorPointerToRead(unsigned lambda, unsigned *pIndexArray)
	{
		unsigned s = pIndexArray[lambda];
		return m_ppProcessorPointer[s];
	}

	KernelProcessor **GetProcessorPointerToRead(unsigned lambda, unsigned l)
	{
		unsigned *pIndexArray = GetIndexArrayPointer(l);
		return GetProcessorPointerToRead(lambda, pIndexArray);
	}

	///recursive update LLRs
	void IterativelyCalcS(unsigned phi, unsigned *pIndexArray);

	/// compute partial sums
	void IterativelyUpdateC(unsigned phi, unsigned *pIndexArray);

	void IterativelyUpdateC(unsigned phi, unsigned l)
	{
		unsigned *pIndexArray = GetIndexArrayPointer(l);
		IterativelyUpdateC(phi, pIndexArray);
		return;
	}

	/// Load input LLRs to kernel processors
	void LoadLLRs(const LLRType *pLLR, unsigned l);

	void Cleanup() { CTVMemoryEngine::Cleanup(); }

	void WriteCodeword(tBit *pCodeword, unsigned l);
};

#endif
