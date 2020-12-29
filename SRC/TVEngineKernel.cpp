/*------------------------------------------------------------------------
Extension of SCL decoding for large kernels

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
#include "TVEngineKernel.hpp"

void CTVEngineKernel::Init(unsigned LogLength, unsigned UpperLayerBoundary, unsigned MemoryLimit)
{
	CTVMemoryEngine::Init(LogLength, UpperLayerBoundary, MemoryLimit);

	m_ppProcessorPointer = new KernelProcessor **[m_NumOfLayers * m_ListSize];
	unsigned currentSize = m_Length;
	for (unsigned lambda = 0; lambda <= m_UpperLayer; lambda++) {
		currentSize /= m_KernelSize;
		for (unsigned s = 0; s < m_ListSize; s++) {
			unsigned index = s * m_NumOfLayers + lambda;
			m_ppProcessorPointer[index] = new KernelProcessor *[currentSize];
			for (unsigned j = 0; j < currentSize; ++j) {
				m_ppProcessorPointer[index][j] = m_pKernelProcessorFactory->GetKernelProcessor();
			}
		}
	}
}

KernelProcessor **CTVEngineKernel::GetProcessorPointerToWrite(unsigned lambda, unsigned *pIndexArray)
{
	unsigned prevIndex = pIndexArray[lambda];
	unsigned prevCount = m_pArrayReferenceCount[prevIndex];
	bool flag = (prevCount > 1 && prevCount <= m_ListSize);

	unsigned s = UpdateIndex(lambda, pIndexArray);

	if (flag) {
		unsigned numOfKernels = m_Length / (unsigned)pow(m_KernelSize, lambda + 1);
		KernelProcessor **srcProcs = m_ppProcessorPointer[prevIndex];
		KernelProcessor **dstProcs = m_ppProcessorPointer[s];
		for (unsigned k = 0; k < numOfKernels; ++k) {
			dstProcs[k]->CopyState(*srcProcs[k]);
		}
	}
	return m_ppProcessorPointer[s];
}

void CTVEngineKernel::IterativelyCalcS(unsigned phi, unsigned *pIndexArray)
{
	unsigned depth = CalculateDepth(phi);

	for (unsigned stage = 0; stage < depth + 1; ++stage) {
		unsigned lambda = m_NumOfLayers - (depth + 1) + stage;
		unsigned numOfKernels = m_Length / (unsigned)pow(m_KernelSize, lambda + 1);

		KernelProcessor **pProcessorsToSend = nullptr;
		if (stage > 0) pProcessorsToSend = GetProcessorPointerToRead(lambda - 1, pIndexArray);

		KernelProcessor **pProcessorsToReceive = GetProcessorPointerToWrite(lambda, pIndexArray);

		for (unsigned k = 0; k < numOfKernels; ++k) {
			KernelProcessor &kernelProcessor = *pProcessorsToReceive[k];
			if (stage > 0) {
				LLRType *pKernelLLRPointerToWrite = kernelProcessor.GetInputLLRPointer();

				// get LLRs from the previous layer
				for (unsigned i = 0; i < m_KernelSize; ++i)
					pKernelLLRPointerToWrite[i] = pProcessorsToSend[numOfKernels * i + k]->GetCurrentLLR();
			}
			kernelProcessor.ProcessNextSymbol();
		}
	}
	return;
}

void CTVEngineKernel::IterativelyUpdateC(unsigned phi, unsigned *pIndexArray)
{
	unsigned depth = CalculateDepth(phi);

	for (unsigned stage = 0; stage < depth; ++stage) {
		unsigned lambda = m_NumOfLayers - 1 - stage;
		unsigned numOfKernels = m_Length / (unsigned)pow(m_KernelSize, lambda + 1);

		KernelProcessor **pProcessorsToSend = GetProcessorPointerToRead(lambda, pIndexArray);
		KernelProcessor **pProcessorsToReceive = GetProcessorPointerToWrite(lambda - 1, pIndexArray);

		for (unsigned k = 0; k < numOfKernels; ++k) {
			KernelProcessor &kernelProcessor = *pProcessorsToSend[k];
			tBit *currentC = kernelProcessor.GetCodeword();
			for (unsigned i = 0; i < m_KernelSize; ++i)
				pProcessorsToReceive[numOfKernels * i + k]->LoadNextSymbol(currentC[i]);
		}
	}
	return;
}

void CTVEngineKernel::LoadLLRs(const LLRType *pLLR, unsigned l)
{
	unsigned *pIndexArray = GetIndexArrayPointer(l);
	KernelProcessor **pProcessorsToReceive = GetProcessorPointerToWrite(0, pIndexArray);

	unsigned NumOfUpperLayerKernels = m_Length / m_KernelSize;
	for (unsigned k = 0; k < NumOfUpperLayerKernels; ++k) {
		LLRType *pKernelLLRPointerToWrite = pProcessorsToReceive[k]->GetInputLLRPointer();
		for (unsigned i = 0; i < m_KernelSize; ++i) {
			pKernelLLRPointerToWrite[i] = pLLR[NumOfUpperLayerKernels * i + k];
		}
	}
}

void CTVEngineKernel::WriteCodeword(tBit *pCodeword, unsigned l)
{
	KernelProcessor **pProcessorsToSend = GetProcessorPointerToRead(0, l);
	unsigned NumOfUpperLayerKernels = m_Length / m_KernelSize;
	for (unsigned k = 0; k < NumOfUpperLayerKernels; ++k) {
		KernelProcessor &kernelProcessor = *pProcessorsToSend[k];
		tBit *currentC = kernelProcessor.GetCodeword();
		for (unsigned i = 0; i < m_KernelSize; ++i) pCodeword[NumOfUpperLayerKernels * i + k] = currentC[i];
	}
	return;
}
