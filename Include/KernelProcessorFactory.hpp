/*------------------------------------------------------------------------
Interface for factory of kernel processors

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

class KernelProcessorFactory
{
public:
	KernelProcessorFactory(const unsigned KernelSize, const tBit *pKernel)
	: m_KernelSize(KernelSize)
	, m_pKernel(pKernel){};
	virtual KernelProcessor *GetKernelProcessor() = 0;
	unsigned GetKernelSize() { return m_KernelSize; }
	const tBit *GetKernel() { return m_pKernel; }
	virtual const char *GetName() const = 0;

protected:
	const unsigned m_KernelSize;
	const tBit *m_pKernel;
};