/*------------------------------------------------------------------------
The interface class of kernel processor for SCL decoder

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

#include "misc.hpp"

class KernelProcessor
{
public:
	//Methods needed to be implemented in kernel processor 

	/// return u*K, K is a kernel
	virtual tBit *GetCodeword() = 0;

	/// process the next input symbol
	virtual void ProcessNextSymbol() = 0;

	/// copy state information to new path in list decoder
	virtual void CopyState(KernelProcessor &dst) = 0;

	//Methods for SCL decoding 

	/// initialization of the Kernel with given LLRs
	void SetInputLLR(LLRType *pLLR)
	{
		m_LastCalledPhase = -1;
		std::copy_n(pLLR, m_KernelSize, m_pLLR);
	}

	/// initialization of the kernel processor, poiner to kernel LLRs is returned
	LLRType *GetInputLLRPointer()
	{
		m_LastCalledPhase = -1;
		return m_pLLR;
	}

	/// return processed LLR
	LLRType GetCurrentLLR() { return m_CurrentLLR; }

	/// send input polarizing transform bit to the kernel
	void LoadNextSymbol(tBit u) { m_CurrentC = u; }

	KernelProcessor(){};

protected:
	unsigned m_KernelSize;

	// input LLRs
	LLRType *m_pLLR;

	LLRType m_CurrentLLR;
	tBit m_CurrentC;

	unsigned m_LastCalledPhase;
};