/*------------------------------------------------------------------------
Some routines for SCL decoding

Copyright 2020 Peter Trifonov

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
#include "TVMemoryEngine.hpp"
#include "Stack.hpp"

void CTVMemoryEngine::InitMemoryEngine(unsigned MaxNumOfPaths)
{
	m_ListSize = MaxNumOfPaths;
	m_pPathPhases = new unsigned[m_ListSize];
	m_pInactivePathIndices = new unsigned[m_ListSize + 1];  // this is a stack
	m_pDFMask = new tDFMask[MaxNumOfPaths];
}

CTVMemoryEngine::CTVMemoryEngine(unsigned MaxNumOfPaths)
: m_ListSize(MaxNumOfPaths)
, m_pPathIndex2ArrayIndex(nullptr)
, m_pInactiveArrayIndices(nullptr)
, m_pArrayReferenceCount(nullptr)
{
	InitMemoryEngine(MaxNumOfPaths);
}

CTVMemoryEngine::~CTVMemoryEngine()
{
	delete[] m_pPathPhases;
	delete[] m_pInactivePathIndices;

	delete[] m_pDFMask;
	delete[] m_pPathIndex2ArrayIndex;
	delete[] m_pArrayReferenceCount;
	delete[] m_pInactiveArrayIndices;
}

void CTVMemoryEngine::Init(unsigned LogLength, unsigned UpperLayerBoundary, unsigned MemorySize)
{

	m_LogLength = LogLength;
	m_UpperLayer = UpperLayerBoundary;

	m_pPathIndex2ArrayIndex = new unsigned[(m_UpperLayer + 1) * m_ListSize];
	m_pInactiveArrayIndices = new unsigned[(m_UpperLayer + 1) * (m_ListSize + 1)];
	m_EmptyIndex = (m_UpperLayer + 1) * m_ListSize;
	m_pArrayReferenceCount = new unsigned[m_EmptyIndex + 1];
	for (unsigned lambda = 0; lambda <= m_UpperLayer; lambda++) {
		for (unsigned s = 0; s < m_ListSize; s++)
			m_pPathIndex2ArrayIndex[s * (m_UpperLayer + 1) + lambda] = m_EmptyIndex;
	}

	Cleanup();
}

void CTVMemoryEngine::Cleanup()
{
	// final initialization
	unsigned *pInactiveArray = m_pInactiveArrayIndices + 0 * (m_ListSize + 1);

	for (unsigned lambda = 0; lambda <= m_UpperLayer; lambda++) {
		pInactiveArray[0] = m_ListSize;   // mark the stack as empty
		pInactiveArray[m_ListSize] = ~0;  // force on-demand initialization of the stack
		pInactiveArray += m_ListSize + 1;
	}

	m_pInactivePathIndices[0] = m_ListSize;   // mark the stack as empty
	m_pInactivePathIndices[m_ListSize] = ~0;  // force on-demand initialization of the stack

	m_pArrayReferenceCount[m_EmptyIndex] = (~0u) / 2u;
}

unsigned CTVMemoryEngine::AssignInitialPath()
{
	unsigned l = Pop(m_pInactivePathIndices);

	unsigned *pIndexArray = GetIndexArrayPointer(l);
	for (unsigned lambda = 0; lambda <= m_UpperLayer; lambda++) pIndexArray[lambda] = m_EmptyIndex;

	m_pDFMask[l] = 0;
	return l;
}

void CTVMemoryEngine::KillPath(unsigned l)
{
	// mark the path index as inactive
	Push(l, m_pInactivePathIndices);

	// disassociate arrays with path index
	const unsigned *pArrayIndex = GetIndexArrayPointer(l);
	unsigned *pInactiveArray = m_pInactiveArrayIndices + 0 * (m_ListSize + 1);

	for (unsigned lambda = 0; lambda <= m_UpperLayer; lambda++) {
		unsigned s = pArrayIndex[lambda];
		assert(m_pArrayReferenceCount[s] > 0);
		if (!--m_pArrayReferenceCount[s]) {
			assert(s != m_EmptyIndex);
			Push(s, pInactiveArray);
		}
		pInactiveArray += (m_ListSize + 1);
	}
}

unsigned CTVMemoryEngine::ClonePath(unsigned l)
{
	unsigned l1 = Pop(m_pInactivePathIndices);
	// make l1 reference same arrays as l

	const unsigned Z = 0;
	const unsigned *pIndexArray = GetIndexArrayPointer(l);
	unsigned *pIndexArray1 = GetIndexArrayPointer(l1);

	for (unsigned lambda = Z; lambda <= m_UpperLayer; lambda++) {
		unsigned s = pIndexArray[lambda];
		pIndexArray1[lambda] = s;
		m_pArrayReferenceCount[s]++;
	}

	m_pDFMask[l1] = m_pDFMask[l];
	return l1;
}
