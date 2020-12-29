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
#pragma once

#include "misc.hpp"

class CTVMemoryEngine
{
protected:

    unsigned m_EmptyIndex;
    unsigned m_UpperLayer;

    /// This contains m+1 stacks of size L
    unsigned* m_pInactiveArrayIndices;
public:
	tDFMask*  m_pDFMask;
    unsigned* m_pArrayReferenceCount;
    /// Stack of size L
    unsigned* m_pInactivePathIndices;
    unsigned* m_pPathIndex2ArrayIndex;
    unsigned  m_ListSize;
    /// The phase of each path
    unsigned* m_pPathPhases;
    unsigned* m_pBlocks;
    unsigned* m_pBlocksMax;
    unsigned  m_LogLength;

    CTVMemoryEngine() {};
    CTVMemoryEngine(unsigned MaxNumOfPaths);
    virtual ~CTVMemoryEngine();

    /// Initialize data structures for given max num of paths
    void InitMemoryEngine(unsigned MaxNumOfPaths);

    /// Initialize the Tal-Vardy  data structures for specific number of layers. 
    /// Only layers up to UpperLayerBoundary will be actually initialized
    void Init(unsigned LogLength, unsigned UpperLayerBoundary, unsigned MemorySize);

    /// cleanup all datastructures
    void Cleanup();

    inline unsigned * GetIndexArrayPointer(unsigned l) const
    {
        return m_pPathIndex2ArrayIndex + l * (m_UpperLayer + 1);
    }

    inline unsigned GetPathIndex(const unsigned *pArrayIndex) const
    {
        return static_cast<unsigned>((pArrayIndex - m_pPathIndex2ArrayIndex) / (m_UpperLayer + 1));
    }

    inline unsigned UpdateIndex(unsigned lambda, unsigned *pIndexArray)
    {
        unsigned Index = pIndexArray[lambda];

        if (m_pArrayReferenceCount[Index] > 1) {
            m_pArrayReferenceCount[Index]--;
            Index = GetNewArrayIndex(lambda);
            m_pArrayReferenceCount[Index] = 1;
            pIndexArray[lambda] = Index;
        }

        return Index;
    }

    /// Get an index of an unused array at layer lambda
    unsigned GetNewArrayIndex(unsigned lambda)
    {
        unsigned* pStack = m_pInactiveArrayIndices + lambda*(m_ListSize + 1);
        assert(pStack[0] != ~0u);
        if (!pStack[0])
            //out of memory
            longjmp(JUMP_BUFFER, 1);
        if (pStack[pStack[0]] == ~0u) {
            //the stack was not properly initialized yet. The memory has not been allocated too
            unsigned s = --pStack[0];
            if (s > 0u)
                pStack[pStack[0]] = ~0u;
            s = s * (m_UpperLayer + 1) + lambda;
            m_pArrayReferenceCount[s] = 1;
            return s;
        } else
            return pStack[pStack[0]--];
    }

    unsigned AssignInitialPath();
    void KillPath(unsigned l);
    unsigned ClonePath(unsigned l);

    unsigned GetMaxNumOfPaths() const
    {
        return m_ListSize;
    }

    unsigned GetNumOfUnusedPaths() const
    {
        return m_pInactivePathIndices[0];
    }
};