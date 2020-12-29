/*------------------------------------------------------------------------
Some defines

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
#include <cassert>
/*A stack  of capacity L is implemented as an array S of size L+1.
S[0] contains the number of elements in the stack
*/
inline void StackInit(const unsigned &capacity, unsigned *pStack)
{
	pStack[0] = capacity;
	pStack[capacity] = ~0u;
}

inline void Push(const unsigned &x, unsigned *pStack) { pStack[++pStack[0]] = x; };

inline unsigned Pop(unsigned *pStack)
{
	assert(pStack[0] > 0);
	if (pStack[pStack[0]] == ~0u) {
		// the stack was not properly initialized yet
		unsigned s = --pStack[0];
		if (s > 0) pStack[pStack[0]] = ~0u;
		return s;
	} else
		return pStack[pStack[0]--];
};

inline unsigned Peek(const unsigned *pStack) { return pStack[pStack[0]] != ~0u ? pStack[pStack[0]] : pStack[0] - 1; }

inline unsigned IsInit(const unsigned *pStack) { return pStack[pStack[0]] != ~0u; }

inline bool IsEmpty(const unsigned *pStack) { return pStack[0] == 0u; }
