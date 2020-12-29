/*------------------------------------------------------------------------
Some defines

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
#ifndef MISC_H
#define MISC_H
#include <algorithm>
#include <cassert>
#include <cstdarg>
#include <exception>
#include <csetjmp>

typedef unsigned char tBit;
typedef float LLRType;
typedef unsigned long long tDFMask;
typedef std::pair<LLRType, unsigned> tPairOfLLRTypeUnsigned;

extern std::jmp_buf JUMP_BUFFER;

#define ERROR_MESSAGE_BUFFER_SIZE 1024
extern char ErrorMessageBuffer[ERROR_MESSAGE_BUFFER_SIZE];

class Exception :public std::exception
{
public:
	Exception(const char* Format, ...)
	{
		va_list args;
		va_start(args, Format);
		vsprintf_s(ErrorMessageBuffer, ERROR_MESSAGE_BUFFER_SIZE-1, Format, args);
		va_end(args);
		*(std::exception*)this=std::exception(ErrorMessageBuffer);
	};

	Exception() {
	}
};

//operation counters
extern unsigned long long SumCount;
extern unsigned long long CmpCount;

#ifdef OPERATION_COUNTING
#define SUM_COUNT SumCount++
#define CMP_COUNT CmpCount++
#define SUM_BLOCK_COUNT(N) SumCount+=N
#define CMP_BLOCK_COUNT(N) CmpCount+=N

class CPairComparator : public std::greater<tPairOfLLRTypeUnsigned>
{

public:
	bool operator()(const tPairOfLLRTypeUnsigned &_Left, const tPairOfLLRTypeUnsigned &_Right) const
	{
		CMP_COUNT;
		return std::greater<tPairOfLLRTypeUnsigned>::operator()(_Left, _Right);
	};
};
#else
#define SUM_COUNT 
#define CMP_COUNT 
#define SUM_BLOCK_COUNT(N)
#define CMP_BLOCK_COUNT(N)

typedef std::greater<tPairOfLLRTypeUnsigned> CPairComparator;
#endif

#endif