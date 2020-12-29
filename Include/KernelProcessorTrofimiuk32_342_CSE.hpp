/*------------------------------------------------------------------------
CSE indices of K_32 kernels

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

/// bit-reversal permutations
const unsigned BR3[] = { 0, 4, 2, 6, 1, 5, 3, 7 };
const unsigned BR4[] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };
/*
phase 8, DW = [5, 6, 7], constraints: []
NumOfBits in phase number =  3 phase = 8
['00000000', '11111111', '11110000', '00001111', '11001100', '00110011', '00111100', '11000011'] 8
['0000', '0011', '1100', '1111'] 4
['00', '11'] 2
*/

const unsigned U8_X5[8][2] = {
	{ 0, 0 }, { 3, 3 }, { 3, 0 }, { 0, 3 }, { 2, 2 }, { 1, 1 }, { 1, 2 }, { 2, 1 },
};

/*
phase 10, DW = [5, 6, 7], constraints: ['9=5+6']
NumOfBits in phase number =  4 phase = 8+2
['00000000|00', '11111111|00', '11110000|11', '00001111|11', '11001100|11', '00110011|11', '00111100|00', '11000011|00'] 8
['0000|0', '0011|0', '1100|0', '1111|0', '0000|1', '0011|1', '1100|1', '1111|1'] 8
['00', '11'] 2
*/
const unsigned U10_X5[8][2] = {
	{ 0, 0 }, { 3, 3 }, { 7, 4 }, { 4, 7 }, { 6, 6 }, { 5, 5 }, { 1, 2 }, { 2, 1 },
};

const unsigned U10_X4[8][2] = {
	{ 0, 0 }, { 0, 1 }, { 1, 0 }, { 1, 1 }, { 0, 0 }, { 0, 1 }, { 1, 0 }, { 1, 1 },
};

/*
phase 16, DW = [12, 13, 14, 15], constraints: []
NumOfBits in phase number =  4 phase = 16
['0000000000000000', '1111111111111111', '1111111100000000', '0000000011111111', '1111000011110000', '0000111100001111',
'0000111111110000', '1111000000001111', '1111000000000000', '0000111111111111', '0000111100000000', '1111000011111111',
'0000000011110000', '1111111100001111', '1111111111110000', '0000000000001111'] 16
['00000000', '00001111', '11110000', '11111111'] 4
['0000', '1111'] 2
['00', '11'] 2
*/
const unsigned U16_X5[16][2] = {
	{ 0, 0 }, { 3, 3 }, { 3, 0 }, { 0, 3 }, { 2, 2 }, { 1, 1 }, { 1, 2 }, { 2, 1 },
	{ 2, 0 }, { 1, 3 }, { 1, 0 }, { 2, 3 }, { 0, 2 }, { 3, 1 }, { 3, 2 }, { 0, 1 },
};

/*
phase 18, DW = [14, 15], constraints: []
NumOfBits in phase number =  5 phase = 16+2
['0000000000000000|00', '1111111111111111|00', '1111111100000000|00', '0000000011111111|00'] 4
['00000000|0', '11111111|0'] 2
['0000', '1111'] 2
['00', '11'] 2
['0', '1'] 2
*/
const unsigned U18_X5[4][2] = {
	{ 0, 0 },
	{ 1, 1 },
	{ 1, 0 },
	{ 0, 1 },
};