/*------------------------------------------------------------------------
Simulation, entry point of application

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
#include <conio.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <random>
#include <string>

#include "ListKernelDecoder.hpp"

using namespace std;

void printUsage();

int main(int argc, char *argv[])
{
	unsigned NumOfArguments = 6;
	if (argc != NumOfArguments) {
		printUsage();
		return 1;
	};
	unsigned Arg = 1;
	const char *pSpecFile = argv[Arg++];         // 1
	unsigned ListSize = atoi(argv[Arg++]);       // 2
	double SNR = atof(argv[Arg++]);              // 3
	unsigned MaxErrors = atoi(argv[Arg++]);      // 4
	unsigned MaxIterations = atoi(argv[Arg++]);  // 5

	unsigned MemoryLimit = 0;

	try {
		ifstream SpecFile(pSpecFile, ios::in);

		CListKernelPolarDecoder Codec(SpecFile, ListSize);

		unsigned N = Codec.GetLength();
		unsigned K = Codec.GetDimension();

		cout << '(' << Codec.GetLength() << ',' << Codec.GetDimension()
		     << ") code with kernel "
		     << Codec.GetKernelName()
		     << " initialized\n";

		unsigned NumOfErrors = 0;
		unsigned NumOfMLError = 0;

		double NoiseStddev = sqrt(0.5 * pow(10, -SNR / 10) * Codec.GetLength() / Codec.GetDimension());
		double NoiseVariance = NoiseStddev * NoiseStddev;

		unsigned seed = 8;
		std::mt19937_64 generator(seed);
		std::uniform_int_distribution<unsigned> bitDistribution(0, 1);
		std::normal_distribution<LLRType> noiseDistribution(0, NoiseStddev);

		tBit *pData = new tBit[K];
		tBit *pEncoded = new tBit[N];
		tBit *pDecoded = new tBit[N];
		LLRType *pTransmitted = new LLRType[N];

		unsigned It = 1;
		for (It; It <= MaxIterations; ++It) {

			// generate random data
			for (unsigned i = 0; i < K; ++i) pData[i] = bitDistribution(generator);

			// encoding
			Codec.Encode(pData, pEncoded);

			// modulation
			for (unsigned i = 0; i < N; ++i) pTransmitted[i] = 2 * pEncoded[i] - 1;

			// add noise
			for (unsigned i = 0; i < N; ++i) {
				LLRType noise = noiseDistribution(generator);
				pTransmitted[i] += noise;
			}

			// demodulation, computing log(P{c_i=0|y_i}/P{c_i=1|y_i})
			for (unsigned i = 0; i < N; ++i) pTransmitted[i] = -2 * pTransmitted[i] / NoiseVariance;

			// decoding
			Codec.Decode(pTransmitted, pDecoded);

			bool isErrorOccured = false;
			// check if error occured
			for (unsigned i = 0; i < N; ++i) {
				if ((pDecoded[i] != 0) != (pEncoded[i] != 0)) {
					isErrorOccured = true;
					break;
				}
			}

			bool isMLErrorOccured = false;
			// check if ML error occured
			if (isErrorOccured) {
				LLRType CorrelationEncoded = 0;
				LLRType CorrelationDecoded = 0;

				for (unsigned i = 0; i < N; ++i)
					CorrelationEncoded += (pEncoded[i]) ? -pTransmitted[i] : pTransmitted[i];
				for (unsigned i = 0; i < N; ++i)
					CorrelationDecoded += (pDecoded[i]) ? -pTransmitted[i] : pTransmitted[i];

				if (CorrelationDecoded >= CorrelationEncoded) isMLErrorOccured = true;
			}

			NumOfErrors += isErrorOccured;
			NumOfMLError += isMLErrorOccured;

			if (NumOfErrors >= MaxErrors) break;

			if (It % 1000 == 0 && It > 0) {
				cout << It << "\t" << double(NumOfErrors) / double(It) << endl;
			}
		}

		cout << "SNR\tStdDev\tFER\tMLFER" << endl;
		cerr.precision(5);
		cerr << scientific << SNR << '\t' << NoiseStddev << '\t' << double(NumOfErrors) / double(It) << '\t'
		     << double(NumOfMLError) / double(It)
#ifdef OPERATION_COUNTING
		     << '\t' << double(SumCount) / double(It) << '\t' << double(CmpCount) / double(It)
#endif
		     << endl;
	} catch (const exception &ex) {
		cerr << ex.what() << endl;
		return 5;
	};
};

void printUsage()
{
	std::cerr << "cmd_args:\n"
	          << "\tCodeSpecFile ListSize \n"
	          << "SNR MaxErrors MaxIterations\n"
	          << std::endl;
}
