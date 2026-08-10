/*By ANDREW FARRELL; updated by SJG Jul2025
 * Overlap.cpp
 * --------------------------------------------------
 * Assembles k-mers containing variation into contigs
 * that represent the variant sequence
 * --------------------------------------------------
 */

#include <algorithm>
#include <bitset>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>
#include <omp.h>

#include "Util.h"

using namespace std;

bool FullOut = false;

struct OverlapArgs {
	string FastqIn;
	float MinPercent;
    long int MinOverlap;
    long int MinCoverage;
	string NameStub;
    long int hashLength;
    long int ACT;  // the alignment count threshold (how many times a kmer must be found in the sequences to be considered for alignment)
	string OverlapStub;
    long int TrimLCcuttoff;
    long int Threads;
	bool verbose = false;
};

/*
	Completely clears and builds the Hashes hash table, which has (numeric hash of) kmer as keys and a vector as the value, which contains
	all indices of full-length reads in "sequences" array that contain that kmer. This needs to be re-populated occassionally to 
	ensure that the hash table is up to date with the latest sequences after some have been collapsed and moved.
*/
int RebuildHashTable(vector<string>& sequences, int Ai, int hashLength, unordered_map<unsigned long, vector<int>>& Hashes, int Threads)
{
	cout << "\nDestroying HashTable\n";
	Hashes.clear();
	cout << "HashTable destroyed\n";
	cout << "Rebuilding HashTable - starting at " << Ai << endl;
	int size = sequences.size();

	#pragma omp parallel for num_threads(Threads) shared(Hashes, sequences)
	for (int i = Ai; i < size; i++) {

		if (i % 10000 > 1 && i % 10000 < Threads) {
			#pragma omp critical (progressOut) 
			{
				cout << "Hashed " << i << " of " << sequences.size() << "\r";
			}
		}
		// Iterate through the sequence and get hashLength sized chunks
		// If the chunk does NOT have an 'N' in it, add it to the hash table
		string Sequence = sequences[i];
		int LoopLimit = Sequence.size() - hashLength;
		for (int j = 0; j < LoopLimit; j++) {
			string hash = Sequence.substr(j, hashLength);
			size_t found = hash.find('N');
			
			if (found == std::string::npos) {	// npos is a constant for "not found"
				unsigned long LongHash = Util::HashToLong(hash);
				unsigned long RevHash = Util::HashToLong(Util::RevComp(hash));

				// Add sequene index to hash table
				#pragma omp critical(updateHash)
				{
					Hashes[LongHash].push_back(i);
					Hashes[RevHash].push_back(i);
				}
			}
		}
	}
	cout << "\nDone Rebulding HashTable size is " << Hashes.size() << endl;
	return 0;
}

/*
	For a single sequence, iterates through each kmer of hashLength N and looks to see what sequences it is 
	contained within (i.e. the indexes of those sequences in "sequences" array).
	Counts how many times each kmer is found in the sequences, and if it is found more than ACT times,
	it is considered a candidate for alignment.
	Stores the results in the array map, which has the index of the sequence as key and a vector of indexes of sequences
	that contain that kmer as value.

	Note: This function is called within a parallel section from main, hence the critical sections.

	Parameters:
	A: sequence to search
	Ai: the index of the sequence in the original list which is not passed here
	hashLength: the window for creating kmers
	ACT: the alignment count threshold (how many times a kmer must be found in the sequences to be considered for alignment)
	Hashes: the hash table containing the kmer as key and a vector of indexes of sequences containing that kmer as value
*/
int PrepareSearchList(string A, int Ai,	unordered_map<unsigned long, vector<int>>& Hashes, int hashLength, int ACT, map<int, vector<int>>& array,bool& hitPosLimit, bool& hitIndexLimit, int& NumberPos, int& NumberIndex) 
{
	int Alength = A.size();
	map<int, int> Positions;
	int added = 0;

	// Iterate through each kmer in sequence A
	for (int i = 0; i < Alength - hashLength; i++) {
		string hash = A.substr(i, hashLength);
		size_t found = hash.find('N');

		if (found == std::string::npos) {
			// Pull out list of sequences that contain this kmer from hash table
			unsigned long LongHash = Util::HashToLong(hash);
			#pragma omp critical(updateHash) 
			{
				int numMatches = Hashes[LongHash].size();
			
				// Iterate through sequences that contain this kmer
				for (vector<int>::size_type j = 0; j < numMatches; j++) {
					int holder = Hashes[LongHash][j];	// Index of sequence which "holds" this kmer
									
					// Avoid redundant comparisons by only looking at sequences that are after the current sequence A
					// I.e. we've already processed any sequences before us in the 'sequences' array
					if (holder > Ai){
						if (Positions.count(holder) > 0) {
							Positions[holder]++;
							added++;
						} else {
							Positions[holder] = 1;
							added++;
						}
					}

					if (added > 100000) {
						hitPosLimit = true;
						NumberPos = added;
						break;
					}
				}
			}
		}
	}

	if (FullOut) {
		#pragma omp critical(progressOut) 
		{ 
			cout << "done Hashing read" << endl;
		}
	}

	// Sort our postions by the number of times a kmer was found in the sequence at that position
	// i.e. the highest value of the Positions map to the lowest value
	NumberPos = added;
	map<int, int>::iterator uspos;
	multimap<int, int> SortedPositions;

	for (uspos = Positions.begin(); uspos != Positions.end(); ++uspos) {
			if (uspos->second > ACT)
			{
				SortedPositions.insert(std::make_pair(uspos->second,uspos->first));
			}
	}


	if (FullOut) {
		#pragma omp critical(progressOut) 
		{ 
			cout << "found - " << Positions.size() << " possible locations" << endl;
		}
	}

	map<double, int>::iterator pos;
	vector<int> indexes;
	int sanity = 0;

	// Check to make sure that we are meeting our minimum alignment count threshold (ACT)
	for (auto pos = SortedPositions.rbegin() ; pos !=  SortedPositions.rend(); pos++)
    {
		if (pos->first >= ACT) {
            indexes.push_back(pos->second + 0);
			sanity++;

			if (sanity > 1000) {
					hitIndexLimit = true;
					NumberIndex = sanity;
					break;
			}
        }

    }

	NumberIndex = sanity;
	#pragma omp critical (array) 
	{ 
		array[Ai] = indexes; 
	}

	if (FullOut) {
		#pragma omp critical(progressOut) 
		{ 
			cout << "			 " << indexes.size() << " locations passed filter" << endl;
		}
	}
	return 1;
}

int Align3(vector<string>& sequenes, string Ap, string Aq, int Ai, int& overlap, int& BestIndex, float minPercent, bool& PerfectMatch, int MinOverlap,vector<int>& indexes, int Threads, int NumReads) 
{
	int QualityOffset = 33;
	bool verbose = false;
	int Alength = Ap.size();
	int bestScore = 0;

	#pragma omp parallel for num_threads(Threads) shared(BestIndex, overlap, bestScore, PerfectMatch, sequenes)
	for (int i = 0; i < indexes.size(); i++) 
	{
		string A = Ap; 
		int AlengthL = A.size(); 
		int j = indexes[i]; 
		
		string B = sequenes[j];
		float score = 0;
		int Blength = B.size();

		int window = -1;
		int longest = -1;
		bool Asmaller = true;
		bool LocalPerfectMatch = false;

		if (Blength > AlengthL) 
		{
			window = AlengthL;
			longest = Blength;
			Asmaller = false;
		} else 
		{
			Asmaller = true;
			window = Blength;
			longest = AlengthL;
		}

		int MM = window - (window * minPercent);
		int LbestScore = -1;
		int LBestIndex = -1;
		int Loverlap = 0;
		int Acount = 0;
		int Bcount = 0;
		int k;

		for (int i = 0; i <= longest - window; i++) 
		{
			score = 0;

			for (k = 0; k < window; k++)	// base compare loop
			{
				if (A.c_str()[k + Acount] == B.c_str()[k + Bcount]) 
				{
					if (A.c_str()[k + Acount] != 'N') 
					{
						score++;
					}
				}

				if ((k - score) > MM) 
				{
					score = -1;
					break;
				}
			}

			if (Asmaller) {
				Acount++;
			} else {
				Bcount++;
			}

			if (verbose) {cout << "	Score = " << score << endl;}

			float percent = score / (window);
			if (percent >= minPercent) 
			{
				if (LbestScore < score) 
				{
						LbestScore = score;
						LBestIndex = j;
					
					if (Asmaller) {
						Loverlap = i * -1;
					} else {
						Loverlap = i;
					}

					if (score == window) 
					{
						LocalPerfectMatch = true;
						break;
					}
				}
			}
		}

		// If we haven't found a perfect match, continue searching for overlaps
		if (LocalPerfectMatch == false) 
		{
			for (int i = window - 1; i >= MinOverlap; i--) 
			{
				if (verbose) {cout << "i = " << i << endl;}
				score = 0;
				for (k = 0; k <= i; k++) 
				{
					if (verbose) {cout << "	k = " << k << " so A = " << AlengthL - i + k<< " /\\ B = " << 0 + k << endl;}
					if (verbose) {cout << "		 A >> " << A.c_str()[AlengthL - i + k - 1] << "="<< B.c_str()[0 + k] << " << B" << endl;}

					if (A.c_str()[AlengthL - i + k - 1] == B.c_str()[0 + k]) 
					{
						if (B.c_str()[0 + k] != 'N') 
						{
							score++;
						}
					}
					if ((k - score) > MM) {
						score = -1;
						break;
					}
				}
				if (verbose) {cout << "	Score = " << score << endl;}
				float percent = score / (k);
				if (percent >= minPercent) 
				{
					if (LbestScore < score) 
					{
							LbestScore = score;
							LBestIndex = j;
							Loverlap = i - AlengthL + 1;
						if (score == i) 
						{
							break;
						}
					}
				}
			}

			for (int i = window - 1; i >= MinOverlap; i--) 
			{
				if (verbose) {cout << "i = " << i << endl;}

				score = 0;
				for (k = 0; k <= i; k++) 
				{
					if (B.c_str()[Blength - i + k - 1] == A.c_str()[0 + k]) 
					{
						if (A.c_str()[0 + k] != 'N') 
						{
							score++;
						}
					}

					if ((k - score) > MM) 
					{
						score = -1;
						break;
					}
				}

				if (verbose) {cout << "	Score = " << score << endl;}

				float percent = score / (k);

				if (percent >= minPercent) 
				{
					if (LbestScore < score) 
					{
							LbestScore = score;
							LBestIndex = j;
							Loverlap = Blength - i - 1;
						if (score == i) {
							break;
						}
					}
				}
			}
		}
		#pragma omp critical (best)
		{
			if (LbestScore > bestScore ||
        		(LbestScore == bestScore && LBestIndex < BestIndex) ||
        		(LbestScore == bestScore && LBestIndex == BestIndex && Loverlap < overlap)) {
					bestScore = LbestScore;
					BestIndex = LBestIndex;
					overlap = Loverlap;
			}
			// Only want to update this logic if we have found a perfect match
			if (LocalPerfectMatch) {
				PerfectMatch = true;
			}
		}
	}
	return bestScore;
}

/* Combines sequences A and B into a single contiguous string. 
 * k is the offset between A and B, where positive k means A is upstream of B.
 * Aq, Bq, Ad, Bd, As, Bs are the quality strings, depth strings, and strand strings for A and B respectively.
 * Returns the combined string.
 * 
 * Merges sequences based on the following logic:
 * 1. If both sequences have the same base at a position → use that base, take the higher quality score, 
 * and sum the depths (capped at 250)
 * 2. If only one sequence has a base at that position → use that sequence's data
 * 3. If sequences disagree → prefer the base with higher depth, or if depths are equal, prefer the one with higher quality
 * 
 * Updates the reference parameters (Bq, Bd, Bs) with the merged quality scores, depth data, and combined sequence information
 */
string ColapsContigs(string A, string B, int k, string Aq, string& Bq,string Ad, string& Bd, string As, string& Bs) {
	bool verbose = false;
	if (verbose) {cout << "Combining; \n" << A << endl << B << endl;}

	int Asize = A.size();
	int Bsize = B.size();
	int Aoffset = 0;
	int Boffset = 0;
	int window;
	string newString = "";
	string newQual = "";
	string newDepth = "";

	if (k > 0) {
		Aoffset = k;
	} else {
		Boffset = abs(k);
	}

	if (verbose) {cout << "K = " << k << " so Aofset = " << Aoffset<< " and Boffset = " << Boffset << endl;}

	for (int i = 0; i < Asize + Bsize; i++) 
	{
		char Abase = 'Z';
		char Bbase = 'Z';
		char Aqualc = '!';
		char Bqualc = '!';
		unsigned char Adepc = 0;
		unsigned char Bdepc = 0;

		if (((i - Aoffset) >= 0) && ((i - Aoffset) < A.size())) 
		{
			Abase = A.c_str()[i - Aoffset];
			Aqualc = Aq.c_str()[i - Aoffset];
			Adepc = Ad.c_str()[i - Aoffset];
		} else {
			Abase = 'Z';
			Aqualc = '!';
			Adepc = 0;
		}

		if (i - Boffset >= 0 && i - Boffset < B.size()) {
			Bbase = B.c_str()[i - Boffset];
			Bqualc = Bq.c_str()[i - Boffset];
			Bdepc = Bd.c_str()[i - Boffset];
		} else {
			Bbase = 'Z';
			Bqualc = '!';
			Bdepc = 0;
		}

		if (verbose) {cout << "I = " << i << " Bi = " << i - Boffset << " Ai = " << i - Aoffset<< " thus " << Abase << "-" << Bbase << endl;}

		if (Abase == Bbase && Abase != 'Z') {

			newString += Abase;
			if (Aqualc >= Bqualc) {
				newQual += Aqualc;
			} else {
				newQual += Bqualc;
			}
			if ((int)Adepc + (int)Bdepc < 250) {
				newDepth += (Adepc + Bdepc);
			} else {
				newDepth += (char)250;
			}
		} else if (Abase == 'Z' && Bbase != 'Z') {
			newString += Bbase;
			newQual += Bqualc;
			newDepth += Bdepc;
		} else if (Abase != 'Z' && Bbase == 'Z') {
			newString += Abase;
			newQual += Aqualc;
			newDepth += Adepc;

		} else if (Abase != 'Z' && Bbase != 'Z') {
			 if (Adepc > Bdepc) {
				newString += Abase;
				newQual += Aqualc;
				newDepth += Adepc;
			} else if (Adepc < Bdepc) {
				newString += Bbase;
				newQual += Bqualc;
				newDepth += Bdepc;
			} else if (Aqualc >= Bqualc) {
				newString += Abase;
				newQual += Aqualc;
				newDepth += Adepc;
			} else {
				newString += Bbase;
				newQual += Bqualc;
				newDepth += Bdepc;
			}

		} else if (Abase == 'Z' && Bbase == 'Z') {
			Bq = newQual;
			Bd = newDepth;
			break;
		}	
	}
	Bs += As;
	Bq = newQual;
	Bd = newDepth;
	return newString;
}

string TrimNends(string S, string& qual) {
	bool base = false;
	string NewS = "";
	string NewQ = "";
	for (int i = S.size() - 1; i >= 0; i--) {

		if (base) {
			NewS = S.c_str()[i] + NewS;
			NewQ = qual.c_str()[i] + NewQ;
		} else if (S.c_str()[i] != 'A' && S.c_str()[i] != 'C' &&
							 S.c_str()[i] != 'G' && S.c_str()[i] != 'T') {
		} else {
			base = true;
			NewS = S.c_str()[i] + NewS;
			NewQ = qual.c_str()[i] + NewQ;
		}
	}

	S = NewS;
	qual = NewQ;
	base = false;
	NewS = "";
	NewQ = "";

	for (int i = 0; i < S.size(); i++) {

		if (base) {
			NewS = NewS + S.c_str()[i];
			NewQ = NewQ + qual.c_str()[i];
		} else if (S.c_str()[i] != 'A' && S.c_str()[i] != 'C' &&
							 S.c_str()[i] != 'G' && S.c_str()[i] != 'T') {
		} else {
			base = true;
			NewS = NewS + S.c_str()[i];
			NewQ = NewQ + qual.c_str()[i];
		}
	}
	
	qual = NewQ;
	return NewS;
} 

// The reason this is done is because if a couple of contigs have a gap in them (because of illumina problem not haplotype)
// we won't collapse them on the first round without trimming these ends
string TrimLowCoverageEnds(string S, string& quals, string& depth, int cutoff) {
	bool base = false;
	string NewS = "";
	string NewD = "";
	string NewQ = "";

	for (int i = S.size() - 1; i >= 0; i--) {
		if (base) {
			NewS = S.c_str()[i] + NewS;
			NewD = depth.c_str()[i] + NewD;
			NewQ = quals.c_str()[i] + NewQ;
		} else if ((int)depth.c_str()[i] > cutoff) {
			base = true;
			NewS = S.c_str()[i] + NewS;
			NewD = depth.c_str()[i] + NewD;
			NewQ = quals.c_str()[i] + NewQ;
		}
	}
 

	string NewS2 = "";
	string NewD2 = "";
	string NewQ2 = "";
	if (NewS.size() > 1) {
		base = false;
		for (int i = 0; i < NewS.size(); i++) {
			if (base) {
				NewS2 = NewS2 + NewS.c_str()[i];
				NewD2 = NewD2 + NewD.c_str()[i];
				NewQ2 = NewQ2 + NewQ.c_str()[i];
			} else if ((int)NewD.c_str()[i] > cutoff) {
				base = true;
				NewS2 = NewS2 + NewS.c_str()[i];
				NewD2 = NewD2 + NewD.c_str()[i];
				NewQ2 = NewQ2 + NewQ.c_str()[i];
			}
		}
	}

	depth = NewD2;
	quals = NewQ2;
	return NewS2;
}


string AdjustBases(string sequence, string qual) {
	int MinQ = 5;
	int QualOffset = 33;
	string NewString = "";

	for (int i = 0; i < sequence.size(); i++) {
		if (qual.c_str()[i] - QualOffset < MinQ) {
			NewString += 'N';
		} else {
			NewString += sequence.c_str()[i];
		}
	}
	return NewString; 
}

bool replace(std::string& str, const std::string& from, const std::string& to) {
	size_t start_pos = str.find(from);
	if (start_pos == std::string::npos) { 
		return false;
	}
	str.replace(start_pos, from.size(), to);
	return true;
}

string FlipStrands(string strand) {
	string NewStrand = "";

	for (int i = 0; i < strand.size(); i++) {
		if (strand.c_str()[i] == '+') {
			NewStrand += "-";
		} else if (strand.c_str()[i] == '-') {
			NewStrand += "+";
		} else if (strand.c_str()[i] == '.'){
			NewStrand += ".";
		}
	}
	return NewStrand;
}
void compressStrand(string S, int& F, int& R) {
	for (int i = 0; i < S.size(); i++) {
		if (S.c_str()[i] == '+')
			F++;
		else if (S.c_str()[i] == '-')
			R++;
	}
	return;
}

bool parse_args(int argc, char* argv[], OverlapArgs& args) {

	// todo: test if this is correct number logic
	if (argc < 10) {
		cout << argc << " arguments provided, but at least 10 are required.\n";
		for (int i = 0; i < argc; i++) {
			cout << "Arg " << i << ": " << argv[i] << endl;
		}
		cout << "Usage: " << argv[0] << " <fastq_file> <MinPercent> <MinOverlap> "
						"<MinCoverage> <ReportStub> <hashLengthSize> <ACT> <OutFile> "
						"<LCendTrimLength> <Threads> [--verbose]\n";
		return false;
	}

	args.FastqIn = argv[1];
	args.MinPercent = stof(argv[2]);
	args.MinOverlap = strtol(argv[3], nullptr, 0);
	args.MinCoverage = strtol(argv[4], nullptr, 0);
	args.NameStub = argv[5];
	args.hashLength = strtol(argv[6], nullptr, 0);
	args.ACT = strtol(argv[7], nullptr, 0);
	args.OverlapStub = argv[8];
	args.TrimLCcuttoff = strtol(argv[9], nullptr, 0);
	args.Threads = strtol(argv[10], nullptr, 0);

	cout << "There were at least 10 args" << endl;
	return true;
}

int main(int argc, char* argv[]) {

	OverlapArgs args;
	if (!parse_args(argc, argv, args)) {
		cout << "Error overlap parsing arguments. Please check the usage." << endl;
		return 1; // Error in argument parsing
	}
	long int Buffer = 100;
	// long int Buffer = 100 * args.Threads; - THIS LEADS TO NON-DETERMINISM

	// Check & open file streams
	ifstream fastq;
	fastq.open(args.FastqIn.c_str());
	if (!fastq.is_open()) {
		cout << "Error, Fastq file could not be opened - " << args.FastqIn << endl;
		return -1;
	}

	ofstream report;
	std::stringstream ss;
	string FirstPassFile = args.FastqIn;
	ss << args.OverlapStub << ".fastq";
	FirstPassFile = ss.str();
	report.open(FirstPassFile.c_str());
	if (!report.is_open()) {
		cout << "Error, Mut-Output file could not be opened - " << FirstPassFile
				 << endl;
		return -1;
	}

	ofstream DepReport;
	FirstPassFile += "d";
    DepReport.open(FirstPassFile.c_str());
	if (!report.is_open()) {
		cout << "Error, Mut-Output depth file could not be opened - " << FirstPassFile
				 << endl;
		return -1;
	}

	ofstream good;
	FirstPassFile = ss.str();
	FirstPassFile += "good.fastq";
	good.open(FirstPassFile.c_str());
	if (!good.is_open()) {
		cout << "Error, Mut-Output good file could not be opened - " << FirstPassFile
				 << endl;
		return -1;
	}

	ofstream bad;
	FirstPassFile = ss.str();
	FirstPassFile += "bad.fastq";
	bad.open(FirstPassFile.c_str());
	if (!bad.is_open()) {
		cout << "Error, Mut-Output bad file could not be opened - " << FirstPassFile
				 << endl;
		return 0;
	}

	string line;
	vector<string> sequenes;	// The array of full-length sequences extracted from the input fastq file
	vector<string> qual;		// An array of the per-nucleotide qualities corresponding to the sequences
	vector<string> depth;		// An array of the per-nucleotide kmer-depths corresponding to the sequences
	vector<string> strand;		// An array of the strands each sequence is located on
	std::unordered_map<unsigned long, vector<int>> Hashes;	// The hash table of kmer hashes to indices of sequences containing that kmer
	int lines = -1;
	int goodlines = 0;
	int dup = 0;
	string L1;
	string L2;
	string L3;
	string L4;
	string L5;
	string L6;
	int Rejects = 0;
	string Fastqd = args.FastqIn;
	size_t found = Fastqd.find(".fastqd");

	// Read in entire fastq file, 6 lines at a time
	// If we're reading in a fastq+depth file, we simply trim off the low coverage ends before starting to process the reads
	if (found != string::npos) {
		int counter = 0;
		cout << "ATTENTION - Fastq+depth input detected, reading in FASTQD file \n";

		while (getline(fastq, L1)) {
			counter++;
			if (counter % 100 == 1) {
				cout << "Read in " << counter << " lines and rejected " << Rejects
						 << " reads\r";
			}

			getline(fastq, L2);
			getline(fastq, L3);
			getline(fastq, L4);
			getline(fastq, L5);
			getline(fastq, L6);
			string depths = "";
			int ReadSize = L2.size();
			bool Multiple = false;
			vector<string> temp = Util::Split(L6, ' ');

			for (vector<string>::size_type i = 0; i < temp.size(); i++) {
				unsigned char C = atoi(temp[i].c_str());
				depths += C;
				if ((int)C > 1) {
					Multiple = true;
				}
			}

			if (Multiple == true) {
				L2 = TrimLowCoverageEnds(L2, L4, depths, args.TrimLCcuttoff);
			}

			if (L2.size() > args.hashLength + 1) {
				lines++;
				sequenes.push_back(L2);
				qual.push_back(L4);
				depth.push_back(depths);
				strand.push_back(L5);
				ReadSize = L2.size();
			} else {
				Rejects++;
				bad << L1 << endl << L2 << endl << L3 << endl << L4 << endl;
			}
		}
	// If we have a fastq file, we (should be) checking for duplicates, trimming Ns, and adjustung
	//  bases based on quality values before processing reads
	} else {
		vector<string> DupCheck;
		cout << "Reading in raw fastq \n";
		int counter = 0;

		while (getline(fastq, L1)) {
			counter++;
			if (counter % 100 == 1) {
				cout << "Read in " << counter << " lines, rejected " << Rejects
						 << " reads with " << dup << " duplicates\r";
			}
			getline(fastq, L2);
			getline(fastq, L3);
			getline(fastq, L4);
			int Ns = 0;

			for (int i = 0; i < L2.size(); i++) {
				if (L2.c_str()[i] == 'N') {
					Ns++;
				}
			}

			lines++;
			int ReadSize = L2.size();
			string depths = "";
			ReadSize = L2.size() - 1;
			bool found = false;
			bool RunDupCheck = true;

			// if (RunDupCheck) {

			// 	// BUG FIX NEEDED
			// 	// NOTE: this is currently NEVER run because nothing added to DupCheck until we've already run the loop
			// 	#pragma omp parallel for num_threads(args.Threads) shared(DupCheck, L2, found)
			// 	for (int i = 0; i < DupCheck.size(); i++) {
			// 		if (L2.size() == DupCheck[i].size()) {
			// 			bool AllBasesMatch = true;

			// 			for (int k = 0; k < L2.size(); k++) {
			// 				if (L2.c_str()[k] == 'N' or DupCheck[i].c_str()[k] == 'N') {
			// 				} else if (L2.c_str()[k] == DupCheck[i].c_str()[k]) {
			// 				} else {
			// 					AllBasesMatch = false;
			// 					break;
			// 				}
			// 			}

			// 			if (AllBasesMatch) {
			// 				#pragma omp critical (found) 
			// 				{ 
			// 					found = true; 
			// 				}
			// 			}
			// 		}
			// 	}

			// 	if ((double)Ns / (double)L2.size() < 0.20) {
			// 		DupCheck.push_back(L2);
			// 	}
			// }

			if (found == false) {
				L2 = AdjustBases(L2, L4);
				L2 = TrimNends(L2, L4);

				if ((double)L2.size() / (double)ReadSize > .6) {
					goodlines++;
					sequenes.push_back(L2);
					qual.push_back(L4);
					strand.push_back("+");
					unsigned char C = 1;

					for (int i = 0; i <= ReadSize; i++) {
						depths += C;
					}

					depth.push_back(depths);
					good << L1 << endl << L2 << endl << L3 << endl << L4 << endl;
				} else {
					Rejects++;
				}
			} else {
				dup++;
				bad << L1 << endl << L2 << endl << L3 << endl << L4 << endl;
			}
		}
		DupCheck.clear();
	}


	good.close();
	bad.close();
	cout << "done reading " << endl;
	int NumReads = sequenes.size(); 
	cout << "\nDone reading in \n		 Read in a total of " << lines
			 << " and rejected " << Rejects << " with " << dup
			 << " duplicate reads detected for a total of " << goodlines
			 << "good reads" << endl;


	// First kmer table build after reading in all of the fastq/d reads
	RebuildHashTable(sequenes, 0, args.hashLength, Hashes, args.Threads);

	clock_t St, Et;
	int FoundMatch = 0;
	struct timeval start, end;
	gettimeofday(&start, NULL);
	int LinesSinceLastBuild = 1;
	int NumberHitPosLimit = 0;
	int NumberHitSanityLimit = 0;
	double AverageFPos = 0.0;
	double AverageFSanity = 0.0;
	double AverageRPos = 0.0;
	double AverageRSanity = 0.0;

	// Outer loop iterating through every sequence in chunks of Buffer size
	for (std::vector<string>::size_type b = 0; b < sequenes.size(); b += Buffer) {
		LinesSinceLastBuild += Buffer;

		// Rebuild hash table every million lines
		if (LinesSinceLastBuild > 1000000) { 
			RebuildHashTable(sequenes, b, args.hashLength, Hashes, args.Threads);
			LinesSinceLastBuild = 0;
		}

		vector<string> ToAddHashes;
		vector<int> ToAddPos;
		map<int, vector<int>> Forwards;
		map<int, vector<int>> Revs;
		int max = b + Buffer; // todo: check off by one errors here

		if (max > sequenes.size()) {
			max = sequenes.size();
		}

		if (FullOut) {
			cout << "Bulding list to align" << endl;
		}

		// For each sequence in chunk, prepare list of potential alignment matches by comparing kmers
		#pragma omp parallel for num_threads(args.Threads) shared(Hashes, Forwards)
		for (int i = b; i < max; i++) 
		{
			string A = sequenes[i];
			bool posLimit = false;
			bool sanityLimit = false;
			int NumPos = 0;
			int NumSanity = 0;
			PrepareSearchList(A, i, Hashes, args.hashLength, args.ACT, Forwards, posLimit, sanityLimit, NumPos, NumSanity);
			if (posLimit) {
				NumberHitPosLimit++;
			}
			if (sanityLimit) {
				NumberHitSanityLimit++;
			}
			AverageFPos =((AverageFPos * (double)b) + (double)NumPos) / ((double)b + 1.0);
			AverageFSanity = ((AverageFSanity * (double)b) + (double)NumSanity) /((double)b + 1.0);
		}

		#pragma omp parallel for num_threads(args.Threads) shared(Hashes, Revs)
		for (int i = b; i < max; i++) 
		{
			string A = Util::RevComp(sequenes[i]);
			bool posLimit = false;
			bool sanityLimit = false;
			int NumPos = 0;
			int NumSanity = 0;
			PrepareSearchList(A, i, Hashes, args.hashLength, args.ACT, Revs, posLimit,sanityLimit, NumPos, NumSanity);
			if (posLimit) {
				NumberHitPosLimit++;
			}
			if (sanityLimit) {
				NumberHitSanityLimit++;
			};
			AverageRPos =((AverageRPos * (double)b) + (double)NumPos) / ((double)b + 1.0);
			AverageRSanity = ((AverageRSanity * (double)b) + (double)NumSanity) /((double)b + 1.0);
		}

		if (FullOut) {
			cout << "Done Bulding List" << endl;
		}

		// Serially iterate through this chunk of sequences
		for (int i = b; i < max; i++) {
			string A, Aqual, Adep, Astr;
			A = sequenes[i];
			Aqual = qual[i];
			Adep = depth[i];
			Astr = strand[i];
			int k;
			bool PerfectMatch = false;
			int bestIndex = -1;
			float Dt;

			if (FullOut) {
				Et = clock();
				Dt = ((double)(Et - St)) / CLOCKS_PER_SEC;
				cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<**************************************************************>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
			}

			if (i % 100 == 5) {
				gettimeofday(&end, NULL);
				float Dt = end.tv_sec - start.tv_sec;
				cout << "done " << i << " of " << NumReads
						 << ", \% = " << ((double)i / (double)NumReads) * 100.00
						 << ", TT= " << Dt << " , S/R = " << Dt / i
						 << ", \% match= " << ((double)FoundMatch / (double)i) * 100.00
						 << " %Plimit = " << NumberHitPosLimit / 2 << ", "
						 << ((double)NumberHitPosLimit / 2.0) / (double)i
						 << " AvP= " << (AverageFPos + AverageRPos) / 2.0
						 << " % ILmit = " << NumberHitSanityLimit << ", "
						 << ((double)NumberHitSanityLimit / 2.0) / (double)i
						 << " AvI= " << (AverageRSanity + AverageFSanity) / 2.0 << "\r";
			}

			// Align the current sequence with the possible options in Forwards
			int bestScore = Align3(sequenes, A, Aqual, i, k, bestIndex, args.MinPercent, PerfectMatch, args.MinOverlap, Forwards[i], args.Threads, NumReads);

			if (FullOut) {
				cout << "best forward score is " << bestScore << " k is " << k
						 << " index = " << bestIndex << endl;
			}

			// If we don't have a perfect match, align the current sequence with the possible options in Revs
			if (!(PerfectMatch)) {
				string revA = Util::RevComp(A);
				string revAqual = Util::RevQual(Aqual);
				string revAdep = Util::RevQual(Adep);
				string revAstr = FlipStrands(Astr);
				int revk = -1;
				int revbestIndex = -1;

				if (FullOut) {
					cout << "Checking Reverse\n";
				}

				int revBestScore =	Align3(sequenes, revA, revAqual, i, revk, revbestIndex, args.MinPercent, PerfectMatch, args.MinOverlap, Revs[i], args.Threads, NumReads);
				
				if (FullOut) {
					cout << "best reverse score is " << revBestScore << " k is " << revk
							 << " index = " << revbestIndex << endl;
				}

				// TODO: here is one part where we need to keep revBestScore AND bestScore if they have equal alignment scores
				if (revBestScore > bestScore) {
					A = revA;
					Aqual = revAqual;
					Adep = revAdep;
					Astr = revAstr;
					k = revk;
					bestScore = revBestScore;
					bestIndex = revbestIndex;
				}
			} else {
				if (FullOut) {
					cout << "Perfect Match Found, Skipping Reverse Search" << endl;
				}
			}

			// Check that we meet our minimum overlap requirement
			if (bestScore < args.MinOverlap) {
				if (FullOut) {
					cout << "No good match found, skipping" << endl;
				}
			} else {
				string B, Bqual, Bdep, Bstr;
				B = sequenes[bestIndex];
				Bqual = qual[bestIndex];
				Bdep = depth[bestIndex];
				Bstr = strand[bestIndex];
				FoundMatch++;

				if (FullOut) {
					if (k > 0) {
						cout << "found match at " << k << endl;

						for (int z = 0; z < k; z++) {
							cout << "+";
						}

						cout << A << endl << B << endl;

						for (int z = 0; z < Bdep.size(); z++) {
							int bam = Bdep.c_str()[z];
							cout << bam;
						}

						cout << endl;

					} else {
						cout << "found match at " << k << endl;
						cout << A << endl;

						for (int z = 0; z < abs(k); z++) {
							cout << "-";
						}
						cout << B << endl;
						for (int z = 0; z < abs(k); z++) {
							cout << "-"; 
						}
						for (int z = 0; z < Bdep.size(); z++) {
							int bam = Bdep.c_str()[z];
							cout << bam;
						}
			
						cout << endl;
					}
				}

				// Collapse the sequences for the best match - again TODO: will need to make this work for multiple equal matches
				string combined = ColapsContigs(A, B, k, Aqual, Bqual, Adep, Bdep, Astr, Bstr);
				if (Bqual.size() != combined.size()) {
					cerr << "Error: something went wrong combining sequences into contigs" << endl;
				}

				qual[bestIndex] = Bqual;
				depth[bestIndex] = Bdep;
				sequenes[bestIndex] = combined;
				strand[bestIndex] = Bstr;
				sequenes[i] = "moved";

				// Update hash table Hashes with updated collapsed info
				// We might have new hashes here from combined sequence
				// #pragma omp parallel for num_threads(1) shared(Hashes)
				// for (int j = 0; j < combined.size() - args.hashLength; j++) {
				// 	string hash = combined.substr(j, args.hashLength);
				// 	size_t foundIdx = hash.find('N');
				// TODO: change this back, this was correct logic
					
				// 	if (foundIdx == std::string::npos) {
				// 		unsigned long forwardHash = Util::HashToLong(hash);
				// 		unsigned long reverseHash = Util::HashToLong(Util::RevComp(hash));
						
				// 		#pragma omp critical(updateHash) 
				// 		{	
				// 			bool foundForwardMatch = false;
				// 			vector<int>& forwardList = Hashes[forwardHash];
				// 			for (int k = 0; k < forwardList.size(); k++) {
				// 				if (forwardList[k] == bestIndex) {
				// 					foundForwardMatch = true;
				// 					break;
				// 				}
				// 			}
				// 			if (!foundForwardMatch) {
				// 				Hashes[forwardHash].push_back(bestIndex);
				// 			}
							
							
				// 			bool foundReverseMatch = false;
				// 			vector<int>& reverseList = Hashes[reverseHash];
				// 			for (int k = 0; k < reverseList.size(); k++) {
				// 				if (reverseList[k] == bestIndex) {
				// 					foundReverseMatch = true;
				// 					break;
				// 				}
				// 			}
				// 			if (!foundReverseMatch) {
				// 				Hashes[reverseHash].push_back(bestIndex);
				// 			}
				// 		}
				// 	}
				// }

				if (FullOut) {
					cout << combined << endl;
					cout << Bqual << endl;

					for (int z = 0; z < Bdep.size(); z++) {
						int bam = Bdep.c_str()[z];
						cout << bam;
					}

					cout << endl;
				}
			}
		}
	}

	cout << "\nRESULTS\n";
	int count = 0;

	// Iterate through each sequence from fastq in order of receipt
	// If the sequence is not "moved" and has sufficient length and coverage, write it to the report files
	// I don't think the index here should change between runs, since the sequences array is serially populated
	for (int i = 0; i < sequenes.size(); i++) {

		if (sequenes[i] != "moved" && sequenes[i].size() >= 95) {
			string rDep = depth[i];
			int maxDep = -1;

			for (int z = 0; z < rDep.size(); z++) {
				unsigned char bam = rDep.c_str()[z];
				if ((int)bam > maxDep) {
					maxDep = (int)bam;
				}
			}

			if (maxDep >= args.MinCoverage) {
				count++;
			 	int F = 0;
				int R = 0;
				compressStrand(strand[i], F, R);
				report << "@NODE_" << args.hashLength << "_" << i << "_L" << sequenes[i].size()<< "_D" << maxDep << ":" << F << ":" << R << ":" << endl;
				report << sequenes[i] << endl;
				report << "+" << endl;
				report << qual[i] << endl;

				DepReport << "@NODE_" << args.hashLength << "_" << i << "_L" << sequenes[i].size()<< "_D" << maxDep << ":" << F << ":" << R << ":" << endl;
				DepReport << sequenes[i] << endl;
				DepReport << "+" << endl;
				DepReport << qual[i] << endl;
				DepReport << strand[i] << endl;
				unsigned char C = depth[i].c_str()[0];
				int bestScore = C;
				DepReport << bestScore;

				for (int w = 1; w < depth[i].size(); w++) {
					C = depth[i].c_str()[w];
					bestScore = C;
					DepReport << ' ' << bestScore;
				}

				DepReport << endl;
			}
		}
	}
	cout << "Wrote " << count << " sequences" << endl;
	report.close();
	DepReport.close();
}
