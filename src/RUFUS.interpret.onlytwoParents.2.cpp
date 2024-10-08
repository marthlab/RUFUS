
/*This vertion imcorporates both copy number and mutation detection in 1
 *  * it needs to be run in two sptes, first the build, then the filter
 *   * it is split up to allow distribution to a cluster */

#include <unistd.h>
#include <ios>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unordered_map>
#include <bitset>
#include <time.h>
#include <sys/resource.h>
#include <stack>
#include <map>
#include <ctime>
#include <cmath>
#include "externals/fastahack/Fasta.h"
#include <fcntl.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>

#define NUMINTS  (1000)
#define FILESIZE (NUMINTS * sizeof(int))

using namespace std;

vector <unordered_map <unsigned long int, int> > ParentHashes;
unordered_map  <unsigned long int, int> MutantHashes; 
unordered_map <unsigned long int, int> ExcludeHashes;
FastaReference Reff;
int HashSize = 25; 
int totalDeleted; 
int totalAdded;
int MaxVarentSize = 1000; 
ofstream VCFOutFile;
ofstream BEDOutFile;
ofstream BEDBigStuff;
ofstream BEDNotHandled;
ofstream Invertions; 
ofstream Translocations; 
ofstream Translocationsbed; 
ofstream Unaligned; 
map <string, int> Hash; 
/////////////////////////
const vector<string> Split(const string& line, const char delim) {
    vector<string> tokens;
    stringstream lineStream(line);
    string token;
    while ( getline(lineStream, token, delim) )
	tokens.push_back(token);
    return tokens;
}
unsigned long HashToLong (string hash)
{
	bitset<64> HashBits;
	for(int i=0; i<hash.length();i++)
	{
		if (hash.c_str()[i] == 'A')
		{
			HashBits[i*2] = 0;
			HashBits[i*2+1] = 0;
		}
		else  if (hash.c_str()[i] == 'C')
		{
			HashBits[i*2] = 0;
			HashBits[i*2+1] = 1;
		}
		else  if (hash.c_str()[i] == 'G')
		{
			HashBits[i*2] = 1;
			HashBits[i*2+1] = 0;
		}
		else  if (hash.c_str()[i] == 'T')
		{
			HashBits[i*2] = 1;
			HashBits[i*2+1] = 1;
		}
		else
		{
			//cout << "ERROR, invalid character - " << hash.c_str()[i] << endl;
		}
	
	}
	return HashBits.to_ulong();
}

int checkPage( char *data, string hash, long int pageSize, string line)
{
	//cout << "checking page" << "with size = " << pageSize<< endl;
	bool firstNew = false;
	for (int i = 0; i <  pageSize; i++)
	{
	//      cout << i << endl;
	//      cout << data[i] << endl;
		if (data[i] == '\n')
		{
	//	      cout << "n found";
			if (firstNew != true)
				firstNew = true;
			else
			{
				vector<string> stuff;
				stuff = Split(line, '\t');
				string PageHash = stuff[0];
	//		      cout << "PageHash " << endl;
				if (hash == PageHash)
				{
	  //			      cout << "found a hash " << hash  << " - " << PageHash;
					return atoi(stuff[1].c_str());
				}
				line = "";
			}
		}
		else if (firstNew == true)
		{
	//	      cout << "first Newline found" << endl;
			line += data[i];
		}

	}
	return 0;
}

void ProcessPage(  char *data, string& PageFirstHash, string& PageLastHash, long int pageSize)
{
	string line = "";
	bool firstNew = false;\
	for (int i = 0; i < pageSize; i++)
	{
		if (data[i] == '\n')
		{
		
			if (firstNew == true)
				break;
			else
				firstNew = true;
		}
		else if (firstNew == true)
			line += data[i];
	}

	vector<string> stuff;
	stuff = Split(line, '\t');
	PageFirstHash = stuff[0];
	firstNew = false;
	line = "";
	for (int i = pageSize-1; i > 0; i+=-1)
	{
		if ( data[i] == '\n')
		{
			if (firstNew == true)
				break;
			else
				firstNew = true;
		}
		else if(firstNew == true)
			line = data[i] + line;
	}
	stuff = Split(line, '\t');
	PageLastHash = stuff[0];

}



int search(long int& fd, string hash, char* fileptr)
{
	//cout << "searching for " << hash << endl;
	char *data;
	struct stat sb;
	fstat(fd, &sb);

	long int pageSize;
	pageSize = sysconf(_SC_PAGE_SIZE);
	long int NumPages = sb.st_size/pageSize;
	//cout << "Number of pages = " << NumPages << endl;
 //       char *fileptr = NULL;

	long int off = 0;
	long int firstPos;
	long int lastPos;
	firstPos = 0;
	fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, firstPos*pageSize);
	data =  fileptr;
	// above should get me first page
	string FirstPageFirstHash;
	string FirstPageLastHash;
	ProcessPage(data, FirstPageFirstHash, FirstPageLastHash, pageSize);
	if (munmap(fileptr, pageSize) == -1) {
		perror("Error un-mmapping the file");
	}
	//cout << FirstPageFirstHash << endl << FirstPageLastHash << endl;
	//quck check to see if on first page
	if (hash >= FirstPageFirstHash and hash <= FirstPageLastHash)
	{
	 //     cout << "found on first page" << endl;
	 	fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, firstPos*pageSize);
		int val =  checkPage(data, hash, pageSize, "");
		if (munmap(fileptr, pageSize) == -1) {
			perror("Error un-mmapping the file");
		}
		return val; 
	}
	if (hash < FirstPageFirstHash)
	{
		cout << "HASH NOT IN FILE " << hash << endl;
		return 0; 
	}
	fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, pageSize*(NumPages-1));
	data = fileptr;
	lastPos = NumPages-1;
	//the above should get me the last two pages, we take two to ensure the last pages isnt just one character or something like that

	string LastPageFirstHash;
	string LastPageLastHash;
	ProcessPage(data, LastPageFirstHash, LastPageLastHash, pageSize);
	if (munmap(fileptr, pageSize) == -1) {
		perror("Error un-mmapping the file");
	}
	//cout << LastPageFirstHash << endl << LastPageLastHash << endl;
	//quck check to see if on last page
	if (hash >= LastPageFirstHash and hash <= LastPageLastHash)
	{
	      	fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, pageSize*(NumPages-1));
	//	cout << "found on last page" << endl;
		int val =  checkPage(data, hash, pageSize, "");
		if (munmap(fileptr, pageSize) == -1) {
			perror("Error un-mmapping the file");
		}
		return val; 
	}
	if (hash > LastPageLastHash)
	{
		cout << "HASH NOT IN FILE " << hash << endl;
		return 0;
	}
	//start the search
	int counter = 0;
	while (true)
	{
	 //     cout << "ON LOOP " << counter << endl << endl;
		counter++;
		long int currentPage = lastPos - ((lastPos-firstPos)/2);
	   //   cout << "checking page " << currentPage << " last = " << lastPos << " and first = " << firstPos << endl;;
		if (currentPage == lastPos or currentPage == firstPos or lastPos - firstPos < 3)
		{
			string extra = "";
		 //     cout << "\nenvoked this" << endl;
			fileptr = (char*)mmap64(NULL, pageSize*5, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, pageSize * (firstPos-1));
			data =  fileptr;
		     // cout << "made it here" << endl;
			int val = checkPage(data, hash, pageSize*5, extra);
			if (munmap(fileptr, pageSize*5) == -1) {
				perror("Error un-mmapping the file");
			}
			return val; 
		}
		//cout << " fileptr = (char*)mmap64(NULL, " << pageSize*2 <<", PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, " << pageSize<<" * " << currentPage <<");"<< endl;
		fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, pageSize*currentPage);
		data =  fileptr;
		string CurrentPageFirstHash;
		string CurrentPageLastHash;
		ProcessPage(data, CurrentPageFirstHash, CurrentPageLastHash, pageSize);
		if (munmap(fileptr, pageSize) == -1) {
			perror("Error un-mmapping the file");
		}
	//	cout << " with " << CurrentPageFirstHash << " and " << CurrentPageLastHash << endl;
		if (hash >= CurrentPageFirstHash and hash <= CurrentPageLastHash)
		{
			fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, pageSize*currentPage);
			int val = checkPage(data, hash, pageSize, "");
			if (munmap(fileptr, pageSize) == -1) {
			perror("Error un-mmapping the file");
			}
			return val;
		}
		else
		{
			if (hash < CurrentPageFirstHash)
			{
	  //		    cout << "hash " << hash << " is greater than " << CurrentPageFirstHash << " looking above" << endl;
				lastPos = currentPage;
				LastPageFirstHash = CurrentPageFirstHash;
				LastPageLastHash = CurrentPageLastHash;
			}
			else if (hash > CurrentPageLastHash)
			{
	    //		  cout << "hash \n" << hash << " is less than \n" << CurrentPageLastHash << " looking below" << endl;
				firstPos = currentPage;
				FirstPageFirstHash = CurrentPageFirstHash;
				FirstPageLastHash = CurrentPageLastHash;
			}
		}

	}
	close(fd); 

}

bool fncomp (char lhs, char rhs) {return lhs<rhs;}

struct classcomp {
  bool operator() (const char& lhs, const char& rhs) const
  {return lhs<rhs;}
};
///////////////////////////
int GetReadOrientation(int flag)
{
	// returns 0 = forward, 1 = reverse	
	bool b[16];
        
       	for (int j = 0;  j < 16;  ++j){
        	b [j] =  0 != (flag & (1 << j));	
	}
	cout << "bits are " << b[0] << " " << b[1] << " " << b[2] << " " << b[3] << " " << b[4] << " " << b[5] << " " << b[6] << " " << b[7] << " " << b[8] << " " << b[9] << " " << b[10] << " " << b[11] << " " << b[12] << endl;
	cout << "flag is " << flag << endl; 
	cout << "orientation is " << b[4] << endl; 
	return b[4];
	/*if (flag == 0 or flag == 2048){
		cout << "0 forward" << endl;
		return 0;}
	else if (flag == 16 or flag == 2064){
		cout << "1 reverse" << endl; 
		return 1;}
	else
		cout <<  "error flag not recognised " << flag << endl;

	return -1;*/
}
string getHash(string seq, int j, int HashSize)
{
	int bases = 0;
	string NewHash = "";
	if (j <seq.size())
	{
		if (seq.c_str()[j] == 'A' or seq.c_str()[j] == 'C' or seq.c_str()[j] == 'G' or seq.c_str()[j] == 'T') //if the first base is not in a real seqeunce, return a blank hash
		{
			while (bases<HashSize and j<seq.size())
			{
				if (seq.c_str()[j] == 'A' or seq.c_str()[j] == 'C' or seq.c_str()[j] == 'G' or seq.c_str()[j] == 'T')
				{
					NewHash+=seq.c_str()[j];
					bases++;
				}
				j++;
			}
			//cout << "grabbed hash " << NewHash << endl;
		}
	}
	return NewHash;
}

string RevComp (string Sequence)
{
	string NewString = "";
	//cout << "Start - " << Sequence << "\n";
	//cout << Sequence.length() << endl; 
	for(int i = Sequence.length()-1; i>=0; i+= -1)
	{
		char C = Sequence.c_str()[i];
	//      cout << C << endl;
		if (C == 'A')
			NewString += 'T';
		else if (C == 'C')
			 NewString += 'G';
		else if (C == 'G')
			 NewString += 'C';
		else if (C == 'T')
			 NewString += 'A';
		else if (C == 'N')
			 NewString += 'N';
		else
		{
			cout << "ERROR IN RevComp - " << C << " " ;
			NewString += C; 
		}

	}
	//cout << "end\n";
	return NewString;
}

void process_mem_usage(double& vm_usage, double& resident_set, double& MAXvm, double& MAXrss)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
	       >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
	       >> utime >> stime >> cutime >> cstime >> priority >> nice
	       >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
	if (vm_usage > MAXvm){MAXvm = vm_usage;}
	if (resident_set > MAXrss){MAXrss = resident_set;}
}

class SamRead
{
	public:
	string name;
	int flag;
	bool FlagBits[16];
	string chr;
	int pos;
	int mapQual;
	int AlignScore; 
	string cigar;
	string seq;
	string qual;
	string RefSeq; 
	string originalSeq;
	string originalQual; 
	string cigarString;
	string strand;  
	float StrandBias; 
	string strands; 
	int forward; 
	int reverse; 
	bool UsedForBigVar; 
	vector<int> alignments; 
	vector<int> Positions; 
	vector<string> ChrPositions;
	int AlignmentSegments; 

	vector<long> MutAltCounts; 
	vector<long> MutRefCounts; 
	vector<long> MutHashListCounts; 

	vector<vector <long>> RefAltCounts; 
	vector<vector <long>> RefRefCounts; 
	vector<string> AltKmers;
	vector<string> RefKmers;	
	bool first; // = true;  
	bool combined; // = false; 
	vector<bool> PeakMap; 

	void createPeakMap(); 
	void parse(string read);
	void getRefSeq();
	void CountAlignmentSegments(); 
	void processCigar();
	void parseInsertions( SamRead B); 
	void parseMutations( char *argv[] );
	void GetModes(int pos, string alt, string reff, int &MutRefMode, int &MutAltMode, vector <int> &ParRefModes, vector <int> &ParAltModes,  vector <int> &HashCounts, vector <int> &HashCountsOG, int &PossibleVarKmer);
	//string ShittyGenotyper(int Alt, int Ref); 	
	int GetSupportingHashCount(int pos, string alt,  string reff);
	void processMultiAlignment(); 
	void write();
	void writeVertical();
	void writetofile(ofstream &out);
	void flipRead(); 
	void LookUpKmers();
	void FixTandemRef();
	int CheckParentCov(int &mode); 
	bool StartsWithAlign(int &pos, string &insert); 
	bool EndsWithAlign(int &pos, string &insert); 
	bool StartsWithAlignAtPeak(int &pos, string &insert, int &Kdepth);
        bool EndsWithAlignAtPeak(int &pos, string &insert, int &Kdepth);
		
	bool CheckEndsAlign(); 
	int CheckBasesAligned();



	vector <string> hashes;
        vector <string> hashesRef;
        vector <bool> varHash;
	vector <bool> candidateHash; 
	vector <vector<int>> parentCounts;
        vector <vector<int>> parentCountsReference;
	vector <int> mutCounts;
        vector <int> mutCountsRef;
	void BuildUpHashCountTable(); 
};



void SamRead::BuildUpHashCountTable()
{
  /////////////////Building up varHash and hash lists ///////////// 
        cout << "Building up varHash" << endl;
        for (int i = 0; i <  seq.size() - HashSize; i++)
        {
                string newHash = "";
                string newHashRef = "";
                newHash += seq.c_str()[i];
                newHashRef += RefSeq.c_str()[i];
                int count = 0;
                ////can i replace this with get hash ?  
                if ((cigarString.c_str()[i] != 'D' and cigarString.c_str()[i] != 'R' and cigarString.c_str()[i] != 'H'))
                {
                        for (int j = 1; j<seq.size() - i and count < HashSize-1 ; j++)
                        {

                                if (cigarString.c_str()[i+j] != 'D' and cigarString.c_str()[i+j] != 'R'  and cigarString.c_str()[i+j] != 'H')
                                {
                                        newHash += seq.c_str()[i+j];
                                        newHashRef += RefSeq.c_str()[i+j];
                                        count++;
                                }
                        }
                }
                hashes.push_back(newHash);
                  hashesRef.push_back(newHashRef);
                if (Hash.count(newHash) > 0 or Hash.count(RevComp(newHash)) > 0)
                        varHash.push_back(true);
                else
                        varHash.push_back(false);
        }
        ///////////////////////////////////////////////////


        ///////////////////building up parent hash counts //////////////////
        cout << "Bulding Par hash counts" << endl;
        //vector <long int> ParentHash;
        for(int pi = 0; pi<ParentHashes.size(); pi++)
        {
                vector <int> counts;
                vector <int> countsRef;
                for(int i = 0; i< hashes.size(); i++)
                {
                        string hash =  hashes[i];
                        string hashRef = hashesRef[i];
                        bool checkHash = true;
                        for (int j = 0; j < HashSize; j++)
                        {
                                if (!(hash[j] == 'A' or hash[j] == 'C' or hash[j] == 'G' or hash[j] == 'T'))
                                {
                                        checkHash = false;
                                        break;
                                }
                        }
                        if (checkHash)
                        {
                                unsigned long int LongHash = HashToLong(hash);
                                if (ParentHashes[pi].count(LongHash) >0)
                                        counts.push_back(ParentHashes[pi][LongHash]);
                                else
                                        counts.push_back(0);
                                 unsigned long int LongHashRef = HashToLong(hashRef);
                                if (ParentHashes[pi].count(LongHashRef) >0)
                                        countsRef.push_back(ParentHashes[pi][LongHashRef]);
                                else
                                        countsRef.push_back(0);
                        }
                        else
                        {
                                counts.push_back(-1);
                                countsRef.push_back(-1);
                        }
                }
                parentCounts.push_back(counts);
                parentCountsReference.push_back(countsRef);
        }
         //////////////////////////////////////////////////////////////////
	                /////////////////////bulid Mut counts/////////////////////////////
        cout << "bulding mut counts" << endl; 
        cout << hashes.size() << endl;
        cout << hashesRef.size() << endl; 
        for(int i = 0; i< hashes.size(); i++)
        {
                cout << i<< endl;
                string hash =  hashes[i];
                cout << "   hash = " << hash << endl;
                string hashRef = hashesRef[i];
                cout << "RefHash = " << hashRef << endl;
                bool checkHash = true;
                for (int j = 0; j < HashSize; j++)
                {
                        if (!(hash[j] == 'A' or hash[j] == 'C' or hash[j] == 'G' or hash[j] == 'T'))
                        {
                                checkHash = false;
                                break;
                        }
                }
                cout << "check hash = " << checkHash << endl;
                if (checkHash)
                {
                        unsigned long int LongHash = HashToLong(hash);
                        if (MutantHashes.count(LongHash) >0)
                        {       mutCounts.push_back(MutantHashes[LongHash]);}
                        else
                        {       mutCounts.push_back(0);}

                        unsigned long int LongHashRef = HashToLong(hashRef);
                        if (MutantHashes.count(LongHashRef) >0)
                        {       mutCountsRef.push_back(MutantHashes[LongHashRef]);}
                        else
                        {       mutCountsRef.push_back(0);}
                }
                else
                {
                        mutCounts.push_back(-1);
                        mutCountsRef.push_back(-1);
                }
        }
        /////////////////////////////////////////////////////////////////////

	 ////////////////////write out vertical table/////////////////////////   
        cout << "writing hashes out vert" << endl;
        for(int i =0; i < hashes.size(); i++)
        {
                cout << i+pos << "\t" <<  i  << "\t" << hashes[i] << "\t" << varHash[i] << "\t" << PeakMap[i] << "\t" << (int) qual.c_str()[i]-33;
                cout << "\t" << "MutVar-" << mutCounts[i]; 
                for (int j = 0; j < parentCounts.size(); j++)
                {
                                cout << "\t" <<  parentCounts[j][i];
                }
                cout << "\t" << "MutRef-" << mutCountsRef[i];
                for (int j = 0; j < parentCountsReference.size(); j++)
                {
                                cout << "\t" <<  parentCountsReference[j][i];
                }
                cout << endl;

        }
        ////////////////////////////////////////////////////////////////
}
int SamRead::GetSupportingHashCount(int pos, string alt,  string reff)
{
				int Count =0; 
                                int lower = pos-HashSize;
                                if (lower < 0){lower =0;}
                                int upper = pos+alt.length()+reff.length();//-1;
                                cout << pos<<" + "<< alt.length() << " + " << reff.length() << endl;
                                if (upper > MutRefCounts.size()){
                                        cout << "this is going to break " << upper << " > " << MutRefCounts.size() << endl;
                                        upper = MutRefCounts.size();
                                }

                                for (int j = lower; j<upper; j++)
                                {
                                        if (Hash.count(AltKmers[j]) > 0 and Hash[AltKmers[j]] > 0)
                                                Count++; 
                                        else if (Hash.count(RevComp(AltKmers[j])) > 0 and Hash[RevComp(AltKmers[j])])
                                                Count++; 
                                }
	return Count; 
}
string ShittyGenotyper(int Alt, int Ref)
{
	if (Alt ==0 and Ref ==0)
		return ".";
	else if (Alt == 0 and Ref > 1)
		return "0/0";
	else if (Alt >0 and Ref ==0)
		return "1/1";
        else if ((double) Alt / ((double) Ref + (double) Alt)  >.85)
        	return "1/1";
        else if ((double) Alt / ((double) Ref + (double) Alt)  <.15)
                return "0/0";
        else
                return "0/1";
}
void SamRead::GetModes(int pos, string alt,  string reff, int &MutRefMode, int &MutAltMode, vector <int> &ParRefModes, vector <int> &ParAltModes,  vector <int> &HashCounts, vector <int> &HashCountsOG, int &PossibleVarKmer)
{
 				int lower = pos-HashSize+1;
        			if (lower < 0){lower =0;}
                                int upper = pos+alt.length()+reff.length()-1;
                                cout << pos<<" + "<< alt.length() << " + " << reff.length() << endl;
                                if (upper > MutRefCounts.size()){
                                        cout << "this is going to break " << upper << " > " << MutRefCounts.size() << endl;
                                        upper = MutRefCounts.size();
                                }

                                //////////////chekcing allele frequencies /////////// 
                                //vector <int> HashCountsOG;
                                vector <int> varMutRefCounts;
                                vector <int> varMutAltCounts;
                                vector<vector <int>> varParRefCounts;
                                vector<vector <int>> varParAltCounts;
                                vector <int> temp;
                                for(int pi = 0; pi<ParentHashes.size(); pi++){
                                        varParRefCounts.push_back(temp);
                                        varParAltCounts.push_back(temp);
                                }
                                string nonspecific = "nonspecific";
                                //vector <float> freqs;
                                cout << "checking NonSpecic Kmers";
				string LastAtlKmer = "boomba";
				for (int j = lower; j<upper; j++)
                                {
					cout <<  AltKmers[j] << endl << RefKmers[j] <<endl;
					if ( AltKmers[j] != RefKmers[j] and (ExcludeHashes[HashToLong(AltKmers[j])]<1 or ExcludeHashes[HashToLong(RevComp(AltKmers[j]))]<1) and AltKmers[j] != LastAtlKmer)
					{
						
						PossibleVarKmer++;
						cout << "Different" << endl;
					} 
					else
						cout << "SAME" << endl; 
					LastAtlKmer = AltKmers[j]; 
                                        if(MutRefCounts[j]>0 and AltKmers[j] != RefKmers[j] ) //and MutRefCounts[j]<200 and (ExcludeHashes[HashToLong(RefKmers[j])]<2 or ExcludeHashes[HashToLong(RevComp(RefKmers[j]))]<2)) //needs to be fixed, should be based on cov not cutoff of 200 
                                                varMutRefCounts.push_back(MutRefCounts[j]);
                                        if (MutAltCounts[j]>0 and AltKmers[j] != RefKmers[j] and MutAltCounts[j]<200 and (Hash.count(AltKmers[j]) > 0 or Hash.count(RevComp(AltKmers[j])) > 0)  and (ExcludeHashes[HashToLong(AltKmers[j])]<1 or ExcludeHashes[HashToLong(RevComp(AltKmers[j]))]<1)) //needs to be fixed, should be based on cov not cutoff of 200 
                                                varMutAltCounts.push_back(MutAltCounts[j]);
                                        
					for (int pi=0; pi < varParRefCounts.size(); pi++){
                                                if (RefRefCounts[pi][j] >0 and AltKmers[j] != RefKmers[j] )//and RefRefCounts[pi][j] <200 and (ExcludeHashes[HashToLong(AltKmers[j])]<2 or ExcludeHashes[HashToLong(RevComp(AltKmers[j]))]<2)) //needs to be fixed, should be based on cov not cutoff of 200 
                                                        varParRefCounts[pi].push_back(RefRefCounts[pi][j]);
                                                if (RefAltCounts[pi][j] >0 and AltKmers[j] != RefKmers[j] and RefAltCounts[pi][j] < 200 and (Hash.count(AltKmers[j]) > 0 or Hash.count(RevComp(AltKmers[j])) > 0)   and (ExcludeHashes[HashToLong(RefKmers[j])]<1 or ExcludeHashes[HashToLong(RevComp(RefKmers[j]))]<1)) //needs to be fixed, should be based on cov not cutoff of 200 
                                                        varParAltCounts[pi].push_back(RefAltCounts[pi][j]);

                                        }

                                        if (Hash.count(AltKmers[j]) > 0 and Hash[AltKmers[j]] > 0 and AltKmers[j] != RefKmers[j] )
                                                HashCountsOG.push_back(Hash[AltKmers[j]]);
                                        else if (Hash.count(RevComp(AltKmers[j])) > 0 and Hash[RevComp(AltKmers[j])] and AltKmers[j] != RefKmers[j] )
                                                HashCountsOG.push_back(Hash[RevComp(AltKmers[j])]);
                                        
					
					if (Hash[AltKmers[j]] > 0 and AltKmers[j] != RefKmers[j])
                                                HashCounts.push_back(Hash[AltKmers[j]]);
                                        else if (Hash[RevComp(AltKmers[j])] > 0 and AltKmers[j] != RefKmers[j])
						 HashCounts.push_back(Hash[RevComp(AltKmers[j])]);
					else
                                                HashCounts.push_back(-1);
                                }
                               // float freq = 0;
                               // if (freqs.size() > 0){
                               //         for (int i =0; i<freqs.size(); i++){
                               //                 freq+=freqs[i];
                               //         }
                               //         freq = freq/freqs.size();
                               //` }
                                //////////////////////////////////////////////

                                cout << "<><><><><>MutRef<><><><><><>" << endl ;
                                for (int s =0; s<varMutRefCounts.size(); s++){
                                        cout << varMutRefCounts[s] << " " ;
                                }
                                cout << endl << "<><><><><>MutAlt<><><><><><>" << endl;
                                for (int s =0; s<varMutAltCounts.size(); s++){
                                        cout << varMutAltCounts[s] << " " ;
                                }
                                cout << endl;

                                sort (varMutRefCounts.begin(), varMutRefCounts.end());
                                sort (varMutAltCounts.begin(), varMutAltCounts.end());

                                cout << "<><><><><>MutRefSorted<><><><><><>" << endl ;
                                for (int s =0; s<varMutRefCounts.size(); s++){
                                        cout << varMutRefCounts[s] << " " ;
                                }
                                cout << endl << "<><><><><>MutAltSorted<><><><><><>" << endl;
                                for (int s =0; s<varMutAltCounts.size(); s++){
                                        cout << varMutAltCounts[s] << " " ;
                                }
                                cout << endl;

                                for(int pi = 0; pi<varParRefCounts.size(); pi++){
                                         sort (varParRefCounts[pi].begin(), varParRefCounts[pi].end());
                                        cout << "<><><><><>Ref" << pi << "<><><><><<><>" << endl;
                                        for (int s =0; s<varParRefCounts[pi].size(); s++){
                                                cout << varParRefCounts[pi][s] << " " ;
                                        }
                                        cout << endl;
                                }
                                ///////////////////////////////////////////////
                                if (varMutRefCounts.size() >1)
                                        MutRefMode = varMutRefCounts[0];
					//MutRefMode = varMutRefCounts[(varMutRefCounts.size())/2]; //// switch this for line above to get the mode, right now were taking the min 
                                else if (varMutRefCounts.size() ==1)
					MutRefMode = varMutRefCounts[0]; 
				else
                                        MutRefMode = 0;

                                if (varMutAltCounts.size() >1)
                                         MutAltMode = varMutAltCounts[0];
					//MutAltMode= varMutAltCounts[(varMutAltCounts.size()-2)/2]; //// switch this for line above to get the mode, right now were taking the min 
                                else if (varMutAltCounts.size() ==1)
					MutAltMode = varMutAltCounts[0];
				else
                                        MutAltMode=0;

                                for(int pi = 0; pi<varParRefCounts.size(); pi++){
                                        if (varParRefCounts[pi].size() >1)
						ParRefModes.push_back(varParRefCounts[pi][0]);
                                                //ParRefModes.push_back(varParRefCounts[pi][((varParRefCounts[pi].size())/2)]);  //// switch this for line above to get the mode, right now were taking the min 
                                        else if (varParRefCounts[pi].size() ==1)
						ParRefModes.push_back(varParRefCounts[pi][0]);
					else
                                                ParRefModes.push_back( 0);
                                }
				for(int pi =0; pi<varParAltCounts.size(); pi++){
					if (varParAltCounts[pi].size()>1)
						ParAltModes.push_back(varParAltCounts[pi][0]);
						// ParAltModes.push_back(varParAltCounts[pi][((varParAltCounts[pi].size())/2)]); //// switch this for line above to get the mode, right now were taking the min 
					else if (varParAltCounts[pi].size()==1)
						ParAltModes.push_back(varParAltCounts[pi][0]);
					else
						ParAltModes.push_back( 0); 
				}
}
void SamRead::CountAlignmentSegments()
{
	AlignmentSegments = 0;
	char last = cigarString.c_str()[0]; 
	for (int i =1; i < cigarString.size(); i++)
	{
	  if (cigarString.c_str()[i] == 'M')
	  {}
	  else if (last == 'M')
	  {
	    AlignmentSegments++; 
	  }
	  last = cigarString.c_str()[i];
	}
	if (last == 'M')
	{
		AlignmentSegments++;
	}
}
int SamRead::CheckParentCov(int &mode)
{
	//vector<vector <long>> RefAltCounts;
	//vector<vector <long>> RefRefCounts;
	bool good = true;
	int lowC = 0;  
	vector<int> cov; 
	for (int pi = 0; pi < RefRefCounts.size(); pi++){
		for (int i = 0; i < RefRefCounts[pi].size(); i++){
			if (RefKmers[i] != ""){
				int ParRef = 0;
				int ParAlt = 0; 
				if (RefAltCounts[pi][i] > 0)
					ParAlt = RefAltCounts[pi][i]; 
				if (RefRefCounts[pi][i] > 0)
					ParRef = RefRefCounts[pi][i]; 
				cov.push_back(ParRef+ParAlt); 
				if (ParRef+ParAlt > 0 && ParRef+ParAlt < 10)
					lowC++; 
			}
		}
	}	
	if(cov.size()>1){
	
		sort (cov.begin(), cov.end());
		mode = cov[cov.size()/2];
	}
	else
		mode = -1; 

	return lowC; 
}

void SamRead::flipRead()
{
	cout <<"FLIPPING reads not on the same strand";
	write(); 
	string FlipSeq = ""  ;
	string FlipQual = ""   ; 
	string FlipRefSeq = "";
	string FlipCigarString = ""; 
	string FlipStrand = "";
	vector<bool> FlipPeakMap ; 
	vector<int> FlipPos;
	vector<string> FlipChrPos; 
	for (int i = seq.size() -1;  i >=0; i--)
	{
	//	FlipSeq += seq.c_str()[i]; 
		FlipQual += qual.c_str()[i];
	//	FlipRefSeq += RefSeq.c_str()[i];
		FlipCigarString += cigarString.c_str()[i];
		FlipStrand += '-'; 
		FlipPos.push_back(Positions[i]);
		FlipChrPos.push_back(ChrPositions[i]);
		FlipPeakMap.push_back(PeakMap[i]);
	}
	FlipSeq = RevComp(seq); 
	FlipRefSeq = RevComp(RefSeq);

	seq = FlipSeq; 
	qual = FlipQual; 
	RefSeq = FlipRefSeq;
	cigarString = FlipCigarString; 
	strand = FlipStrand; 
	Positions = FlipPos; 
	ChrPositions = FlipChrPos; 
	PeakMap=FlipPeakMap; 
	write(); 
	
}

void SamRead::processMultiAlignment()
{
	//check if this is a mis-joined contig

}
void SamRead::write()
{
	cout << name << endl;
	cout << "   flag = " << flag << endl;
	cout << "   mapQual = " << mapQual << endl;
	cout << "   Strand = " << GetReadOrientation(flag) << endl;
	cout << "   Alignments = " << alignments.size() << endl;
	cout << "   chr " << chr << " - " << pos << " qual = " << mapQual <<  endl;
	cout << "   cigar  = " << cigar << endl;
	cout << "   Seq    = " << seq << endl;
	cout << "   Qual   = " << qual << endl; 
	cout << "   Cigar  = " << cigarString << endl; 
	cout << "   RefSeq = " << RefSeq << endl;
	cout << "   strand = " << strand << endl;
	cout << "  PeakMap = "; 
	for (int i =0; i < PeakMap.size(); i++)
		{cout << PeakMap[i]; }
	cout << endl;
	cout << "   RefPositions: ";
	for (int i =0; i < Positions.size(); i++)
		cout << Positions[i] << " \t"; 
	cout << endl;
	cout << "   RefChromoso: "; 
	for (int i =0; i < ChrPositions.size(); i++)
	       	cout << ChrPositions[i] << " \t";
	cout << endl;
	
}
void SamRead::writetofile(ofstream &out)
{

	out << name << endl;
	out << "   flag = " << flag << endl;
	out << "   mapQual = " << mapQual << endl; 
        out << "   Strand = " << GetReadOrientation(flag) << endl;
	out << "   Alignments = " << alignments.size() << endl;
	out << "   chr " << chr << " - " << pos << " qual = " << mapQual <<  endl;
	out << "   AlignScore = " << AlignScore << endl; 
	out << "   cigar  = " << cigar << endl;
	out << "   Seq    = " << seq << endl;
	out << "   Qual   = " << qual << endl;
	out << "   Cigar  = " << cigarString << endl;
	out << "   RefSeq = " << RefSeq << endl;
	out << "   PeakMap= "; 
	for (int i =0; i < PeakMap.size(); i++)
	{
		out <<  PeakMap[i] ;
	}
	out << endl; 
	out << "   PMSize = " << PeakMap.size() << endl;
}

void SamRead::writeVertical()
{
	cout << "ParentHashes size = " <<  ParentHashes.size() << "RefAltCounts size " << RefAltCounts.size() << endl;
	cout << name << endl;
	cout << "   flag = " << flag << endl;
	cout << "   chr " << chr << " - " << pos << " qual = " << mapQual <<  endl;
	cout << "   cigar  = " << cigar << endl;
	cout << "   Alignments = " << alignments.size() << endl; 
	cout << "   AlignScore = " << AlignScore << endl; 
	for (int i =0; i < seq.size(); i++){

		cout << seq.c_str()[i] << "\t" << qual.c_str()[i] << "\t" << cigarString.c_str()[i] << "\t" << RefSeq.c_str()[i] << "\t" << ChrPositions[i] << "\t" << Positions[i] << "\t" << MutAltCounts[i] << "\t" << MutRefCounts[i] << "\t" << MutHashListCounts[i] << "\t" ;
		cout << "\tParents";
		for (int pi=0; pi < RefAltCounts.size(); pi++){
			cout << "\t" << RefAltCounts[pi][i] << "\t" << RefRefCounts[pi][i];
		}
		cout << "\t" << RefKmers[i] << "\t" << AltKmers[i]; 
		cout<< endl;
	}
}
string compressVar(string line, int start, string& StructCall)
{      
	cout << "compressing var" << endl; 
	char current = line.c_str()[0];
	int currentCount = 1;
	string CV = ""; 
	for (int i = 1; i< line.size(); i++)
	{       
		cout << current << endl;
		if (line.c_str()[i] == current)
		{       
			currentCount++;
		}
		else    
		{       
			if (currentCount > 2)
			{       
				ostringstream convert;
				convert << currentCount; 
				CV+=convert.str();
				CV+=current;
				
				ostringstream convertEND; 
				int end = currentCount+start; 
				convertEND << end; 

				
				if (current == 'Y')
				{
					cout << "YAAAY STRUCT" << endl; 
					StructCall = "SVTYPE=DUP;END="; 
					StructCall += convertEND.str(); 
					StructCall += ";SVLEN="; 
					StructCall += convert.str(); 
					StructCall += ";";
					cout << StructCall << endl; 
				}
			}
			else if (currentCount ==2)
			{       
				CV+=current;
				CV+=current;
			}
			else if (currentCount == 1)
			{       
				CV+=current;
			}
			else    
			{       
				cout << "ERROR in compress " << current << " " << currentCount << endl;
			}
			
			current = line.c_str()[i];
			currentCount = 1;
		}
	}
	if (currentCount > 2)
	{       
		ostringstream convert;
		convert << currentCount;
		CV+=convert.str();
		CV+=current;
		
		ostringstream convertEND;
		int end = currentCount+start;
		convertEND << end;


		if (current == 'Y')
		{
			cout << "YAAAY STRUCT" << endl;
			StructCall = "SVTYPE=DUP:TANDEM;END=";
			StructCall += convertEND.str();
			StructCall += ";SVLEN=";
			StructCall += convert.str();
			StructCall += ";";
			cout << StructCall << endl;
		}
	}
	else if (currentCount ==2)
	{       
		CV+=current;
		CV+=current;
	}
	else if (currentCount == 1)
	{       
		CV+=current;
	}
	else
	{
		cout << "ERROR in compress " << current << " " << currentCount << endl;
	}
	
	return CV;
}
void SamRead::createPeakMap()
{
	vector<bool> tempPeakMap;
	int i =0; 
	int max = -1; 
	int maxSpot = -1; 
	int last = 0; 
	for (int i =0; i< qual.size(); i++)
	{

		if (qual[i] <='!')
                {
                        //cout << "!"; 
                        tempPeakMap.push_back(0);
                }
                else
                {
                        //cout <<endl << qual[i];
                        int j = i; 
                        max = qual[j]; 
                        while ( j < qual.size() and qual[j] > '!' )
                        {
                        //      cout << qual[j] << "-" << j ; 
                                if (max < qual[j])
                                {max = qual[j];}
                                j++;
                        //      cout << "max = " << (char) max << endl;
                        }
                        cout << endl;
                        j = j-1;
                        for ( int k = i; k < qual.size() and k <= j ; k++)
                        {
                        //      cout << qual[k]; 
                                if (qual[k]==max and cigarString[k] != 'H')
                                        tempPeakMap.push_back(1);
                                else
                                        tempPeakMap.push_back(0);
                        }
                        //cout << endl; 
                        //not sure why I need this, figure it out 
                        //tempPeakMap.push_back(0);
                        i = j;
                }
		/*if (qual[i] <='!')
		{
			//cout << "!"; 
			tempPeakMap.push_back(0);
		}
		else
		{
			//cout <<endl << qual[i];
			int j = i; 
			max = qual[j]; 
			while ( j < qual.size() and qual[j] > '!' )
			{
			//	cout << qual[j] << "-" << j ; 
				if (max < qual[j])
				{max = qual[j];}
				j++;
			//	cout << "max = " << (char) max << endl;
			}
			cout << endl;
			j = j-1;
			for ( int k = i; k < qual.size() and k <= j ; k++)
			{
			//	cout << qual[k]; 
				if (qual[k]==max and cigarString[k] != 'H')
					tempPeakMap.push_back(1);
				else
					tempPeakMap.push_back(0);
			}
			//cout << endl; 
			//not sure why I need this, figure it out 
			//tempPeakMap.push_back(0);
			i = j;
		}*/
	}

	// I hate one time corrections, but here on is to correct if ther is a del 
	for (int i =0; i< qual.size(); i++)
        {
		if (seq[i] == '-'){
			tempPeakMap[i] == tempPeakMap[i-1];
		}
	}

	PeakMap.clear();	
	PeakMap = tempPeakMap;
	//cout << "done with peak map " << endl;  
}
/*void SamRead::createPeakMap()
{
//	cout << "crateing PeakMap" << endl;
	vector<bool> tempPeakMap;
	int i =0; 
	int max = -1; 
	int maxSpot = -1; 
	while (i<qual.size())
	{
		int j;
		for ( j = i+1; j < qual.size()-1; j++)
		{
			if(qual.c_str()[j+1] < qual.c_str()[j])
			{
				max =  qual.c_str()[j];
				while (qual.c_str()[j+1] < qual.c_str()[j] && j < qual.size()-1)
				{
					j++;
				}
				break;
			}
		}
		int k;
		for ( k = i; k <=j; k++)
		{
			if (qual.c_str()[k] == max)
			       tempPeakMap.push_back(1);
			else
				tempPeakMap.push_back(0);
		}
		i = i +k;
	}
	for (int i = tempPeakMap.size()-1; i < qual.size(); i++)
	{
		tempPeakMap.push_back(0); 
	}
	PeakMap.clear();	
	PeakMap = tempPeakMap; 
//	for (int i =0; i< PeakMap.size(); i++)
//		cout << i << "-" << PeakMap[i] << endl;

}*/
void SamRead::parseMutations( char *argv[])
{
	cout << "In Parsing Mutations " << endl;
	BuildUpHashCountTable();
	createPeakMap();
	write(); 
	
	
	string StructCall = ""; 
	string reff = "";
	string alt = "";
	string varType = "";
	for(int i = 0; i<cigarString.size(); i++)
	{
		reff = ""; 
		alt = "";
		varType = "";
		//find the first variant base
		if ((cigarString.c_str()[i] == 'X' or cigarString.c_str()[i] == 'I' or cigarString.c_str()[i] == 'D' or cigarString.c_str()[i] == 'Y'/*  or cigarString.c_str()[i] == 'S' *or cigarString.c_str()[i] == 'H'*/) and RefSeq.c_str()[i] != 'N')
		{
			cout << "found a " << cigarString.c_str()[i] << endl;
			cout << "at pos " << i << " so pos " << i+pos << endl; 
			int size = -1; 
			int startPos = i; 
			bool  AnyBasesOver0  = false; 
			string Denovo = "inherited"; 
			
			if (qual.c_str()[i] > '!')
				AnyBasesOver0 = true;
			if (PeakMap[i] == 1)
				Denovo = "DeNovo";
			
			for(int j = 0; j< cigarString.size() - i; j++)
			{
				if(cigarString.c_str()[i+j] == 'X' or cigarString.c_str()[i+j] == 'D' or cigarString.c_str()[i+j] == 'I' or cigarString.c_str()[i+j] == 'Y' /*or cigarString.c_str()[i+j] == 'S' or cigarString.c_str()[i+j] == 'H'*/)
				{
					size = j; 
					if (qual.c_str()[i+j] > '!')
						AnyBasesOver0 = true;
					if (PeakMap[i+j] == 1)
						Denovo = "DeNovo";
	
				}
				else //if (qual.c_str()[i+j] == '!')
					break; 
			}
			cout << "size =" << size<< endl;
			
			if (AnyBasesOver0)  //enabling this will only report varites covered by hashes 
			{
				
				if ( cigarString.c_str()[i] == 'I' or cigarString.c_str()[i] == 'D' or cigarString.c_str()[i] == 'Y' /*or cigarString.c_str()[i] == 'S' or cigarString.c_str()[i] == 'H'*/)
				{
					for (int k = 1; i-k >= 0; k++)
					{
						if (ChrPositions[i-k] == "nope")
						{}
						else
						{
							reff+=RefSeq.c_str()[i-k]; 
							alt+=seq.c_str()[i-k];
							startPos = i-k;
							break;
						}
					}  
				}
				
				/////////build the alleles and var type/////////
				for(int j = 0; j<= size; j++)
				{
					if (RefSeq.c_str()[i+j] == 'A' or RefSeq.c_str()[i+j] == 'C' or RefSeq.c_str()[i+j] == 'G' or RefSeq.c_str()[i+j] == 'T')
						reff+=RefSeq.c_str()[i+j]; 
					if (seq.c_str()[i+j] == 'A' or seq.c_str()[i+j] == 'C' or seq.c_str()[i+j] == 'G' or seq.c_str()[i+j] == 'T')
						alt+=seq.c_str()[i+j]; 
					varType += cigarString.c_str()[i+j]; 
				}
				//***********check that the alese are only baess**************
				bool good = true; 
				for (int j = 0; j<reff.size(); j++)
				{	if (reff.c_str()[j] != 'A' or reff.c_str()[j] != 'C' or reff.c_str()[j] != 'G' or reff.c_str()[j] != 'T'){good = false;}}
				for (int j = 0; j<alt.size(); j++)
				{	if (alt.c_str()[j] != 'A' or alt.c_str()[j] != 'C' or alt.c_str()[j] != 'G' or alt.c_str()[j] != 'T'){good = false;}}
				if (good = false)
				{
					cout <<"ERROR in SNP detect"<< endl;
					cout <<endl<< chr << "\t" << pos+i << "\t" << reff << "\t" << alt << endl;
					write();
				}
				 //***********check that the alese are only baess done**************
			        ////////Starting Checks /////////////     
                                int MutRefMode;
                                int MutAltMode;
                                vector <int> ParRefModes;
                                vector <int> ParAltModes;
                                vector <int> HashCounts;
                               	vector <int> HashCountsOG;
                                int PossibleAltKmer=0; 
				GetModes(i, alt, reff, MutRefMode, MutAltMode, ParRefModes, ParAltModes, HashCounts, HashCountsOG, PossibleAltKmer);
                        	int SupportingHashes = GetSupportingHashCount(i, alt,  reff);
				string Genotype = ShittyGenotyper(MutAltMode, MutRefMode); 
				string CompressedVarType = compressVar(varType, Positions[startPos], StructCall); 
					cout <<  chr << "\t" << pos+i << "\t" << CompressedVarType /*"."*/ << "\t" << reff << "\t" << alt << "\t" << SupportingHashes << "\t" << varType << "\t" << "." << "\t" << "." << "\t" << "." << endl;
				////////////////generatre parent genotypes and check///////////////////////
				vector <string> ParGenotypes; 
				for (int p = 0; p< ParRefModes.size(); p++)
                                {
                                        ParGenotypes.push_back(ShittyGenotyper(ParAltModes[p], ParRefModes[p]) );
                                }
				
                                cout << endl;
				////////////////check that parents have enough coverage////////////////////
				cout << "PAR LOW COV CHECK" << endl; 
				int NumLowCov = 0; 
				int low = i-HashSize-10; 
				if (low < 0)
				  low = 0; 

				for(int k = low ; k <= i+10 and k < hashes.size(); k++)
        			{
               				for (int j = 0; j < 2 /*parentCounts.size()*/; j++)
                			{
						int sum = 0; 
						if (hashesRef[k] == hashes[k])
                                                {sum = parentCountsReference[j][k];cout <<hashesRef[k]<< "\t" << hashes[k] << "\tParSame" << sum; }
                                                else
                                                {sum = parentCounts[j][k] + parentCountsReference[j][k]; cout <<hashesRef[k]<< "\t" << hashes[k] << "\tParDiff" << sum;  }
						
						if (sum <= 7 and parentCounts[j][k] + parentCountsReference[j][k] > 0 )
						{
							NumLowCov++;
							 cout << "\tLOWCOV" << NumLowCov ;
						}
					}
                			cout << endl;

        			}
				//////////////////////////////check if the parenst contain any of mut hashes//////////////////////////////////////////
				bool LowCov = false; 
                        	int lowCount = 0;
                        	low = i - HashSize ; 
                        	if (low < 0){low = 0;}
                        	cout << "checking bases " << low << " to " << i+size+5 << endl;
                        	for(int j = low; j <= i+size and j < hashes.size(); j++)
                        	{
					if (hashesRef[j] != hashes[j])
					{
					  for (int k = 0; k < 2 /*parentCounts.size()*/; k++)
                               	          {
                                	                cout << "Checking Par Hash " << hashes[j] << "\t" << parentCounts[k][j]  << "\t" << hashesRef[j] << "\t" << parentCountsReference[k][j];
							float varFreq = 1; 
							if (parentCountsReference[k][j] > 0)
							{
							  varFreq = (double)parentCounts[k][j]/((double)parentCountsReference[k][j] + (double)parentCounts[k][j]); 
							}
							cout << "\tvarFreq=" << varFreq<<endl; 
                                	                if (parentCounts[k][j] >= 1 and parentCounts[k][j] <= 10 and varFreq > .02 )//and parentCountsReference[k][j]<150 ) //if (parentCounts[k][j] <= 5 and parentCounts[k][j] > 0 )
                                	                {
                                	                  if ((ExcludeHashes[HashToLong(hashes[j])]<1 and ExcludeHashes[HashToLong(RevComp(hashes[j]))]<1))
                                	                  {
                               	                     		cout << "      LC HASH FOUND" << endl;
                               	                     		LowCov = true;
                                	                    	lowCount++;
                                	                  }
                                	                  else
                                	                  { cout << "      HASH FOUND IN REF " << hashes[j]  << " - " << ExcludeHashes[HashToLong(hashes[j])]  << " - " << ExcludeHashes[HashToLong(RevComp(hashes[j]))] << endl;}
                                	                }
                                	  }
                                	}
                        	}
				///////////////////final filter check/////////////////////////////////////////
				string Filter = ".";
				string InfoFilter = ""; 
 				if (Genotype.find("1") == std::string::npos) {
                                        Denovo = "Mosaic";
                                } 
				if (AlignmentSegments > 10)
				{
				  	Denovo = "PoorAlignment"; 
					stringstream ss;
					ss << AlignmentSegments; 
					Denovo += ss.str();
					if (Filter == ".")
						Filter = "";  
					Filter+="PA";
					InfoFilter+="PA";
					InfoFilter+=ss.str();
					InfoFilter+=",";
					Filter+=";";
				}
				if (NumLowCov > 3)
				{
					Denovo = "ParLowCovRegion";
					cout << "ParLowCov " << NumLowCov << endl;
					if (Filter == ".")
                                                Filter = "";
					stringstream ss;
					ss << NumLowCov; 
					Filter+="PLC";
					InfoFilter+="PLC"; 
                                        InfoFilter+=ss.str();
                                        Filter+=";";  
					InfoFilter+=",";
				}
				//if (LowCov)
				if (lowCount >=2)
				{
					cout << "LOW COVERAGE" << endl;
					Denovo = "LowCov";
					stringstream ss;
					ss << lowCount;
					Denovo += ss.str(); 
					if (Filter == ".")
                                                Filter = "";
					Filter+="LCH";
					InfoFilter+="LCH";
                                        InfoFilter+=ss.str();
                                        Filter+=";";
					InfoFilter+=",";
				}
				else
				   	cout << "GOOD COVERAGE" << endl;
				if (StrandBias >= 0){

                                        if (StrandBias >0.9 or StrandBias < 0.1)
                                        {        
						Denovo = "StrandBias";
						stringstream ss;
                                        	ss << StrandBias; 
						if (Filter == ".")
							Filter = "";
						Filter+="SB";
						InfoFilter+="SB";
                                        	InfoFilter+=ss.str();
                                        	Filter+=";";
						InfoFilter+=",";
					}
                                }	 
				if (Denovo == "DeNovo" and Filter == ".")
					Filter = "PASS";
				if (InfoFilter=="")
					InfoFilter="PASS"; 

				for (int p = 0; p< ParRefModes.size(); p++)
				{
					//if (ParGenotypes[p].find("1") != std::string::npos) {
						 //Denovo = "PresentInParents"; }
				}
				/////////////////////////////////////////////////////// 
				 cout << "startpos = " << startPos << " chrsize = " << ChrPositions.size() << endl; 	
				 cout << ChrPositions[startPos] << "\t" << endl;
				 cout << Positions[startPos] << "\t"  << endl;
				 cout << CompressedVarType <<"-"  << endl;
				 cout << Denovo /*"."*/  << "\t"  << endl;
				 cout << reff << "\t"  << endl;
				 cout << alt << "\t"  << endl;
				 cout << SupportingHashes << "\t"  << endl;
				 cout << Filter << "\t"  << endl;
				 cout << StructCall  << endl;
				 cout <<"RN=" << name  << endl;
				 cout << ";MQ=" << mapQual  << endl;
				 cout << ";cigar=" << cigar  << endl;
				 cout << ";" << "CVT=" << CompressedVarType << ";HD=" << endl;

				double Score = ((double)SupportingHashes/(double)PossibleAltKmer) * 100.0; 
			////////////////////////Writing var out to file/////////////////////////
				cout       << ChrPositions[startPos] << "\t" <<Positions[startPos] << "\t" << CompressedVarType <<"-" <<Denovo /*"."*/  << "\t" << reff << "\t" << alt << "\t" << SupportingHashes << "\t" << Filter << "\t" << StructCall <<"FEX=" << InfoFilter << ";RN=" << name << ";MQ=" << mapQual << ";cigar=" << cigar << ";" << "CVT=" << CompressedVarType << ";HD="; 
					
				VCFOutFile << ChrPositions[startPos] << "\t" <<Positions[startPos] << "\t" << CompressedVarType <<"-" << Denovo /*"."*/  << "\t" << reff << "\t" << alt << "\t" << Score << /*SupportingHashes << "/" << PossibleAltKmer << */ "\t" << Filter << "\t" << StructCall << "FEX=" << InfoFilter << ";FS=" << SupportingHashes << "/" << PossibleAltKmer << ";RN=" << name << ";MQ=" << mapQual << ";cigar=" << cigar << ";SB=" << StrandBias << ";" << "AS=" << AlignmentSegments << ";" << "CVT=" << CompressedVarType << ";HD="; 

				for (int j = 0; j < HashCounts.size(); j++) 
				{	
					cout       << HashCounts[j] << "_";
					VCFOutFile << HashCounts[j] << "_"; 
				}
				VCFOutFile << ";AO=" << MutAltMode;
				cout       << ";VT=" <<  varType << "\t" ;
				VCFOutFile <<  ";VT=" <<  varType << "\t" ;
				cout       << "GT:DP:RO:AO" << "\t" << Genotype << ":" << MutRefMode + MutAltMode << ":" << MutRefMode << ":" << MutAltMode ;
				VCFOutFile << "GT:DP:RO:AO" << "\t" << Genotype << ":" << MutRefMode + MutAltMode << ":" << MutRefMode << ":" << MutAltMode ; 
				int ParentMode; 
				int lowC = CheckParentCov(ParentMode); 
				//      cout << ":" << lowC << ":" << ParentMode << ":" << StrandBias;; 
				//VCFOutFile << ":" << lowC << ":" << ParentMode << ":" << StrandBias;;
				for (int p = 0; p< ParRefModes.size(); p++)
				{
					VCFOutFile << "\t" << ParGenotypes[p] << ":" << ParAltModes[p] + ParRefModes[p] << ":" << ParRefModes[p] << ":" << ParAltModes[p]; 
				}
				cout << endl; 
				VCFOutFile << endl;
					
				BEDOutFile << chr << "\t" << pos+i << "\t" <<  pos+i+size << "\t" << chr << ":" << pos+i << ":" << (int)(reff.length() - alt.length()) << ":" << SupportingHashes << endl;
				i+=size;
				

				cout << "\nModes\t" << ChrPositions[startPos] << "\t" << Positions[startPos] << "\t" << CompressedVarType <<"-" <<Denovo /*"."*/  << "\t" << reff << "\t" << alt <<"\t" << MutRefMode <<"\t" << MutAltMode;
                                for(int pi = 0; pi<ParentHashes.size(); pi++){
                                         cout << "\t" << ParRefModes[pi];
                                         }
                                cout << endl;
			} 
		}
	}
		cout << "Out parsingMutations" << endl;
}
void SamRead::processCigar()
{
	string num = ""; 
	for (int i = 0; i < cigar.length(); i++)
	{
		if (cigar.c_str()[i] >= 48 and cigar.c_str()[i] <= 57)
			num = num + cigar.c_str()[i]; 
		else 
		{
			int number = atoi(num.c_str()); 
			for(int j = 0; j < number; j++)
			{cigarString += cigar.c_str()[i];}
			num = ""; 
		}
	}
	
}
void SamRead::FixTandemRef()
{
	cout << "FOUND TANDEM" << endl;
	write();
	//writeVertical();  
	string lastChr = "nope"; 
	int lastPos = -1; 
	string NewRef = ""; 
	for (int i = 0; i<seq.size(); i++){
		cout << i << " " << seq[i] << " " <<  cigarString[i] << " " << RefSeq[i] << endl;
		if (cigarString[i]=='Y'){
			NewRef+= toupper(Reff.getSubSequence(lastChr, lastPos+1, 1).c_str()[0]); 
			lastPos++; 
		}
		else{
			NewRef+=RefSeq[i]; 
			lastChr = ChrPositions[i]; 
			lastPos = Positions[i];
		}
	}
	RefSeq = NewRef; 
	cout << "done tandtem fix" << endl;
	write();
}
void SamRead::getRefSeq()
{
	originalSeq = seq; 
	originalQual = qual; 
	RefSeq = ""; 
	string NewSeq = "";
	string NewQual = "";
	string NewCigar = "";
	string NewStrand = ""; 
	vector<int> NewPositions; 
	vector<string> NewChromosome;
	int InsOffset = 0;
	
	if ( Reff.sequenceNameStartingWith(chr) == "") //come back to, need to check if chr is in reference 
	{
		cout << "ERROR chr " << chr << " not found\n"; 
		return; 
	}

   //correct star position of the read to account for Hard and soft clipped bases as we are counting those now
	for (int i = 0; i< cigarString.size(); i++)
	{
		if (cigarString.c_str()[i]== 'H')
		{}
		else
		{
			pos = pos-i;
			break;
		}
	}
	for (int i = 0; i< cigarString.size(); i++)
      	{
		if (cigarString.c_str()[i]== 'S')
		{}
		else
		{
			pos = pos-i;
			break;
		}
	}
	



	int Roffset = 0; 
	int Coffset = 0; 
	for (int i =0; i<cigarString.length(); i++)
	{
		if (cigarString.c_str()[i] == 'M')
		{
			//cout << "yay cigar = M" << endl;
			RefSeq += toupper(Reff.getSubSequence(chr, i+pos-1+Roffset, 1).c_str()[0]);
			NewSeq += seq.c_str()[i-Coffset]; 
			NewQual += qual.c_str()[i-Coffset];
			NewPositions.push_back(pos+i-InsOffset);
			NewChromosome.push_back(chr);
			if (   toupper(Reff.getSubSequence(chr,i+pos-1+Roffset, 1).c_str()[0]) == seq.c_str()[i-Coffset])
				NewCigar +='M'; 
			else
				NewCigar +='X';

		}
		else if(cigarString.c_str()[i] == 'I')
		{
			// cout << "yay cigar = I" << endl;
			RefSeq+='-'; 
			Roffset+= -1;
			NewSeq += seq.c_str()[i-Coffset];
			NewQual += qual.c_str()[i-Coffset];
			NewCigar += 'I';
			InsOffset++; 
			NewPositions.push_back(pos+i-InsOffset);
			NewChromosome.push_back(chr);
		}
		else if(cigarString.c_str()[i] == 'D')
		{
			// cout << "yay cigar = D" << endl;
			NewSeq += '-';
			NewQual += ' '; 
			Coffset++;
			RefSeq += toupper(Reff.getSubSequence(chr, i+pos-1+Roffset, 1).c_str()[0]);
			NewCigar += 'D';
		 	NewPositions.push_back(pos+i-InsOffset);
			NewChromosome.push_back(chr);
		}

		else if(cigarString.c_str()[i] == 'H')
		{
			RefSeq += 'H';
			NewSeq += 'H';
			NewQual += ' '; 
			Coffset++; 
			NewCigar+='H';
		//	Roffset+= -1;
			NewPositions.push_back(-1);
			NewChromosome.push_back("nope");
		}
		else if(cigarString.c_str()[i] == 'S')
		{
			//I have an error here, when a read starts with S its ofset wrong
			RefSeq += '-';
			NewSeq += seq.c_str()[i-Coffset];
			NewQual += qual.c_str()[i-Coffset]; 
			NewCigar+='S';
		//	Roffset += -1; 
			NewPositions.push_back(pos+i-InsOffset); //NewPositions.push_back(-1);
			NewChromosome.push_back(chr); //NewChromosome.push_back("nope");
		}
		//{ cout << "yay cigar = H" << endl;}
		else{cout << "well shit" << endl;}
		//cout << "yay" << endl;
	}
	seq = NewSeq;
	cigarString = NewCigar; 
	qual.clear();
	char lastQ = ' ';

	//cout << "qual = " << qual; 
	for (int i =0; i<NewQual.size(); i++)
	{       
	 	if (NewQual.c_str()[i] == ' ')
		{
			if (NewCigar.c_str()[i] == 'D')
				qual += lastQ; //'!';
			else
				qual += '!';
		} 
		else 
		{
			qual+=NewQual.c_str()[i];
			lastQ = NewQual.c_str()[i];
		}
	}
	for (int i = 0; i<qual.size(); i++)
	{
		strand+="+"; 
	}
	Positions = NewPositions;
	ChrPositions = NewChromosome;
	//************************Lookup Kmer counts ***************************//
	LookUpKmers();
	vector <long> blank;
	
			

	//**********************************************************************//
	CountAlignmentSegments(); 	
	cout << "After getRefSeq";
	//FullOutwriteVertical(); 
}

void SamRead::LookUpKmers() 
{
	cout << "SeqSize = " << seq.size() << " RefSize = " << RefSeq.size() << endl;
	vector <long> blank;
        RefAltCounts.clear();
	RefRefCounts.clear();
	MutHashListCounts.clear();
	MutAltCounts.clear();
	MutRefCounts.clear();
	RefKmers.clear();
	AltKmers.clear();
	for (int pi = 0; pi < ParentHashes.size(); pi++){
                RefAltCounts.push_back(blank);
                RefRefCounts.push_back(blank);
        }
        for (int j = 0; j<seq.size(); j++)
        {
                //hash table already has revcomp counted in it
                string hash = getHash(seq, j, HashSize);
                string Refhash = getHash(RefSeq, j, HashSize);
		RefKmers.push_back(Refhash); 
		AltKmers.push_back(hash);
                if (hash != ""){

                        if (hash == Refhash){
                                MutAltCounts.push_back(0);
                                for (int pi = 0; pi < ParentHashes.size(); pi++){
                                        RefAltCounts[pi].push_back(0);
                                }
                        }
                        else{
                                if (MutantHashes.count(HashToLong(hash)) > 0 ){
                                        MutAltCounts.push_back(MutantHashes[HashToLong(hash)]);
                               	} 
                                else
                                        MutAltCounts.push_back(-1);
				
                                for (int pi = 0; pi < ParentHashes.size(); pi++){
                                        if (ParentHashes[pi].count(HashToLong(hash)) > 0){
                                                RefAltCounts[pi].push_back(ParentHashes[pi][HashToLong(hash)]);
                                        }
                                        else
                                                RefAltCounts[pi].push_back(-1);
                                }
                        }

			if (Hash.count(hash) > 0){
                                        MutHashListCounts.push_back(Hash[hash]);
                                }
                                else
                                        MutHashListCounts.push_back(-1);


                }
                else{
                        MutAltCounts.push_back(-3);
                        MutHashListCounts.push_back(-3); 
			for (int pi = 0; pi < ParentHashes.size(); pi++){
                                RefAltCounts[pi].push_back(-3);
                        }
                }

                if(Refhash != ""){
                        if (MutantHashes.count(HashToLong(Refhash)) > 0 )
                        {
                                MutRefCounts.push_back(MutantHashes[HashToLong(Refhash)]);
                        }
                        else
                                MutRefCounts.push_back(-1);

                        for (int pi = 0; pi < ParentHashes.size(); pi++){
                                if (ParentHashes[pi].count(HashToLong(Refhash)) > 0){
                                                RefRefCounts[pi].push_back(ParentHashes[pi][HashToLong(Refhash)]);
                                }
                                else
                                        RefRefCounts[pi].push_back(-1);
                        }

                }
                else{

                        MutRefCounts.push_back(-3);
                        for (int pi = 0; pi < ParentHashes.size(); pi++){
                                RefRefCounts[pi].push_back(-3);
                        }
                }

        }
	cout << "donezo" << endl;
}
void SamRead::parse(string read)
{
  	
	cout << "parsing " << read << endl;
	vector <string> temp = Split(read, '\t');
	cout << "boom" << endl;
	name = temp[0];
	cout << "name " << name; 
	flag = atoi(temp[1].c_str());
	chr = temp[2];
	pos = atoi(temp[3].c_str());
	mapQual = atoi(temp[4].c_str());
	cigar = temp[5];
	seq = temp[9];
	if (temp[10] == "*")
	{
		cout << "correcting missing qualtiy" << endl; 
		string newQual = ""; 
		for (int i = 0; i < seq.size(); i++)
		{
			newQual+='5'; 
		}
		temp[10] = newQual; 
	}
	qual = temp[10];
	alignments.clear(); 
	
	UsedForBigVar = false;UsedForBigVar = false;  
	first = true; 
	combined = false; 
	string NewSeq = ""; 
	for (int i = 0; i < seq.size(); i++)
		NewSeq+=toupper(seq.c_str()[i]);

	seq = NewSeq; 
	
	processCigar();
	cout << "working on strand bias" << endl; 
	cout << temp[0] << endl; 
	vector <string> temp2 = Split(name, ':');
	if (temp2.size() >= 2){
		cout << "break it down " << endl;
		strands = temp2[1];
		cout << "strands = " << strands << endl;
		forward = 0;
		reverse = 0; 
		forward = atoi(temp2[1].c_str()); 
		reverse = atoi(temp2[2].c_str()); 
		cout << "forward = " << forward << endl; 
		cout << "reverse = " << reverse << endl; 
		StrandBias = ((float)forward)/((float)forward+(float)reverse); 
		cout << "strand bias = " << StrandBias << endl; 
	}
	else{
		cout << "no strand data " << temp2.size() << endl; 
		strands = ""; 
		StrandBias = -1; 
		forward = -1; 
		reverse = -1;
	}
	AlignScore = 0;
        for (int i = 11; i< temp.size(); i++){
          vector<string> astemp = Split(temp[i], ':');
          if (astemp[0] == "AS"){
            AlignScore = atoi(astemp[2].c_str());
          }
        }
	cout << "getting flag bits " << endl;
	for (int j = 0;  j < 16;  ++j){
		FlagBits [j] =  0 != (flag & (1 << j));
	}
	cout << "Read Pared = " << FlagBits[0] << endl;
	cout << "read mapped in proper pair = " << FlagBits[1] << endl; 
	cout << "read unmapped = " << FlagBits[2] << endl; 
	cout << "mate unmapped = " << FlagBits[3] << endl; 
	cout << "read reverse strand = " << FlagBits[4] << endl; 
	cout << "mate referse strand = " << FlagBits[5] << endl; 
	cout << "first in pair =" << FlagBits[6] << endl; 
	cout << "second in pair =" << FlagBits[7] << endl; 
	cout << "not primary alignment =" << FlagBits[8] << endl; 
	cout << "read fails platform or vendor quality checks =" << FlagBits[9] << endl; 
	cout << "read is PCR or optical duplicate =" << FlagBits[10] << endl; 
	cout << "supplementary alignment =" << FlagBits[11] << endl; 
}
int findBreak(SamRead& read)
{
	char Afirst = read.cigarString.c_str()[0];


	cout << "Afirst = " << Afirst  << endl;
	cout << "starting A check " << endl;
	if (Afirst == 'H' or Afirst == 'S')
	{
		cout << "forward" << endl;
		for (int i =0; i < read.seq.size(); i++)
		{
			if  (read.cigarString.c_str()[i] == 'H' or read.cigarString.c_str()[i] == 'S')
			{
				//keep going
				 cout << i <<  " == " << read.cigarString.c_str()[i] << ' ' << read.seq.c_str()[i]<< endl;
			}
			else
			{
				cout << i <<  " == " << read.cigarString.c_str()[i]  << ' ' << read.seq.c_str()[i] << endl;
				cout << "fond break at " << i << endl;
				return i;
			}
		}
	}
	else
	{
		cout << "reverse" << endl;
		for (int i = read.seq.size()-1; i >= 0; i += -1)
		{
			if  (read.cigarString.c_str()[i] == 'H' or read.cigarString.c_str()[i] == 'S')
			{
				cout << i <<  " == " << read.cigarString.c_str()[i]  << ' ' << read.seq.c_str()[i] << endl;
				//keep going
			}
			else
			{
				cout << i <<  " == " << read.cigarString.c_str()[i]  << ' ' << read.seq.c_str()[i] << endl;
				cout << "fond break at " << i << endl;
				return i;
		       }
		}
	}


}
SamRead BetterWay(vector<SamRead> reads)
{
	cout << "In BetterWay" << endl;
	int A = 0; 
	int B = 1;
	//cout << "B = " << B << endl;
	cout << "Working on " << reads[A].name << ", size = " << reads[B].pos -reads[A].pos << endl;
	vector<vector<int>> AlignmentPos;
	vector<vector<string>> AlignmentChr;
	int ALastGoodRef = -1;
	string ALastGoodChr = "nope"; 
	int BLastGoodRef = -1; 
	string BLastGoodChr = "nope";


	vector <string> NewSeqs; 
	vector <string> NewQuals; 
	vector <string> NewRefs; 
	vector <string> NewCigars; 

	for(int i = 0; i<reads.size(); i++)
	{
		NewSeqs.push_back("");
		NewQuals.push_back(""); 
		NewRefs.push_back("");
		NewCigars.push_back("");
	}
	
	
	for (int i = 0; i < reads.size(); i++)
	{
		cout << "Pre Lining up reads " << i << endl;
		reads[i].write();
	}
	int Acount = 0; 
	int Bcount = 0;
	 if (B == 1 and GetReadOrientation(reads[A].flag) != GetReadOrientation(reads[B].flag))
	{
		cout <<"FLIPPING reads not on the same strand"; 
		reads[B].flipRead(); 
	}
	while (Acount < reads[A].cigarString.size() and Bcount < reads[B].cigarString.size())
	{
		vector<int> currentPos;
		vector<string> currentChr;
		//need to get all the reads lined up with the same number of bases, taking account of I's and D's that change length//
		if (reads[A].cigarString.c_str()[Acount] == 'D' and reads[B].cigarString.c_str()[Acount] != 'D')
		{
			currentPos.push_back(reads[A].Positions[Acount]);
			currentPos.push_back(-1);

			currentChr.push_back(reads[A].ChrPositions[Acount]);
			currentChr.push_back("nope");
		
			ALastGoodRef = reads[A].Positions[Acount];
			ALastGoodChr = reads[A].ChrPositions[Acount];
			
			NewSeqs[A]+=reads[A].seq.c_str()[Acount];
			NewQuals[A]+=reads[A].qual.c_str()[Acount];
			NewRefs[A]+=reads[A].RefSeq.c_str()[Acount];
			NewCigars[A]+=reads[A].cigarString.c_str()[Acount];

			Acount++; 

			NewSeqs[B]+='-';
			NewQuals[B]+='!';
			NewRefs[B]+="-";
			NewCigars[B]+='R';
		}
		else if (reads[A].cigarString.c_str()[Acount] != 'D' and reads[B].cigarString.c_str()[Acount] == 'D')
		{
			currentPos.push_back(-1);
			currentPos.push_back(reads[B].Positions[Bcount]);
			
			currentChr.push_back("nope");
			currentChr.push_back(reads[B].ChrPositions[Bcount]);

			BLastGoodRef = reads[B].Positions[Bcount];
			BLastGoodChr = reads[B].ChrPositions[Bcount];

			NewSeqs[B]+=reads[B].seq.c_str()[Bcount];
			NewQuals[B]+=reads[B].qual.c_str()[Bcount];
			NewRefs[B]+=reads[B].RefSeq.c_str()[Bcount];
			NewCigars[B]+=reads[B].cigarString.c_str()[Bcount];
			
			Bcount++;

			NewSeqs[A]+='-';
			NewQuals[A]+='!';
			NewRefs[A]+='-';
			NewCigars[A]+='R'; 
		}
		else
		{
			if (reads[A].cigarString.c_str()[Acount] == 'H' or reads[A].cigarString.c_str()[Acount] == 'S')
			{
				currentPos.push_back(-1); 
				currentChr.push_back("nope");

				NewSeqs[A]+=reads[A].seq.c_str()[Acount];
				NewQuals[A]+=reads[A].qual.c_str()[Acount];
				NewRefs[A]+=reads[A].RefSeq.c_str()[Acount];
				NewCigars[A]+=reads[A].cigarString.c_str()[Acount];
				
				Acount++; 
		
			}
			else if (reads[A].cigarString.c_str()[Acount] == 'M' or reads[A].cigarString.c_str()[Acount] == 'X' or reads[A].cigarString.c_str()[Acount] == 'D') 
			{
				currentPos.push_back(reads[A].Positions[Acount]);
				currentChr.push_back(reads[A].ChrPositions[Acount]);
				ALastGoodRef = reads[A].Positions[Acount]; 
				ALastGoodChr = reads[A].ChrPositions[Acount]; 
				
				NewSeqs[A]+=reads[A].seq.c_str()[Acount];
				NewQuals[A]+=reads[A].qual.c_str()[Acount];
				NewRefs[A]+=reads[A].RefSeq.c_str()[Acount];
				NewCigars[A]+=reads[A].cigarString.c_str()[Acount];
				
				Acount++; 
			}
			else if (reads[A].cigarString.c_str()[Acount] == 'I')
			{
				currentPos.push_back(ALastGoodRef); 
				currentChr.push_back(ALastGoodChr);

				NewSeqs[A]+=reads[A].seq.c_str()[Acount];
				NewQuals[A]+=reads[A].qual.c_str()[Acount];
				NewRefs[A]+=reads[A].RefSeq.c_str()[Acount];
				NewCigars[A]+=reads[A].cigarString.c_str()[Acount];
				
				Acount++; 
			}
			else 
				cout << "WTF, cigar = " << reads[A].cigarString.c_str()[Acount]; 
			
			
			
			if (reads[B].cigarString.c_str()[Bcount] == 'H' or  reads[B].cigarString.c_str()[Bcount] == 'S')
			{
				currentPos.push_back(-1);
				currentChr.push_back("nope");
				
				NewSeqs[B]+=reads[B].seq.c_str()[Bcount];
				NewQuals[B]+=reads[B].qual.c_str()[Bcount];
				NewRefs[B]+=reads[B].RefSeq.c_str()[Bcount];
				NewCigars[B]+=reads[B].cigarString.c_str()[Bcount];
			
				Bcount++; 
			}
			else if (reads[B].cigarString.c_str()[Bcount] == 'M' or reads[B].cigarString.c_str()[Bcount] == 'X' or reads[B].cigarString.c_str()[Bcount] == 'D')
			{
				currentPos.push_back(reads[B].Positions[Bcount]);
				currentChr.push_back(reads[B].ChrPositions[Bcount]);
				BLastGoodRef = reads[B].Positions[Bcount];
				BLastGoodChr = reads[B].ChrPositions[Bcount]; 
			
				NewSeqs[B]+=reads[B].seq.c_str()[Bcount];
				NewQuals[B]+=reads[B].qual.c_str()[Bcount];
				NewRefs[B]+=reads[B].RefSeq.c_str()[Bcount];
				NewCigars[B]+=reads[B].cigarString.c_str()[Bcount];
				
				Bcount++; 
			}
			else if (reads[B].cigarString.c_str()[Bcount] == 'I')
			{
				currentPos.push_back(BLastGoodRef);
				currentChr.push_back(BLastGoodChr);
	
				NewSeqs[B]+=reads[B].seq.c_str()[Bcount];
				NewQuals[B]+=reads[B].qual.c_str()[Bcount];
				NewRefs[B]+=reads[B].RefSeq.c_str()[Bcount];
				NewCigars[B]+=reads[B].cigarString.c_str()[Bcount];
	
				Bcount++; 
			}
			else
				cout << "WTF, cigar = " << reads[B].cigarString.c_str()[Bcount];  
		
		}
		AlignmentPos.push_back(currentPos);
		AlignmentChr.push_back(currentChr);
	}

	//write  
	for (int i =0; i< reads.size(); i++)
	{
		reads[i].Positions.clear(); 
		reads[i].ChrPositions.clear(); 
		reads[i].seq = NewSeqs[i];
		reads[i].qual = NewQuals[i];
		reads[i].RefSeq = NewRefs[i]; 
		reads[i].cigarString = NewCigars[i];
		for (int j = 0; j<AlignmentPos.size(); j++)
		{reads[i].Positions.push_back(AlignmentPos[j][i]);}
		for (int j = 0; j<AlignmentPos.size(); j++)
		{reads[i].ChrPositions.push_back(AlignmentChr[j][i]);}
		
	}

	bool deletion = true;
	string NewCigar = "";
	string NewSeq = "";
	string NewQual = "";
	string NewRef = "";
	vector<int>NewPos; 
	vector<string> NewChr; 
	
	char LastAlignedQ = ' ';
	int LastAlignedPos = -1; 
	string LastAlignedChr = "nope"; 


	//set LastAlignedPos to the first base with an aligned base
	bool notfound = true; 
	int base = 0; 
	while (notfound) 
	{
		for (int i = 0; i< reads.size(); i++) 
		{
			if (reads[i].Positions[base] > -1)
			{	
				LastAlignedPos = reads[i].Positions[base];
				LastAlignedChr = reads[i].ChrPositions[base];
				notfound = false;
				break; 
			}
		}
		if (notfound)
			base++;
	}
	for (int i =0; i < base; i++)
	{

		NewCigar += reads[A].cigarString.c_str()[i];
					NewSeq +=reads[A].seq.c_str()[i];
					NewQual += reads[A].qual.c_str()[i];
					NewRef+= reads[A].RefSeq.c_str()[i];
					NewPos.push_back(reads[A].Positions[i]);
					NewChr.push_back(reads[A].ChrPositions[i]);
	}
	//corect qualites so everyone has the same ones, H will produce no quality 
	cout << "checking quals" << endl;
	string bestQual = reads[0].qual; 
	for (int i =0; i<reads.size(); i++)
	{	
		bool h = false; 
		cout << "read " << reads[i].name << endl;
		for (int j = 0; j < reads[i].RefSeq.size(); j ++)
		{
			if (reads[i].RefSeq.c_str()[j] == 'H')
			{
				h = true;
			}
		}
		if (h)
		{
			cout << "contains H" << reads[i].RefSeq << endl;
		}
		else
		{
			cout << "does not contain H " << reads[i].RefSeq << endl;
			bestQual = reads[i].qual; 
		}
		 
		
						
	}
	for (int i =0; i<reads.size(); i++)
	{
		reads[i].qual = bestQual; 
	}
		
	for (int i =0; i< reads.size(); i++)
	{
		reads[i].createPeakMap();
	}

	cout << "Post adjustment" << endl;
	reads[A].write();
	reads[B].write();
	//for (int i =0; i< reads[A].PeakMap.size(); i++)
	 //     cout << i << "-" << reads[A].PeakMap[i] << "-" << reads[B].PeakMap[i] << endl;
	cout << "************INTO**********" << endl;
	if (B == 1 and GetReadOrientation(reads[A].flag) == GetReadOrientation(reads[B].flag)) //if they are on the same strand, and there are only two split reads  
	{
	  	cout << "same strand and only split reads" << endl;
		for (int i =base; i < reads[A].seq.size(); i++)
		{
			if (reads[A].Positions[i] > -1) //if this base is aligned in A 
			{
				if (reads[A].Positions[i] - LastAlignedPos > 1) //indicates a deletion
				{
					if (LastAlignedQ == '!' or reads[A].qual.c_str()[i] == '!' )
					{
						cout << "well fuck this shit A" << endl;
						return reads[A]; 
					}
					//if(reads[A].ChrPositions[i] == LastAlignedChr and abs(reads[A].Positions[i] -LastAlignedPos ) < MaxVarentSize ) //must be on the same chromosome
					if(reads[A].chr == reads[B].chr and abs(reads[A].Positions[i] -LastAlignedPos ) < MaxVarentSize ) //must be on the same chromosome
					{
						cout << "wel this dosnt make any sense" << endl; //reads are in order in the bam so A should always be downstream of B, theus the deletion shoould be detected in B
						BEDBigStuff << reads[A].chr << "\t" << LastAlignedPos << "\t" << reads[A].Positions[i] << "\t" << "Deletion" << endl;
						for (int j = LastAlignedPos; j<reads[A].Positions[i]-1; j++)
						{		
							NewCigar += 'D';
							NewSeq += '-';
							NewQual += LastAlignedQ;
							NewRef+= toupper(Reff.getSubSequence(reads[A].chr, i, 1).c_str()[0]);
							NewPos.push_back(j);
							NewChr.push_back(reads[A].ChrPositions[i]);
						}
					}
					else
					{
						//if( reads[A].ChrPositions[i] == LastAlignedChr and abs(reads[A].Positions[i] -LastAlignedPos ) >= MaxVarentSize )
						if( reads[A].chr == reads[B].chr and abs(reads[A].Positions[i] -LastAlignedPos ) >= MaxVarentSize )
						{
							int Abreak = findBreak(reads[A]);
							int Bbreak = findBreak(reads[B]);
							cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
							if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
							{
								cout <<  "INVERSION written to file"  << endl;
								Translocations << "Too Big, Same strand and chr "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
								reads[A].writetofile(Translocations);
								reads[B].writetofile(Translocations);
								Translocations << endl << endl;
								Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl; 
						       }
							else
							{
								cout <<  "INVERSION skipped"  << endl;
							}
						//	if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
						//	{
						//		Translocations << "TOO BIG " << abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
						//		reads[A].writetofile(Translocations);
						//		reads[B].writetofile(Translocations);
						//		Translocations << endl << endl;	
						//	}
						}
						else
						{
							if ((reads[A].chr == "hs37d5" and reads[B].chr != "hs37d5" ) or  (reads[A].chr != "hs37d5" and reads[B].chr == "hs37d5" ))
							{
								int Abreak = findBreak(reads[A]);
								int Bbreak = findBreak(reads[B]);
								cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
								if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0 )
								{
									cout <<  "INVERSION written to file"  << endl;
									Translocations << "Possible mob event "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
									reads[A].writetofile(Translocations);
									reads[B].writetofile(Translocations);
									Translocations << endl << endl;
									Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl;
								}
								else
								{
									cout <<  "INVERSION skipped"  << endl;
								}
						//		if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
						//		{
						//			Translocations << "mobil elemnt " << endl;
						//			reads[A].writetofile(Translocations); 
						//			reads[B].writetofile(Translocations); 
						//			Translocations << endl << endl;
						//		}
							}
							else
							{
							int Abreak = findBreak(reads[A]);
							int Bbreak = findBreak(reads[B]);
							cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
							if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
							{
								cout <<  "INVERSION written to file"  << endl;
								Translocations << "Translocataion, same strand "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
								reads[A].writetofile(Translocations);
								reads[B].writetofile(Translocations);
								Translocations << endl << endl;
								Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl; 
						       }
							else
							{
								cout <<  "INVERSION skipped"  << endl;
							}	
							//	if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
							//	{
							//		Translocations << "we got a translocation" << endl;
							//		reads[A].writetofile(Translocations);
							//		reads[B].writetofile(Translocations);
							//		Translocations << endl << endl;
							//	}
							}	
						}
					}
				}
				else if (reads[A].Positions[i] -LastAlignedPos < 0 and abs(reads[A].Positions[i] -LastAlignedPos ) < MaxVarentSize) // indicates a possible insertion or tandem duplication
				{
					 if (LastAlignedQ == '!'  or reads[A].qual.c_str()[i] == '!')
					{
						cout << "well fuck this shit B" << endl;
						return reads[A];
					}
					cout << "this could be one A, last = " << LastAlignedPos << " Current = " << reads[A].Positions[i] << " chr = " << LastAlignedChr << " and " << reads[A].ChrPositions[i] << endl;
//					if (reads[A].ChrPositions[i] == LastAlignedChr )
					if (reads[A].chr == reads[B].chr )
					{
						cout << "This is an insertion A at base " << i <<  endl;
						BEDBigStuff << reads[A].chr << "\t" << reads[A].Positions[i] << "\t" << LastAlignedPos << "\tTandemDup" << endl;
						cout << "tadem dup" << endl;
						int j = 0; 
						for( j = i; j<reads[A].seq.size() and reads[A].Positions[j] < LastAlignedPos ; j++)
						{
							NewCigar += 'Y'; //'I';
							NewSeq += reads[A].seq.c_str()[j];
							NewQual += reads[A].qual.c_str()[j];
							NewRef+= '-';
							NewPos.push_back(reads[A].Positions[i]);
							NewChr.push_back(reads[A].ChrPositions[i]);
						}
						i=j;
						//find the last base that was aligned, htere can be novel insertion stuff so you cant jsut take the last base
						int k; 
						for (k =reads[A].Positions.size()-1; k >=0; k+= -1)
						{
							if (reads[A].Positions[k]+1 > 1)
								break;
						}
						//cout << "j= " << reads[A].Positions[k]+1 << " < " <<  LastAlignedPos << " - " << endl;
						for( j = reads[A].Positions[k]+1; j < LastAlignedPos; j++)
						{
							NewCigar += 'Y'; //'I';
							NewSeq += toupper(Reff.getSubSequence(reads[A].chr, j, 1).c_str()[0]);
							NewQual += '!';
							NewRef+= '-';
							NewPos.push_back(j);
							NewChr.push_back(reads[A].chr);
						}
						cout << "yaya finished" << endl;
						 
					}
					else
					{
						if ((reads[A].chr == "hs37d5" and reads[B].chr != "hs37d5" ) or  (reads[A].chr != "hs37d5" and reads[B].chr == "hs37d5" ))
						{
							int Abreak = findBreak(reads[A]);
							int Bbreak = findBreak(reads[B]);
							cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
							if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
							{
								cout <<  "INVERSION written to file"  << endl;
								Translocations << "Possible mob event "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
								reads[A].writetofile(Translocations);
								reads[B].writetofile(Translocations);
								Translocations << endl << endl;
								Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl;
							}
							else
							{
								cout <<  "INVERSION skipped"  << endl;
							}
							//if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
							//{
							//	Translocations << "mobil elemnt " << endl;
							//	reads[A].writetofile(Translocations);
							//	reads[B].writetofile(Translocations);
							//	Translocations << endl << endl;
							//}
						}
				 		int Abreak = findBreak(reads[A]);
						int Bbreak = findBreak(reads[B]);
					       	cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
						if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
						{
							cout <<  "INVERSION written to file"  << endl;
							Translocations << "Translocation, same strand "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
							reads[A].writetofile(Translocations);
							reads[B].writetofile(Translocations);
							Translocations << endl << endl;
				 			Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl;
					       	}
						else
						{
							cout <<  "INVERSION skipped"  << endl;
						}
				//		if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
				//		{
				//			Translocations << "we got a translocation" << endl;
				//			reads[A].writetofile(Translocations);
				//			reads[B].writetofile(Translocations);
				//			Translocations << endl << endl;
				//		}
					}
	
				}
				else if( reads[A].Positions[i] -LastAlignedPos < 0  and abs(reads[A].Positions[i] -LastAlignedPos ) >= MaxVarentSize )
				{
					 if (LastAlignedQ == '!'  or reads[A].qual.c_str()[i] == '!')
					{
						cout << "well fuck this shit C" << endl;

						return reads[A];
					}
					int Abreak = findBreak(reads[A]);
					int Bbreak = findBreak(reads[B]);
					cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
					if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
					{
						cout <<  "INVERSION written to file"  << endl;
						if (reads[A].chr == reads[B].chr)
							Translocations << "TOO BIG 3 "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
						else 
							Translocations << "Translocation 3 "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
						reads[A].writetofile(Translocations);
						reads[B].writetofile(Translocations);
						Translocations << endl << endl;
						Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl; 
				       }
					else
					{
						cout <<  "INVERSION skipped"  << endl;
					}
			//		if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
			//		{
			//			Translocations << "TOO BIG " << abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
			//			reads[A].writetofile(Translocations);
			//			reads[B].writetofile(Translocations);
			//			Translocations << endl << endl;
			//		}
				}

				
				if (i<reads[A].cigarString.size())
				{
					NewCigar += reads[A].cigarString.c_str()[i];
					NewSeq += reads[A].seq.c_str()[i];
					NewQual += reads[A].qual.c_str()[i];
					NewRef+= reads[A].RefSeq.c_str()[i];
					NewPos.push_back(reads[A].Positions[i]);
					NewChr.push_back(reads[A].ChrPositions[i]);
		       	 		LastAlignedQ = reads[A].qual.c_str()[i];
					LastAlignedPos = reads[A].Positions[i];
					LastAlignedChr = reads[A].ChrPositions[i];
					//cout <<	 "A " << reads[A].seq.c_str()[i] << " " << reads[A].cigarString.c_str()[i] << " " << reads[A].RefSeq.c_str()[i] <<" " << reads[A].ChrPositions[i] << " " << reads[A].Positions[i] << NewSeq << endl;
				}	
			}
			else if (reads[B].Positions[i] > -1)
			{
				if (reads[B].Positions[i] - LastAlignedPos > 1)
				{
					 if (LastAlignedQ == '!'  or reads[B].qual.c_str()[i] == '!')
					{
						
						cout << "well fuck this shit D" << endl;

						return reads[A];
					}
					//if(reads[B].ChrPositions[i] == LastAlignedChr  and abs(reads[B].Positions[i] - LastAlignedPos ) < MaxVarentSize )
					if(reads[B].chr == reads[A].chr  and abs(reads[B].Positions[i] - LastAlignedPos ) < MaxVarentSize )
					{
						cout << "striahgtup deletion, size = " <<  abs(reads[B].Positions[i] -LastAlignedPos ) << " at base " << i << " from Position " << LastAlignedPos << " to " << reads[B].Positions[i] << endl; 
						BEDBigStuff << reads[B].chr << "\t" << LastAlignedPos << "\t" << reads[B].Positions[i]  << "\t" << "Deletion" << endl;
//						cout << "Inserting reff sequence from " << LastAlignedPos+1 << " to " << reads[B].Positions[i] << " = " << reads[B].Positions[i]-LastAlignedPos << endl;
						for (int j = LastAlignedPos; j< reads[B].Positions[i]-1; j++)
			       	 		{
							//cout << j << " - " << j - LastAlignedPos<< endl;
			       				NewCigar += 'D';
			       				NewSeq += '-';
			       		 		NewQual += LastAlignedQ;
			       		 		NewRef+= toupper(Reff.getSubSequence(reads[B].chr, j, 1).c_str()[0]);
							NewPos.push_back(j);
							NewChr.push_back(reads[B].ChrPositions[i]);

	//						char tmp = toupper(Reff.getSubSequence(reads[B].chr, j, 1).c_str()[0]); 
	  //		      			cout << "R " << '-' << " " << 'D' << " " << tmp << " " << reads[B].chr << " " << j << " " << NewSeq << endl;
	//						cout << "R " << '-' << " " << 'D' << " " << tmp << " " << reads[B].chr << " " << j << " " << NewRef << endl << endl;
						}
					}
					else
					{
						// if( reads[A].ChrPositions[i] == LastAlignedChr and abs(reads[B].Positions[i] -LastAlignedPos ) >= MaxVarentSize )
						if( reads[A].chr == reads[B].chr and abs(reads[B].Positions[i] -LastAlignedPos ) >= MaxVarentSize )
						{
							int Abreak = findBreak(reads[A]);
							int Bbreak = findBreak(reads[B]);
							cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
							if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
							{
								cout <<  "INVERSION written to file"  << endl;
								Translocations << "TOO BIG 2 "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
								reads[A].writetofile(Translocations);
								reads[B].writetofile(Translocations);
								Translocations << endl << endl;
								Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl;
							}
							else
							{
								cout <<  "INVERSION skipped"  << endl;
							}
						//	if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
						//	{
						//		Translocations << "TOO BIG " << abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
						//		reads[A].writetofile(Translocations);
						//		reads[B].writetofile(Translocations);
						//		Translocations << endl << endl;
						//	}
						}
						else
						{
							if ((reads[A].chr == "hs37d5" and reads[B].chr != "hs37d5" ) or  (reads[A].chr != "hs37d5" and reads[B].chr == "hs37d5" ))
							{
								int Abreak = findBreak(reads[A]);
								int Bbreak = findBreak(reads[B]);
								cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
								if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
								{
									cout <<  "INVERSION written to file"  << endl;
									Translocations << "Possible mob event "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
									reads[A].writetofile(Translocations);
									reads[B].writetofile(Translocations);
									Translocations << endl << endl;
									Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl;
								}
								else
								{
									cout <<  "INVERSION skipped"  << endl;
								}
							//	if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
							//	{
							//		Translocations << "mobil elemnt " << endl;
							//		reads[A].writetofile(Translocations);
							//		reads[B].writetofile(Translocations);
							//		Translocations << endl << endl;
							//	}
							}
							else	
							{
								int Abreak = findBreak(reads[A]);
								int Bbreak = findBreak(reads[B]);
								cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
								if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
								{
									cout <<  "INVERSION written to file"  << endl;
									Translocations << "Translocation 2 "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
									reads[A].writetofile(Translocations);
									reads[B].writetofile(Translocations);
									Translocations << endl << endl;
									Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl; 
							       }
								else
								{
									cout <<  "INVERSION skipped"  << endl;
								}

							
							//	if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
							//	{
							//		Translocations << "we got a translocation" << endl;
							//		reads[A].writetofile(Translocations);
							//		reads[B].writetofile(Translocations);
							//		Translocations << endl << endl;
							//	}
							}
						}
					}

				}
				else if (reads[B].Positions[i] -LastAlignedPos < 0 and abs(reads[B].Positions[i] - LastAlignedPos ) < MaxVarentSize) // indicates a possible insertion or tandem duplication
				{
					 if (LastAlignedQ == '!' or reads[B].qual.c_str()[i] == '!' )
					{
						cout << "well fuck this shit E" << endl;

						return reads[A];
					}
					 // cout << "this could be one B, last = " << LastAlignedPos << " Current = " << reads[B].Positions[i] << " chr = " << LastAlignedChr << " and " << reads[B].ChrPositions[i] << endl;
					//if (reads[B].ChrPositions[i] == LastAlignedChr  
					if (reads[B].chr == reads[A].chr  )
					{
						cout << "This is an insertion B at base " << i << endl;

						BEDBigStuff << reads[B].chr << "\t" << reads[B].Positions[i] << "\t" << LastAlignedPos << "\tTandemDup" << endl;
					//	cout << "tadem dup" << endl;
						int j = 0; 
						for(j = i; j < reads[B].seq.size() and reads[B].Positions[j] <= LastAlignedPos; j++)
						{
							NewCigar += 'Y';
							NewSeq += reads[B].seq.c_str()[j];
							NewQual += reads[B].qual.c_str()[j];
							NewRef+= '-';
							NewPos.push_back(reads[B].Positions[i]);
							NewChr.push_back(reads[B].ChrPositions[i]);
						}
						i=j;
						//need to find last base that was alinged, insetion can mess this up so you can just take the last base pos 
						int k;
						for (k =reads[B].Positions.size()-1; k >=0; k+= -1)
						{
							if (reads[B].Positions[k]+1 > 1)
								break;
						}
						//cout << "j= " << reads[B].Positions[k]+1 << " < " <<  LastAlignedPos << " - " << endl;
						for( j = reads[B].Positions[k]+1; j < LastAlignedPos; j++) 
						{
							NewCigar += 'Y';
							NewSeq += toupper(Reff.getSubSequence(reads[B].chr, j, 1).c_str()[0]);
							NewQual += '!';
							NewRef+= '-';
							NewPos.push_back(j);
							NewChr.push_back(reads[B].chr);
						}
						
					}
					else
					{
						cout << "we got a translocation" << endl;
						if ((reads[A].chr == "hs37d5" and reads[B].chr != "hs37d5" ) or  (reads[A].chr != "hs37d5" and reads[B].chr == "hs37d5" ))
						{
							cout << "mobil elemnt " << endl;
							//reads[A].write();
							//reads[B].write();
						}
					}

				}
				else if( reads[B].Positions[i] -LastAlignedPos < 0  and abs(reads[B].Positions[i] -LastAlignedPos ) >= MaxVarentSize )
				{
					 if (LastAlignedQ == '!' or reads[B].qual.c_str()[i] == '!')
					{
						cout << "well fuck this shit F" << endl;
						return reads[A];
					}
				 	int Abreak = findBreak(reads[A]);
					int Bbreak = findBreak(reads[B]);



		       		 	cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
					if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
					{
						cout <<  "INVERSION written to file"  << endl;
						Translocations << "TOO BIG 1 "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
						reads[A].writetofile(Translocations);
						reads[B].writetofile(Translocations);
						Translocations << endl << endl;
	 					Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl;
			       		}
					else
					{
				 		cout <<  "INVERSION skipped"  << endl;
					}

				//	if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
				//       	{
				//		Translocations << "TOO BIG " << abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
				//		reads[A].writetofile(Translocations);
				//		reads[B].writetofile(Translocations);
				//		Translocations << endl << endl;
				//	}
				}


				if(i<reads[B].cigarString.size())
				{
					NewCigar += 'M'; //reads[B].cigarString.c_str()[i];
					NewSeq += reads[B].seq.c_str()[i];
					NewQual	+= reads[B].qual.c_str()[i];
					NewRef+= reads[B].RefSeq.c_str()[i];
					NewPos.push_back(reads[B].Positions[i]);
					NewChr.push_back(reads[B].ChrPositions[i]);

					LastAlignedQ = reads[B].qual.c_str()[i];
					LastAlignedPos = reads[B].Positions[i];
					LastAlignedChr = reads[B].ChrPositions[i];
					//	 cout << "B " << reads[B].seq.c_str()[i] << " " << reads[B].cigarString.c_str()[i] << " " << reads[B].RefSeq.c_str()[i] <<" " << reads[B].ChrPositions[i] << " " << reads[B].Positions[i] << NewSeq << endl;
				}
			}
			else
			{
				if (reads[A].cigarString.c_str()[i] == 'S')
				{
					NewCigar += reads[A].cigarString.c_str()[i];
					NewSeq +=reads[A].seq.c_str()[i];
					NewQual += reads[A].qual.c_str()[i];
					NewRef+= reads[A].RefSeq.c_str()[i];
					NewPos.push_back(reads[A].Positions[i]);
					NewChr.push_back(reads[A].ChrPositions[i]);

				}
				else if (reads[B].cigarString.c_str()[i] == 'S')
				{
					NewCigar += reads[B].cigarString.c_str()[i];
					NewSeq +=reads[B].seq.c_str()[i];
					NewQual += reads[B].qual.c_str()[i];
					NewRef+= reads[B].RefSeq.c_str()[i];
					NewPos.push_back(reads[B].Positions[i]);
					NewChr.push_back(reads[B].ChrPositions[i]);

				}
				else
					cout << "no base works :(" << endl;
				
			}
		}

		//fix inernal S bases
		int First = -1; 
		int Last = -1; 
		string NewNewCigar = ""; 
		for (int i =0; i < NewCigar.size(); i++)
		{
			if(NewCigar.c_str()[i] != 'S' and NewCigar.c_str()[i] != 'H')
			{
				First = i; 
				break; 
			}
		}
		for (int i = NewCigar.size() -1; i>=0; i--)
		{
			if(NewCigar.c_str()[i] != 'S' and NewCigar.c_str()[i] != 'H')
			{
				Last = i; 
				break; 
			}
		}

		for (int i = 0; i<NewCigar.size(); i++)
		{
			if (i > First and i < Last)
			{
				if (NewCigar.c_str()[i] == 'S' or NewCigar.c_str()[i] == 'H')
					NewNewCigar+= 'I';
				else 
					NewNewCigar+=NewCigar.c_str()[i]; 
			}
			else 
				NewNewCigar+=NewCigar.c_str()[i];
					
		}
		NewCigar = NewNewCigar;
	
		int UnalignedCount = 0;
		for (int i =0; i<NewCigar.size(); i++)
		{
			if (NewCigar.c_str()[i] == 'H' || NewCigar.c_str()[i] == 'S')
				UnalignedCount++;
		}
		cout << "Unaligned Baises = " << UnalignedCount << " %= " << (double)UnalignedCount/(double)NewCigar.size() << endl;
		if (UnalignedCount < 150)
		{
			NewCigar = NewNewCigar; 
			//cout << NewSeq << endl << NewQual << endl << NewCigar << endl << NewRef << endl;
			reads[A].first = true;
//			cout <<"redoing alignemnt number" << endl;
			int temp = reads[A].alignments[0];
			reads[A].alignments.clear();
			reads[A].alignments.push_back(temp);
//			cout << "alignemtns should be " << reads[A].alignments.size() << endl;
	
	       		reads[A].cigarString = NewCigar;
			reads[A].seq = NewSeq;
			reads[A].qual = NewQual;
			reads[A].RefSeq = NewRef;
			reads[A].Positions.clear();
			reads[A].ChrPositions.clear();
			reads[A].Positions = NewPos; 
			reads[A].ChrPositions = NewChr; 
			reads[A].combined = true;
			reads[B].combined = true;
		
		//reads[A].write(); 		
		//{
		//reads[A].write(); 
	       	//reads[A].writeVertical();
		//}		
	
		}
		else
		{
			cout << "Skipping" << endl;
		}
		cout <<	 "Done with " << reads[A].name << endl;
	}
	else //if (3==1)
	{
	  	cout << "here somehow" << endl; 
		//this needs to be the invertion stuff 
		//need to add here looing at every base 
		if( reads[A].chr == reads[B].chr )//and abs(reads[A].Positions[i] -LastAlignedPos ) <= MaxVarentSize )
		{
			//need to check each base and find the breakpoint 
			cout << "found possible Inversion " << endl;
			


			int Abreak = findBreak(reads[A]); 
			int Bbreak = findBreak(reads[B]);
			
			

			cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
			//if( /*1==1 or*/ reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1 or reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0
			if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
			{	
				cout <<  "INVERSION written to file"  << endl;
				Translocations << "INVERSION"  << endl;
				reads[A].writetofile(Translocations);
				reads[B].writetofile(Translocations);
				Translocations << endl << endl;
				Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl;
			}
			else
			{
				 cout <<  "INVERSION skipped"  << endl;
			}
			string Acig = "";
			string Bcig = "";	
			for (int i = 0; i < reads[A].seq.size(); i++)
			{
				char Ab = reads[A].cigarString.c_str()[i]; 
				char Bb = reads[B].cigarString.c_str()[i]; 		
				if ((reads[A].cigarString.c_str()[i] == 'M' or reads[A].cigarString.c_str()[i] == 'X') and (reads[B].cigarString.c_str()[i] == 'S' or reads[B].cigarString.c_str()[i] == 'H'))
					Bb = 'U';
				if ((reads[B].cigarString.c_str()[i] == 'M' or reads[B].cigarString.c_str()[i] == 'X') and (reads[A].cigarString.c_str()[i] == 'S' or reads[A].cigarString.c_str()[i] == 'H'))
			      		Ab = 'U'; 
				Acig += Ab; 
				Bcig += Bb; 
				   	
			}
			reads[A].cigarString = Acig; 
			reads[B].cigarString = Bcig; 
			cout << "invertion adjust string"; 
			reads[A].write(); 
			reads[B].write(); 	
		}
		else if ((reads[A].chr == "hs37d5" and reads[B].chr != "hs37d5" ) or  (reads[A].chr != "hs37d5" and reads[B].chr == "hs37d5" ))
		{
			 int Abreak = findBreak(reads[A]);
			int Bbreak = findBreak(reads[B]);



			cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
			//if( /*1==1 or*/ reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1 or reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0
			 if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
			{
				Translocations << "mobil elemnt inverted" << endl;
				reads[A].writetofile(Translocations);
				reads[B].writetofile(Translocations);
				Translocations << endl << endl;
				Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl;	
			}
		}
		else
		{
			 int Abreak = findBreak(reads[A]);
			int Bbreak = findBreak(reads[B]);



			cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
			//if( /*1==1 or*/ reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1 or reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0
			 if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
			{
				Translocations << "we got a translocation and invertion" << endl;
				reads[A].writetofile(Translocations);
				reads[B].writetofile(Translocations);
				Translocations << endl << endl;
				Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl;
			}	
		}

			
		cout << "*********SKIPPING************\ndifference strands" << endl;
		BEDNotHandled << "Different strands" << endl;
		BEDNotHandled << reads[A].chr << "\t" << reads[A].pos << "\t" << reads[A].pos+reads[A].seq.size() << "\t" << reads[A].name << "\t" << reads[A].cigar << endl;
		reads[A].writetofile(BEDNotHandled);
		BEDNotHandled << reads[B].chr << "\t" << reads[B].pos << "\t" << reads[B].pos+reads[B].seq.size() << "\t" << reads[B].name << "\t" << reads[B].cigar << endl;	
		reads[B].writetofile(BEDNotHandled);
		
		BEDNotHandled << endl << endl;


		Invertions << reads[A].chr << "\t" << reads[A].pos << "\t" << reads[B].pos << "\t" << reads[B].pos- reads[A].pos << endl;

	}

	//**********************Adding K-mer lookup stuff here *********************************
	
	//**************************************************************************************
	//reads[A].FixTandemRef(); 
	reads[A].LookUpKmers();
	cout << "ReAdjustedKmers" <<endl;
	reads[A].write(); 
	//FullOutreads[A].writeVertical(); 
	return reads[A];
	cout << "Out Of SamReadBetterWay " << endl;
}
int SamRead::CheckBasesAligned()
{
  int longest = 0;
  int count = 0; 
 for ( int j = 0; j< cigarString.size(); j++)
 {
	if (cigarString.c_str()[j] != 'H' and cigarString.c_str()[j] != 'S')
	{count++; }
	else
	{
	  if (count > longest){longest = count;}
	  count = 0; 
	}
   
 }
 if (count > longest){longest = count;}
 return longest; 
}
bool SamRead::CheckEndsAlign()
{
  int StartAlign = 0;
  int j;
  for ( j = 10; j< cigarString.size(); j++)
  {
    if (cigarString.c_str()[j] != 'H' and cigarString.c_str()[j] != 'S')
    {StartAlign++;}
    else
    {break;}
  }
  int EndAlign = 0;
  int i;
  for ( i = cigarString.size()-10; i>=0; i--)
  {
    if (cigarString.c_str()[i] != 'H' and  cigarString.c_str()[i] != 'S')
    {EndAlign++;}
    else
    {break;}
  }

  if (StartAlign > 20 or EndAlign > 20)
  {
    return true; 
  }

return false; 
}


bool SamRead::StartsWithAlignAtPeak(int &pos, string &insert, int &Kdepth)
{
  cout << "starst with " << endl;
  ///////////////////
  int EndClip = 0;
  int i;
  int PeakBases = 0; 
  for ( i = cigarString.size()-1; i>=0; i--)
  {
    if (cigarString.c_str()[i] == 'H' or cigarString.c_str()[i] == 'S')
    {EndClip++;if (PeakMap[i]==1){PeakBases++;}}
    else
    {break;}
  }
  //get the inserted sequence 
  for (int s = i; s<cigarString.size(); s++)
  {
    insert = insert+seq[s];
  }
  /////////////////////
  int StartAlign = 0;
  int j;
  for ( j = 0; j< cigarString.size(); j++)
  {
    if (cigarString.c_str()[j] != 'H' and cigarString.c_str()[j] != 'S')
    {StartAlign++;}
    else
    {break;}
  }
  /////////////////////
  cout << "pos = " << j << " - " << Positions[j] << endl;
  pos = Positions[j];
  Kdepth = (int) qual.c_str()[i]-33; 
  if (EndClip > 40 && StartAlign > 40)
  {
    if ((PeakMap[i-1] or PeakMap[i] or PeakMap[i+1]) and (PeakMap[j-1] or PeakMap[j] or PeakMap[j+1])or PeakBases > 10)
    {return true;}
  }


  return false;
 }  
bool SamRead::StartsWithAlign(int &pos, string &insert)
{
  cout << "starst with " << endl; 
  ///////////////////
  int EndClip = 0; 
  int i; 
  for ( i = cigarString.size()-1; i>=0; i--)
  {
    if (cigarString.c_str()[i] == 'H' or cigarString.c_str()[i] == 'S')
    {EndClip++;}
    else
    {break;}
  }
  //get the inserted sequence 
  for (int s = i; s<cigarString.size(); s++)
  {
    insert = insert+seq[s]; 
  }
  /////////////////////
  int StartAlign = 0; 
  int j; 
  for ( j = 0; j< cigarString.size(); j++)
  {
    if (cigarString.c_str()[j] != 'H' and cigarString.c_str()[j] != 'S')
    {StartAlign++;}
    else
    {break;}
  } 
  /////////////////////
  cout << "pos = " << j << " - " << Positions[j] << endl;
  pos = Positions[j]; 
  if (EndClip > 40 && StartAlign > 40)
  {
    
    {return true;}
  }


  return false; 
 }
bool SamRead::EndsWithAlignAtPeak(int &pos, string &insert, int &Kdepth)
{
  cout << "Running EndsWithAlign" << endl;
  ////////////////////
  int EndAlign = 0;
  int i;
  for ( i = cigarString.size()-1; i>=0; i--)
  {
    if (cigarString.c_str()[i] != 'H' and  cigarString.c_str()[i] != 'S')
    {EndAlign++;}
    else
    {break;}
  }
  pos = Positions[i+1]; 
  cout << "EndAlign = " << EndAlign << endl;
  ///////////////////
  int StartClip = 0;
  int j;
  int PeakBases = 0; 
  for ( j = 0; j< cigarString.size(); j++)
  {
    if (cigarString.c_str()[j] == 'H' or cigarString.c_str()[j] == 'S')
    {StartClip++; if (PeakMap[j]==1){PeakBases++;}}
    else
    {break;}
  }
  cout << "StartClip = " << StartClip << endl;
  
  for (int s = 0; s<j; s++)
  {insert = insert+seq[s];}
  Kdepth = (int) qual.c_str()[i]-33;
  if (EndAlign > 40 && StartClip > 40)
  {
    if ((PeakMap[i-1] or PeakMap[i] or PeakMap[i+1]) and (PeakMap[j-1] or PeakMap[j] or PeakMap[j+1])or PeakBases > 10 )
    {return true;}
  }


  return false;
 }
bool SamRead::EndsWithAlign(int &pos, string &insert)
{
  cout << "Running EndsWithAlign" << endl;
  ////////////////////
  int EndAlign = 0;
  int i;
  for ( i = cigarString.size()-1; i>=0; i--)
  {
    if (cigarString.c_str()[i] != 'H' and  cigarString.c_str()[i] != 'S')
    {EndAlign++;}
    else
    {break;}
  }
   pos = Positions[i+1]; 
  cout << "EndAlign = " << EndAlign << endl;
  ///////////////////
  int StartClip = 0;
  int j;
  for ( j = 0; j< cigarString.size(); j++)
  {
    if (cigarString.c_str()[j] == 'H' or cigarString.c_str()[j] == 'S')
    {StartClip++; }
    else
    {break;}
  }
  cout << "StartClip = " << StartClip << endl; 
  for (int s = 0; s<j; s++)
  {insert = insert+seq[s];}

  if (EndAlign > 40 && StartClip > 40)
  {
    
    {return true;}
  }


  return false;
 }
int main (int argc, char *argv[])
{

cout << "###########################RUNNING THIS ONE#########################" << endl; 
cout << "Modes	chr	pos	type	reff	alt	MutRef	MutAlt	Par1Ref	Par2Ref" << endl;	
//	ifstream testthis [100]; 
//	testthis[0].open("./test.txt"); 
//	string boom2; 
//	while (getline(testthis[0], boom2))
//	{
//		cout << boom2 << endl;
//	}
//	return 0;  
	//************************************************
	//my arg parser

	string helptext;       
	helptext = \
"\
RUFUS.interpret: converts RUFUS aligned contigs into a VCF \n\
By Andrew Farrell\n\
   The Marth Lab\n\
\n\ 
options:\
  -h [ --help ]  Print help message\n\
  -sam  argPath to input SAM file, omit for stdin\n\
  -r    argPath to reference file \n\
  -hf   argPath to HashFile from RUFUS.build\n\
  -hS   argHash Size\n\
  -o    argOutput stub\n\
  -m    argMaximum varient size: default 1Mb\n\
(Sorry it has to be a num, no 1kb, must be 1000\n\
  -c    argPath to sorted.tab file for the parent sample\n\
  -s    arg Path to sorted.tab file for the subject sample\n\
  -cR   argPath to the sorted.tab file fo the parnt sample hashes in the reference\n\
  -sR   argPath to the sorted.tab file fo the subject sample hashes in the reference\n\
  -mQ   argMinimum map quality to consider varients in\n\
  -mod argPath to the model file from RUFUS.model\n\
  -e    arg             Path to Kmer file to exlude from LowCov check\n\
";
	
	string MutHashFilePath = "" ;
	string MutHashFilePathReference = "";
	MaxVarentSize = 1000000;
	string RefFile = ""; 
	string HashListFile = "" ; 	
	string samFile = "stdin"; 
	string outStub= "";
	string ModelFilePath = "";
	string ExcludeFilePath = ""; 
	int MinMapQual = 0; 
	for(int i = 1; i< argc; i++)
	{
		cout << i << " = " << argv[i] << endl; 
	}
	cout <<"****************************************************************************************" << endl;
	vector <int> ParentHashFilePaths; 
	vector <int> ParentHashFilePathsReference;
	for(int i = 1; i< argc; i++)
	{
		string p = argv[i];
		cout << i << " = " << argv[i]<< endl;
		if( p == "-h")
		{
			//print help 
			cout << helptext << endl;
			return 0; 
		}
		else if (p == "-r")
		{
			RefFile = argv[i+1];
			i=i+1;
			 cout << "YAAAY added RefFile = " << RefFile << endl;	
		}
		else if (p == "-sam")
		{
			samFile =  argv[i+1];
			i++;
		}
		else if (p == "-o")
		{
			outStub =  argv[i+1];
			i++;
		}
		else if (p == "-hf")
		{
			HashListFile =  argv[i+1];
			i++;
		}
		else if (p == "-hs")
                {
                        HashSize =  atoi(argv[i+1]);
                        i++;
                }
		else if (p == "-m")
		{
			MaxVarentSize =  atoi(argv[i+1]);
			i++;
			cout << "YAAAY added MaxVarSize = " << MaxVarentSize << endl;
		}
		else if (p == "-c")
		{
			cout << "Par Hash = " << argv[i+1] << endl;
			ParentHashFilePaths.push_back(i+1);
			i=i+1;
		}
		else if (p == "-cR")
		{
			cout << "Par Ref Hash = " << argv[i+1] << endl;
			ParentHashFilePathsReference.push_back(i+1);
			i=i+1;
		}
		else if (p == "-s")
		{
			cout << "Sub Hash = " << argv[i+1] << endl;
			MutHashFilePath = argv[i+1];
			i+=1;
		}
		else if (p == "-sR")
		{
			cout << "Sub Hash = " << argv[i+1] << endl;
			MutHashFilePathReference = argv[i+1];
			i+=1;
		}
		else if (p == "-mod")
                {
                        cout << "model file = " << argv[i+1] << endl;
                        ModelFilePath = argv[i+1];
                        i+=1;
                }
		else if (p == "-mQ")
		{
			cout << "Min Mapping Qualtiy = " << argv[i+1] << endl;
			MinMapQual = atoi(argv[i+1]);
			i+=1;
		}
		else if(p == "-e")
		{
		  	cout << "Exclue File Path = " << argv[i+1] << endl; 
			ExcludeFilePath = argv[i+1]; 
			i+=1; 
		}
		else
		{
			cout << "ERROR: unkown command line paramater -" <<  argv[i] << "-"<< endl;
			return 0; 
		}
		
	}
	//check values 
	if (RefFile == "")
	{
		cout << "ERROR Reference required" << endl;
		return 0; 
	}
	if (HashListFile == "")
	{
		cout << "Error HashList required" << endl;
		return 0; 
	}
	if (outStub == "")
	{
		if (samFile != "stdin")
			outStub = samFile; 
		else
		{
			cout << "ERROR out file stub required " << endl;
			return -1; 
		}
	}

	for (int i = 0; i < ParentHashFilePaths.size(); i++)
	{
		ifstream reader; 
		reader.open (argv[ParentHashFilePaths[i]]); 
		string line = "";
		unordered_map <unsigned long int, int> hl;  
		while (getline(reader, line))
		{
			vector <string> temp = Split(line, ' '); 
			unsigned long hash = HashToLong(temp[0]); 
			hl[hash] = atoi(temp[1].c_str());
			hash = HashToLong(RevComp(temp[0])); 
			hl[hash] = atoi(temp[1].c_str()); 
		}
		ParentHashes.push_back(hl); 
		reader.close(); 
	}
	for (int i = 0; i < ParentHashFilePathsReference.size(); i++)
	{
		ifstream reader;
		reader.open (argv[ParentHashFilePathsReference[i]]);
		string line = "";
		unordered_map <unsigned long int, int> hl;
		while (getline(reader, line))
		{
			vector <string> temp = Split(line, ' ');
			unsigned long hash = HashToLong(temp[0]);
			ParentHashes[i][hash] = atoi(temp[1].c_str());
			hash = HashToLong(RevComp(temp[0]));
			ParentHashes[i][hash] = atoi(temp[1].c_str());
		}
		reader.close();
	}
	cout << "check parent thing" << endl; 
	for(int i =0; i < ParentHashes.size(); i++)
	{
		cout << "sample " << i << endl;
	}
	
	ifstream reader;
	reader.open (MutHashFilePath);
	string line = "";
	while (getline(reader, line))
	{
		
		vector <string> temp = Split(line, ' ');
		unsigned long hash = HashToLong(temp[0]);
		MutantHashes[hash] = atoi(temp[1].c_str());
		hash = HashToLong(RevComp(temp[0]));
		MutantHashes[hash] = atoi(temp[1].c_str());
	}
	reader.close(); 
	
	reader.open (MutHashFilePathReference);
	line = "";
	while (getline(reader, line))
	{

		vector <string> temp = Split(line, ' ');
		unsigned long hash = HashToLong(temp[0]);
		MutantHashes[hash] = atoi(temp[1].c_str());
		hash = HashToLong(RevComp(temp[0]));
		MutantHashes[hash] = atoi(temp[1].c_str());
	}
	reader.close();			
	
	reader.open(ExcludeFilePath); 
	while (getline(reader, line))
	{
		vector <string> temp = Split(line, ' ');
		unsigned long hash = HashToLong(temp[0]);
		cout << "adding " << temp[0] << " with C=" << temp[1] << endl; 
		ExcludeHashes[hash] =  atoi(temp[1].c_str());
		//hash = HashToLong(RevComp(temp[0]));
		//ExcludeHashes[hash] =  atoi(temp[1].c_str());
	}
	reader.close();
	//***********************************************
	//cout << "Call is Reference Contigs.fa OutStub HashList MaxVarientSize" << endl;
	double vm, rss, MAXvm, MAXrss;
	MAXvm = 0;
	MAXrss = 0;
	process_mem_usage(vm, rss, MAXvm, MAXrss);
   	cout << "VM: " << vm << "; RSS: " << rss << endl;	

	
	int BufferSize = 1000;

	Reff.open(RefFile);

	ifstream ModelFile; 
	ModelFile.open (ModelFilePath); 
	if (ModelFile.is_open())
	{ cout << "ModelFile is open";}
	else
	{
		cout << "Error no model file given, not worring abou this now" << endl;
		//return -1;
	}

      	ifstream HashList;
      	HashList.open (HashListFile);
      	if ( HashList.is_open())
      	{	 cout << "HashList Open " << HashListFile << endl;}   //cout << "##File Opend\n";
      	else
      	{
      		cout << "Error, HashList could not be opened";
      		return -1;
      	}
	line = "";
	getline(HashList, line);
	cout << "line = " << line << endl; 
	char seperator = '\t';
	vector<string> temp = Split(line, seperator);
	if (temp.size() ==1){
		cout << "separator is not tab" << endl; 
		seperator = ' '; 
		 temp = Split(line, seperator);
	}
	else
		cout << "separator is tab" << endl;

	cout << "split = " << temp[0] << " and " << temp[1] << endl;
	HashSize = temp[0].size(); 
	if (temp.size() ==4)
	{ 
		HashSize = temp[3].length(); 
		Hash.insert(pair<string, int>(temp[3], atoi(temp[2].c_str())));
	
		while ( getline(HashList, line))
		{
			vector<string> temp = Split(line, seperator);
			Hash.insert(pair<string, int>(temp[3], atoi(temp[2].c_str())));
			Hash.insert(pair<string, int>(RevComp(temp[3]), atoi(temp[2].c_str()))); 
			//cout << "added pair " << temp[3] << "\t" << temp[2] << endl;
		}
		HashList.close(); 
		cout << "done with HashList" << endl;
	}
	else if (temp.size() ==2)
	{
		HashSize = temp[0].length();
		Hash.insert(pair<string, int>(temp[0], atoi(temp[1].c_str())));
		Hash.insert(pair<string, int>(RevComp(temp[0]), atoi(temp[1].c_str())));
		while ( getline(HashList, line))
		{
			vector<string> temp = Split(line, seperator);
			Hash.insert(pair<string, int>(temp[0], atoi(temp[1].c_str())));
			//cout << "added pair " << temp[3] << "\t" << temp[2] << endl;
		}
		HashList.close();
		cout << "done with HashList" << endl;
	}
	/*else if (temp.size() ==1)
	{
		vector<string> temp = Split(line, ' ');
		HashSize = temp[0].length();
		Hash.insert(pair<string, int>(temp[0], atoi(temp[1].c_str())));
		Hash.insert(pair<string, int>(RevComp(temp[0]), atoi(temp[1].c_str())));
		while ( getline(HashList, line))
		{
			vector<string> temp = Split(line, ' ');
			Hash.insert(pair<string, int>(temp[0], atoi(temp[1].c_str())));
			Hash.insert(pair<string, int>(RevComp(temp[0]), atoi(temp[1].c_str())));
			//cout << "added pair " << temp[3] << "\t" << temp[2] << endl;
		}
		HashList.close();
		cout << "done with HashList" << endl;
	}*/
	//map<string, int>::iterator it;
	//for ( it = Hash.begin(); it != Hash.end(); it++ )
	//{
	//	cout << "-"<<it->first<<"-" << "\t" << it->second << endl;
	//}

	ifstream SamFile;
       	if (samFile == "stdin")
	{
		cout << "Sam File is STDIN" << endl;
		SamFile.open ("/dev/stdin");
	}
	else
	{
		cout << "Sam File is " << samFile << endl;
		SamFile.open (samFile);
	}
	if ( SamFile.is_open())
	{      cout << "Sam File Opend\n";}
	else
	{
		cout << "Error, SamFile could not be opened";
		return 0;
	}

	string boom = outStub;
	VCFOutFile.open(boom+ ".vcf"); 
	BEDOutFile.open(boom+ ".vcf.bed"); 
	boom = "TempOverlap/" + boom; 
	BEDBigStuff.open(boom+ ".vcf.Big.bed");
	BEDNotHandled.open(boom+ ".vcf.NotHandled.bed");
	Invertions.open(boom+".vcf.invertions.bed");
	Translocations.open(boom+ ".vcf.Translocations");
	Translocationsbed.open(boom+ ".vcf.Translocations.bed");
	Unaligned.open(boom+"vcf.Unaligned");

	//write VCF header
	VCFOutFile << "##fileformat=VCFv4.1" << endl;
	VCFOutFile << "##fileDate=" << time(0) << endl;
	VCFOutFile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
	VCFOutFile << "##FORMAT=<ID=AK,Number=1,Type=Integer,Description=\"Alternte Kmer Count\">" << endl;
	VCFOutFile << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Kmer depth across the variant\">" << endl;
	VCFOutFile << "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Mode of reference kmer counts\">" << endl;
	VCFOutFile << "##FORMAT=<ID=AO,Number=1,Type=Integer,Description=\"Mode of alt kmer counts\">" << endl;
	VCFOutFile << "##INFO=<ID=FEX,Number=1,Type=String,Description=\"Filters failed and value\">" << endl;
	VCFOutFile << "##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand Bias of the aassembled contig\">" << endl;
	VCFOutFile << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV detected\">" << endl;
	VCFOutFile << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of SV detected\">" << endl; 
	VCFOutFile << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"END of SV detected\">" << endl; 
	VCFOutFile << "##INFO=<ID=AO,Number=1,Type=Integer,Description=\"Alternate allele observations, with partial observations recorded fractionally\">"<<endl;
	VCFOutFile << "##INFO=<ID=HD,Number=.,Type=String,Description=\"Hash counts for each k-mer overlapping the vareint, -1 indicates no info\">"<< endl;
	VCFOutFile << "##INFO=<ID=RN,Number=1,Type=String,Description=\"Name of contig that produced the call\">"<< endl;
	VCFOutFile << "##INFO=<ID=FS,Number=1,Type=String,Description=\"Full score, supporthing kmers \\ possible varient kmers based on sequence\">"<< endl;
	VCFOutFile << "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Mapping quality of the contig that created the call\">"<< endl;
	VCFOutFile << "##INFO=<ID=cigar,Number=1,Type=String,Description=\"Cigar string for the contig that created the call\">"<< endl;
	VCFOutFile << "##INFO=<ID=VT,Number=1,Type=String,Description=\"Varient Type\">"<< endl;	
	VCFOutFile << "##INFO=<ID=CVT,Number=1,Type=String,Description=\"Compressed Varient Type\">"<< endl;
	VCFOutFile << "##INFO=<ID=NR,Number=1,Type=Integer,Description=\"Number of total reads in target region\">" << std::endl;
	VCFOutFile << "##INFO=<ID=NH,Number=1,Type=Integer,Description=\"Number of alu heads in target region\">" << std::endl;
	VCFOutFile << "##INFO=<ID=NT,Number=1,Type=Integer,Description=\"Number of polyA tails in target region\">" << std::endl;
	VCFOutFile << "##INFO=<ID=LT,Number=1,Type=Integer,Description=\"Longest polyA tail in target region\">" << std::endl;
	VCFOutFile << "##INFO=<ID=TB,Number=1,Type=Integer,Description=\"Is tail left bound, right bound, or double bound\">" << std::endl;
	VCFOutFile << "##INFO=<ID=AS,Number=1,Type=Integer,Description=\"Number of alignment segments in the contig\">" << std::endl;
	VCFOutFile << "##FILTER=<ID=PA,Description=\"PoorAlignment\">" << std::endl;	
	VCFOutFile << "##FILTER=<ID=PLC,Description=\"Parents are at low coverage in this region, cannt be sure of genotype\">" << std::endl;
	VCFOutFile << "##FILTER=<ID=LCH,Description=\"Parents have hashes showing variant at low coverage, likely inherited\">" << std::endl;
	VCFOutFile << "##FILTER=<ID=SB,Description=\"Contig fails string bias filter\">" << std::endl; 
	VCFOutFile << "##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">" << std::endl;
	VCFOutFile << "##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">" << std::endl;

	VCFOutFile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
	
	string samplename = outStub.substr(0, outStub.find(".generator")); 
	VCFOutFile << samplename; 
	//VCFOutFile << outStub; 
	for(int i =0; i<ParentHashFilePaths.size(); i++)
	{
		string ParPath = argv[ParentHashFilePaths[i]]; 
		int startpos  = ParPath.find("overlap.asembly.hash.fastq.");
		int endpos = ParPath.find(".generator.Jhash");
		string Par = ParPath.substr(startpos+27, endpos-(startpos+27)); 
		VCFOutFile << "\t" << Par; 
		//VCFOutFile << "\tParent" << i; //argv[ParentHashFilePaths[i]]; 
	}
	VCFOutFile<< endl;

	int lines = 0;

	line="";
 
	unsigned long LongHash;

	cout << "Reading in Sam File"  << endl;
	map <string, int> Names; 
	vector <SamRead> reads; 
	int counter = 0; 
	while (getline(SamFile, line))
	{
		if (line.c_str()[0] == '@')
		{
			cout << " HEADER LINE = " << line << endl; 
			  vector <string> temp = Split(line, '\t');
			cout << temp[0] << endl;
			if (temp[0] == "@SQ")
			{
				cout << temp[1] << endl;
				vector <string> chr = Split(temp[1], ':');
				vector <string> len = Split(temp[2], ':');
	
				cout << "##contig=<ID=" <<  chr[1]<<",length=" << len[1] << ">"<< endl;
				VCFOutFile <<"##contig=<ID=" <<  chr[1]<<",length=" << len[1] << ">" << endl;
			}
		}
		else
		{
			cout << line << endl;
			counter ++; 
			SamRead read; 
			cout << "parse " << endl; 
			read.parse(line);
			//if (read.mapQual > 0)
			if (read.FlagBits[2] != 1)//read.flag != 4)
			{
				cout << "RefSeq" << endl;
				read.getRefSeq();
				cout << "peak" << endl;
				read.createPeakMap();
				cout << "ummm" << endl;
				int a; 
				string b;
				cout << "Aligned bases = " << read.CheckBasesAligned() << endl; 
				if (read.CheckBasesAligned() > 50 or read.CheckEndsAlign())
				{reads.push_back(read);}
				else 
				{cout << "SKIPPING Alignment" << endl; read.write();}
				if (counter%100 == 0)
					cout << "read " << counter << " entries " << char(13); 
			}
			//else do I want to track unaliged alignments? 
		}
	}
	cout << endl;
	cout << "Read in " << reads.size() << " reads " << endl;


	//  VCFOutFile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
	//VCFOutFile << outStub << endl;

	
	cout << "procesing split reads" << endl;
	for (int i = 0; i < reads.size(); i++) 
	{
		cout << "processing read " << reads[i].name << endl;
		if (reads[i].alignments.size() == 0)
		{
			reads[i].alignments.push_back(i);
			int count = 0; 
			for (int j = i+1; j < reads.size(); j++)
			{
				count++; 
				if (strcmp(reads[i].name.c_str(), reads[j].name.c_str()) == 0 && reads[i].pos  !=reads[j].pos)
				{
					cout << "found mate " << reads[j].name << endl;
					
					reads[i].alignments.push_back(j); 
					reads[j].first = false;
					
				}
				if (count > 100000)
					break; 
			}
		}
		for (int j = 0; j<reads[i].alignments.size(); j++)
			{reads[reads[i].alignments[j]].alignments = reads[i].alignments;} 

	}

	for (int i = 0; i < reads.size(); i++)
	{
		SamRead read = reads[i]; 
		cout << "starting work on  " << endl; 
		cout << read.name  << endl;
		cout << "alignments = " << read.alignments.size()  << endl; 
		
		if (strcmp(read.chr.c_str(), "*") == 0)
		{
			read.writetofile(Unaligned);
			Unaligned << endl;
		}
		//if  (read.first)
		{	
			//if it looks like a simple split read alignment colaps it into a single read
			
			if (read.first  and read.alignments.size() > 1)
			{
			  	cout << "picking two best alignments" << endl; 
				map <float,int> alignScores; 
				for (int j = 0; j < read.alignments.size(); j++){
				  	float score = (float) reads[read.alignments[j]].AlignScore; 
					while (not (alignScores.find(score) == alignScores.end())){
					  score = score * 1.0001; 
					}
					alignScores[score] = j; 
				}
				vector <int> goodPos; 
				std::map<float, int>::reverse_iterator it;	
				for ( it = alignScores.rbegin(); it != alignScores.rend(); it++ )
				{
				 	cout << it->first << " - " << it->second << endl; 
					goodPos.push_back(it->second); 
				}
				cout << "atempting colaps" << endl;
				cout << read.name << endl;
				vector<SamRead> R; 

				for(int j =0; j < read.alignments.size(); j++) //read.alignments.size(); j++)
				{
					//these better be sorted by position
					//if (reads[read.alignments[j]].chr == read.chr)
					if (j == goodPos[0] or j == goodPos[1])
					{
						R.push_back(reads[read.alignments[j]]);
						cout << reads[read.alignments[j]].name << endl;
					}
				}

				if (R.size() ==2 & /*R[0].chr == R[1].chr & */ R[0].mapQual > 0 or R[1].mapQual > 0)
				{
					read = BetterWay(R); 
				}
			}
			else if(read.first and read.alignments.size() >2)
			{
				BEDNotHandled << "too many alignments" << endl;
				BEDNotHandled << read.chr << "\t" << read.pos << "\t" << read.pos+read.seq.size() << "\t" << read.name << "\t" << read.cigar << endl;
				cout << "too many alignments" << endl;
				cout << read.chr << "\t" << read.pos << "\t" << read.pos+read.seq.size() << "\t" << read.name << "\t" << read.cigar << endl;
				for (int j = 0; j< read.alignments.size(); j++)
				{
					SamRead mate = reads[read.alignments[j]]; 
					cout << j << "\t" << mate.chr << "\t" << mate.pos << "\t" << mate.pos+mate.seq.size() << "\t" << mate.name << "\t" << mate.cigar << endl;
					 BEDNotHandled << mate.chr << "\t" << mate.pos << "\t" << mate.pos+mate.seq.size() << "\t" << mate.name << "\t" << mate.cigar << endl;
					mate.writetofile(BEDNotHandled);
				}
				BEDNotHandled << endl << endl;
			}
			if (read.mapQual > MinMapQual and read.alignments.size() <=2)
			{
				read.parseMutations(argv); 
			}
		}
	}
	cout << "lets start this" << endl; 
	//find big insertionsf
	for (int i = 0; i < reads.size()-1; i++)
	{
	 	int pos, pos2, kdep, kdep2;
		string InsStart, InsEnd; 
		if (reads[i].StartsWithAlign(pos, InsStart))
		{
					                        
			cout << "chekingBig " << endl; 
			reads[i].write(); 
			reads[i+1].write(); 
			if (reads[i].name != reads[i+1].name ){cout << "names are different" << endl; }else{cout << "NAMES ERROR" << endl;}
			if (reads[i].alignments.size() == 1){cout << "primary alignemts are 1" << endl;}else{cout << "PRIAMRY ALIGNEMT ERRPR " << reads[i].alignments.size() << endl; }
			if (reads[i].StrandBias < 0.9 and reads[i].StrandBias > 0.1) { cout << "Primary strand bias is good " << reads[i].StrandBias << endl; }else{cout << "STRAND BIAS ERROR " << reads[i].StrandBias << endl; }
			if (reads[i+1].StrandBias < 0.9 and reads[i+1].StrandBias > 0.1) { cout << "Primary strand bias is good " << reads[i+1].StrandBias << endl; }else{cout << "STRAND BIAS ERROR " << reads[i+1].StrandBias << endl; }
			if (reads[i+1].alignments.size() == 1){cout << "second alignemts are 1" << endl;}else{cout << "second ALIGNEMT ERRPR " << reads[i+1].alignments.size() << endl; }
			if (reads[i+1].EndsWithAlign(pos2, InsEnd)){cout << "second ends with clip" << endl;}else{cout << "SECOND NOPT CLOPPED" << endl;}
			if (reads[i].mapQual > 0 or reads[i+1].mapQual > 0){cout << "passed map qual" << endl; }else{cout << "failed map qual" << endl; }
			if (reads[i].StartsWithAlignAtPeak(pos, InsStart, kdep)){cout << "primary starts with clip at peak " << endl;}else{cout << "PRIMARY DOES NOT START AT PEAK" << endl; }
			if (reads[i+1].EndsWithAlignAtPeak(pos2, InsEnd, kdep2)){cout << "secondary align end with clip at peak" << endl;}else{cout << "SECONDARY DOES NOT START AT PEAK" << endl;}
			
			if (reads[i].name != reads[i+1].name and reads[i].alignments.size() == 1 and reads[i+1].alignments.size() == 1 and  reads[i].StartsWithAlign(pos, InsStart) and reads[i+1].EndsWithAlign(pos2, InsEnd) and (reads[i].StartsWithAlignAtPeak(pos, InsStart, kdep) or reads[i+1].EndsWithAlignAtPeak(pos2, InsEnd, kdep2)) and (reads[i].mapQual > 0 or reads[i+1].mapQual > 0) and abs(pos -pos2) <100 )
			{
			  	string D = "StrandBias"; 
				if ((reads[i+1].StrandBias < 1.0 and reads[i+1].StrandBias > 0.0) or (reads[i].StrandBias < 1.0 and reads[i].StrandBias > 0.0))
				  D = "DeNovo";
				cout       << reads[i].chr << "\t" << pos << "\t" << "LargeInsert" <<"-" << D /*"."*/  << "\t" << InsStart.c_str()[0] << "\t" << InsStart << "NNNNNNNNNNNNNNNNNNNN" << InsEnd << "\t" << kdep << "-" << kdep2 << "\t" << "." << "\t" << "INS" <<"RN=" << reads[i].name << ";MQ=" << reads[i].mapQual << ";cigar=" << reads[i].cigar << ";" << "RN=" << reads[i+1].name << ";MQ=" << reads[i+1].mapQual << ";cigar=" << reads[i+1].cigar << ";" <<"CVT=" << endl ;
			 	VCFOutFile  << reads[i].chr << "\t" << pos << "\t" << "LargeInsert" <<"-" << D /*"."*/  << "\t" << InsStart.c_str()[0] << "\t" << InsStart << "NNNNNNNNNNNNNNNNNNNN" << InsEnd << "\t" <<  kdep << "-" << kdep2 << "\t" << "." << "\t" << "INS" <<"RN=" << reads[i].name << ";MQ=" << reads[i].mapQual << ";cigar=" << reads[i].cigar << ";" << "RN=" << reads[i+1].name << ";MQ=" << reads[i+1].mapQual << ";cigar=" << reads[i+1].cigar << ";" << "CVT=" << ";Pos=" << pos << "-" << pos2 << "-" << pos2-pos;
			  
			 	string Genotype = "0/1"; 
			 	int MutRefMode = 1; 
			 	int MutAltMode = 1; 
			 	int LP = 1;
				int PC = 1; 
				int SB = 1; 
				VCFOutFile << "\tGT:DP:RO:AO:LP:PC:SB" << "\t" << Genotype << ":" << MutRefMode + MutAltMode << ":" << MutRefMode << ":" << MutAltMode << ":" << LP << ":" << PC << ":" << SB << endl;
	
			  	cout << "found INSERT " << endl; reads[i].write(); reads[i+1].write(); reads[i].writeVertical(); reads[i+1].writeVertical();
			}
		}
	}
	VCFOutFile.close(); 
	BEDOutFile.close();
	BEDBigStuff.close();
	BEDNotHandled.close();
	Invertions.close(); 
	cout << "finishing RUFUS.Interpret for " << outStub << std::endl;
	return 0; 
}
	
