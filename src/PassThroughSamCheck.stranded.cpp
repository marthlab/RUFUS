
/*This vertion imcorporates both copy number and mutation detection in 1
 *  * it needs to be run in two sptes, first the build, then the filter
 *   * it is split up to allow distribution to a cluster */

#include <unistd.h>
#include <ios>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <unordered_map>

using namespace std;

unordered_map <string, vector<string> > reads;  


const vector<string> Split(const string& line, const char delim) {
  vector<string> tokens;
  stringstream lineStream(line);
  string token;
  while ( getline(lineStream, token, delim) )
    tokens.push_back(token);
  return tokens;
}

int main (int argc, char *argv[])
{
  ifstream SamIn;         
  SamIn.open ("/dev/stdin");

  ofstream ChrOut;
  ChrOut.open (argv[1]);
  if (ChrOut.is_open())
    {}
  else
    {
      cout << "ERROR, Output file could not be opened -" << argv[1] << endl;
      return 0;
    }

  // ofstream mate1;
  // string m1n = argv[2]; 
  // m1n = m1n + ".mate1.fastq"; 
  // mate1.open (m1n);
  ofstream mate1(argv[2]);
  mate1.setf(std::ios::unitbuf); // Force flush after every insertion
  if (mate1.is_open())
    {}
  else
    {
      cout << "ERROR, Output file could not be opened -" << argv[1] << endl;
      return 0;
    }

  // ofstream mate2;
  // string m2n = argv[2]; 
  // m2n = m2n + ".mate2.fastq";
  // mate2.open (m2n);
  ofstream mate2(argv[3]);
  mate2.setf(std::ios::unitbuf);
  if (mate2.is_open())
    {}
  else
    {
      cout << "ERROR, Output file could not be opened -" << argv[1] << endl;
      return 0;
    }

  
  string L1;
  string current = "notachr";
  //string chr = "";
  

  while (getline(SamIn, L1))
    {
      //string name = "";
      int i = 0; 
      
      const char* L1_array = L1.c_str(); 
      int start =i;
      const char* tmpName = L1_array + i ; 
      //field 1
      while (L1_array[i] != '\t')
	{
	  //name+=L1_array[i]; 
	  ++i;
	}
      int lenName = i-start; 
      ++i;
      //string flag = "";
      //field 2
      int flagStart = i; 
      while (L1_array[i] != '\t')
	{
	  ++i;
	}
      int flagEnd = i;  
      ++i;
      //string chr = "";
      start = i; 
      const char* tmpChr = L1_array + i ;
      //field 3
      while (L1_array[i] != '\t')
	{
	  //chr+=L1_array[i];
	  ++i;
	}
      int lenChr = i-start; 
      ++i;
      //string pos = "";
      //field 4
      while (L1_array[i] != '\t')
	{
	  ++i;
	}
      ++i;
      //string something = "";
      //field 5
      while (L1_array[i] != '\t')
	{
	  ++i;
	}
      ++i;
      //string cigar = "";
      //field 6
      while (L1_array[i] != '\t')
	{
	  ++i;
	}
      ++i;
      //string something2 = "";
      //field 7
      while (L1_array[i] != '\t')
	{
	  ++i;
	}
      ++i;
      //string something3=  ""; 
      //field 8
      while (L1_array[i] != '\t')
	{
	  ++i;
	}
      ++i;
      //string something4=  "";
      //field 9
      while (L1_array[i] != '\t')
	{
	  ++i;
	}
      ++i;
      //string seq = "";
      start = i; 
      int seqStart = i; 
      const char* tmpSeq=L1_array + i ;
      //field 10
      while (L1_array[i] != '\t')
	{
	  //seq += L1_array[i];
	  ++i;
	}
      int seqEnd = i; 
      int lenSeq = i-start; 
      ++i;
      //string qual = "";
      start = i;
      int qualStart = i; 
      const char* tmpQual = L1_array + i ;
      //field 11
      while (L1_array[i] != '\t')
	{
	  //       qual += L1_array[i];
	  ++i;
	}
      int qualEnd = i; 
      int lenQual = i - start; 
      

      
      if (strncmp( tmpChr,  current.c_str(), lenChr) != 0 or lenChr != current.size() )
	{
	  ChrOut << current  << endl;
	  
	  current = string(tmpChr, lenChr); 
	}
      //cout << "@" << name << endl << seq << endl << "+" << endl << qual << endl;
      
      string flagstring = "";

      stringstream name; 
      name.write(tmpName, lenName);
      string namestr = name.str(); 
      if (reads.count(namestr) > 0)
      {
      	for (int i = flagStart; i< flagEnd; i++)
	{
	  flagstring+=L1_array[i];
	}
      
      	if ( 0 != (atoi(flagstring.c_str()) & (1 << 4))){
		mate1.write("@", 1);
		mate1.write(tmpName, lenName);
		mate1 << endl;
		for (int j =seqEnd-1; j>=seqStart; j--){
		  switch (L1_array[j]){
		  case 'A' : mate1 << 'T'; break;
		  case 'C' : mate1 << 'G'; break;
		  case 'G' : mate1 << 'C'; break; 
		  case 'T' : mate1 << 'A'; break; 
		  case 'N' : mate1 << 'N'; break;
		  }
		}
		mate1 << endl; 
		mate1 << "+" << endl;
		for (int j =qualEnd-1; j>=qualStart; j--){
		  mate1 << L1_array[j]; 
		}
		mate1 << endl; 
      	}	
      	else
	{
	  mate1.write("@", 1);
	  mate1.write(tmpName, lenName);
	  mate1 << endl; 
	  mate1.write(tmpSeq, lenSeq);
	  mate1 << endl << "+" << endl; 
	  mate1.write(tmpQual, lenQual); 
	  mate1 << endl;  
	}
        mate2.write("@", 1);
	//mate2.write("MATE", 4);
	mate2.write(tmpName, lenName);
	mate2 << endl; 
	mate2 << reads[namestr][0];
	mate2 << endl << "+" << endl;
	mate2 << reads[namestr][1]; 
	mate2 << endl; 
    	reads.erase(namestr); 
    }
    else
    {
        for (int i = flagStart; i< flagEnd; i++)
        {
          flagstring+=L1_array[i];
        }
        stringstream seq;
        stringstream qual; 
        if ( 0 != (atoi(flagstring.c_str()) & (1 << 4))){
                for (int j =seqEnd-1; j>=seqStart; j--){
                  switch (L1_array[j]){
                  case 'A' : seq << 'T'; break;
                  case 'C' : seq << 'G'; break;
                  case 'G' : seq << 'C'; break;
                  case 'T' : seq << 'A'; break;
                  case 'N' : seq << 'N'; break;
                  }
                }
                for (int j =qualEnd-1; j>=qualStart; j--){
                  qual << L1_array[j];
                }
               vector <string> temp; 
               temp.push_back(seq.str());
	       temp.push_back(qual.str()); 
	       reads[namestr] = temp; 
	}
        else
        {
          seq.write(tmpSeq, lenSeq);
          qual.write(tmpQual, lenQual);
          vector <string> temp;
	  temp.push_back(seq.str());
	  temp.push_back(qual.str());
	  reads[namestr] = temp;
	}
    }
  }
  ChrOut << current << endl;
  SamIn.close();
  return 0; 
}
