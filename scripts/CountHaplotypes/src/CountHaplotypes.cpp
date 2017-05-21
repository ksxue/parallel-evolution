//============================================================================
// Name        : CountHaplotypes.cpp
// Version     : 2.0
// Description : 2.0 Implement haplotype inference from paired-end
//               reads listed sequentially in a SAM format file.
//           1.0 Given a BAM file, reference sequence,
//               and ordered list of 1-indexed sites of interest,
//               iterate through reads and
//               output genotypes at all sites of interest
//               as recorded by a single read, i.e. with linkage.
//               Treats all reads as single-end reads.
//============================================================================

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include <map>
#include <cstring>

using namespace std;

// RUN PARAMETERS
string VERSION="2.0";
string SAM="";
string QUERY="";
string OUTFILE="";
string CHR="";
int BASEQTHRESHOLD=20;
int MAPQTHRESHOLD=20;
int LEFTTRIM=0;
int RIGHTTRIM=0;
bool HEADER=false;

struct SAMRead_t {
	string QName;
	int Flag;
	string Chr;
	int Pos;
	int MapQ;
	string Cigar;
	int TLen;
	string Seq;
	string Quality;
	string ExpandedCigar;
};

// FUNCTIONS
int ArgsParse(int argc, char *argv[]);
void PrintUsage();
void PrintParameters();
void SetDebug();
int ReadMultiFasta(string filename,
		vector<string> *sequencenames,
		vector<string> *sequences);
vector<string> StringSplit(string s, char c);
string ExpandCIGAR(string cigar);
SAMRead_t ReadSAM(string line);

int main(int argc, char *argv[]) {

	//==================================================
	// Parse command-line arguments.
	//==================================================

	if(ArgsParse(argc, argv) != 0){
		PrintUsage();
		return 1;
	}

	PrintParameters();

	//==================================================
	// Read in query sites.
	//==================================================

	printf("Reading queries.\n");
	vector<int> QuerySites;
	// Open the file.
	ifstream fq(QUERY.c_str(), ios::in);

	// Check that file exists.
	if(!fq){
		printf("Error: query file does not exist.\n");
		return 1;
	}

	// Read in the file line by line.
	// Record sites of interest as zero-indexed positions along the chromosome.
	string line;
	while(getline(fq, line)){
		QuerySites.push_back(atoi(line.c_str())-1);
	}
	int NumQueries=QuerySites.size();

	// Close the file.
	fq.close();


	//==================================================
	// Read in BAM file and tally reads.
	//==================================================

	printf("Reading SAM file.\n");

	// Open the file.
	ifstream fin(SAM.c_str(), ios::in);
	ofstream fout(OUTFILE.c_str(), ios::out);

	if(fin && fout){

		// If the header option is turned on,
		// then print a header with the tab-delimited sites of interest.
		// These sites are one-indexed, as in the input query file.
		if(HEADER){
			for(unsigned int i=0; i<QuerySites.size(); i++){
				if(i<QuerySites.size()-1){
					fout << QuerySites[i]+1 << "\t";
				}
				else{
					fout << QuerySites[i]+1 << endl;
				}
			}
		}

		// Read in the file line by line.
		string line;

		// To take in pairs of reads at a time,
		// store the current read ID and a vector format to store reads.
		string CurrentReadID="";
		vector<SAMRead_t> ReadPair;
		int NumReads=0;

		while(getline(fin, line)){

			// Read the line as a SAM format line and store in the appropriate object.
			SAMRead_t NewRead=ReadSAM(line);

			// If the new read matches the current read ID,
			// then save it and continue.
			if(strcmp(NewRead.QName.c_str(),CurrentReadID.c_str())==0){
				ReadPair.push_back(NewRead);
				NumReads++;
				continue;
			}
			// If the new read does not match the previous read,
			// then save the read in the read vector,
			// and reset the read ID and read vector.
			else{

				// Determine whether the read pair is worth parsing.
				bool ParsePair=true;

				// Run through the criteria for parsing the read.
				for(unsigned int i=0; i<ReadPair.size(); i++){

					// Exclude reads that did not map.
					if(ReadPair[i].Cigar=="*"){
						ParsePair=false;
						continue;
					}

					// Exclude reads with more than one alignment.
					if(ReadPair[i].Flag>256){
						ParsePair=false;
						continue;
					}

					// Consider only reads that map to the specified chromosome
					// and exceed the specified minimum mapping quality.
					if(ReadPair[i].Chr!=CHR || ReadPair[i].MapQ<MAPQTHRESHOLD){
						ParsePair=false;
						continue;
					}

					// Exclude reads that contain indels.
					// Determine from the CIGAR string whether indels are present.
					bool indel=false;
					string ForbiddenOperations="IDP";
					for(unsigned int j=0; j<ReadPair[i].Cigar.size();j++){
						if(strchr(ForbiddenOperations.c_str(), ReadPair[i].Cigar[j])!=NULL){
							ParsePair=false;
							continue;
						}
					}

					// Consider only reads that, based on mapping position
					// and length of mapping region, cover at least some of the
					// region of interest.
					if(!(ReadPair[i].Pos+ReadPair[i].TLen>QuerySites[0] &&
							ReadPair[i].Pos<QuerySites[NumQueries-1])){
						continue;
					}

				}

				// If the read pair fails any of the criteria above,
				// then move on to the next line of the file.
				if(!ParsePair){
					// Save the new read ID,
					// reset the number of saved reads and the read array,
					// and save the new read.
					CurrentReadID=NewRead.QName;
					ReadPair.clear();
					ReadPair.push_back(NewRead);
					NumReads=1;
					continue;
				}

				// If the read is worth parsing,
				// then create an array to store the haplotype.
				// Initialize it with 'N'.
				bool HaplotypeNonEmpty=false;
				char Haplotype[NumQueries];
				for(unsigned int i=0; i<NumQueries; i++){
					Haplotype[i]='N';
				}

				// Expand the CIGAR string for each read in the pair.
				for(unsigned int i=0; i<ReadPair.size(); i++){
					ReadPair[i].ExpandedCigar=ExpandCIGAR(ReadPair[i].Cigar);
					// Verify that the expanded CIGAR string matches the read length.
					if(ReadPair[i].ExpandedCigar.size()!=ReadPair[i].Seq.size()){
						printf("CIGAR parsing error.\n");
						return 1;
					}
				}

				// Iterate through the sites of interest and output the genotypes
				// at those sites in this read.
				// Include only sites in the read that exceed the specified quality score.
				for(unsigned int i=0; i<NumQueries; i++){
					char genotype='N';
					for(unsigned int j=0; j<ReadPair.size(); j++){
						if(QuerySites[i]>ReadPair[j].Pos &&
								QuerySites[i]<ReadPair[j].Pos+ReadPair[j].Seq.size()){
							// Iterate along the read until you reach the site of interest.
							// Do not count sites that have been soft-clipped.
							int RefPos=ReadPair[j].Pos;
							for(unsigned int k=0; k<ReadPair[j].Seq.size(); k++){
								if(ReadPair[j].ExpandedCigar[k]=='M'){
									if(RefPos==QuerySites[i] &&
											ReadPair[j].Quality[k]>=BASEQTHRESHOLD){

										// Record the genotype once you reach the site.
										genotype=ReadPair[j].Seq[k];

										// Check that the genotypes of the reads
										// in the pair are concordant.
										// Otherwise, output 'N' at that site.
										if(Haplotype[i]=='N'){
											Haplotype[i]=genotype;
											HaplotypeNonEmpty=true;
										}
										else if(Haplotype[i]!=genotype){
											Haplotype[i]='N';
										}
									}
									RefPos++;
								}
							}
						}
					}
				}

				// If the haplotype is non-empty,
				// output it in tab-delimited form.
				if(HaplotypeNonEmpty){
					for(unsigned int i=0; i<NumQueries; i++){
						fout << Haplotype[i] << "\t";
					}
					fout << "\n";
				}

				// Save the new read ID,
				// reset the number of saved reads and the read array,
				// and save the new read.
				CurrentReadID=NewRead.QName;
				ReadPair.clear();
				ReadPair.push_back(NewRead);
				NumReads=1;
			}
		}
	}
	else{
		printf("Error: SAM file does not exist.\n");
		return 1;
	}

	// Close the file.
	fin.close();
	fout.close();

	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
	return 0;
}

// ArgsParse
// Parses command-line arguments.
// Returns 1 if any argument conditions are violated.
int ArgsParse(int argc, char *argv[]){

	// If the only argument is debug,
	// set all parameters to the debug state.
	if(argc==2 && strcmp(argv[1],"debug")==0){
		SetDebug();
		return 0;
	}

	// Ensure that there are an even number of arguments,
	// leaving aside the program name.
	// The only exception is if the only argument is debug,
	// which sets all parameters to the debug state.
	if((argc - 1) % 2 != 0){
		printf("Invalid number of arguments.\n");
		return 1;
	}
	// Check the structure of arguments.
	for(int i=1; i<argc; i++){
		// Verify that every other argument is a flag.
		if(i%2 != 0){
			if(argv[i][0] != '-' || strlen(argv[i])!=2){
				printf("Invalid use of argument flags.\n");
				return 1;
			}
		}
	}

	// Parse each pair of arguments.
	for(int i=0; i<(argc-1)/2; i++){

		string flag=argv[2*i+1];
		string arg=argv[2*i+2];

		// Parse the flag string.
		switch(flag[1]){
		// -s input SAM file
		case 's':
			SAM = arg;
			break;
		// -q list of query sites
		case 'i':
			QUERY = arg;
			break;
		// -c chromosome
		case 'c':
			CHR = arg;
			break;
		// -o output summary file
		case 'o':
			OUTFILE = arg;
			break;
		// -Q base quality threshold
		case 'Q':
			BASEQTHRESHOLD = atoi(arg.c_str());
			if(BASEQTHRESHOLD > 40 || BASEQTHRESHOLD < 0){
				printf("Invalid -Q base quality threshold.\n");
				return 1;
			}
			break;
		// -q mapping quality threshold
		case 'q':
			MAPQTHRESHOLD = atoi(arg.c_str());
			break;
		// -l left read trim
		case 'l':
			LEFTTRIM = atoi(arg.c_str());
			break;
		// -r right read trim
		case 'r':
			RIGHTTRIM = atoi(arg.c_str());
			break;
		// -h header
		case 'h':
			HEADER=true;
		}
	}

	// Check that the required arguments exist.
	if(SAM==""){
		printf("Invalid arguments. Specify SAM file.\n");
		return 1;
	}
	if(QUERY==""){
		printf("Invalid arguments. Specify query sites.\n");
		return 1;
	}
	if(CHR==""){
		printf("Invalid arguments. Specify chromosome of interest.\n");
		return 1;
	}
	if(OUTFILE==""){
		printf("Invalid arguments. Specify output file.\n");
		return 1;
	}
	return 0;
}

// PrintParameters
// When called, prints the parameters for the run.
void PrintParameters(){
	cout << "CountHaplotypes version " << VERSION << endl;
	cout << "RUN PARAMETERS" << endl;
	cout << "SAM file: " << SAM << endl;
	cout << "query: " << QUERY << endl;
	cout << "chromosome: " << CHR << endl;
	cout << "output file: " << OUTFILE << endl;
	cout << "header: " << HEADER << endl;
	cout << "base quality threshold: " << BASEQTHRESHOLD << endl;
	cout << "mapping quality threshold: " << MAPQTHRESHOLD << endl;
	cout << "left read trimming: " << LEFTTRIM << endl;
	cout << "right read trimming: " << RIGHTTRIM << endl;
	cout << endl;
}

// PrintUsage
// When called, prints the usage statement for this program.
void PrintUsage(){
	printf("\n\n");
	printf("Usage: CountHaplotypes -i query.txt -s input.sam -c chromosome -o out.haplotypes\n");
	printf("Given a SAM format file, a reference sequence,"
			"a list of sites of interest, \nand the chromosome"
			"on which they are located,"
			"iterate through the reads and establish haplotypes.\n");
	printf("\n");
	printf("  -i FILE\tordered list of 1-indexed sites of interest, one per line\n");
	printf("  -s FILE\tSAM-format file of reads, sorted so read pairs are adjacent to each other\n");
	printf("  -c STRING\tname of chromosome of interest\n");
	printf("  -o FILE\toutput list of haplotypes, one per line\n");
	printf("options (defaults in parentheses):\n");
	printf("  -Q INT\tminimum base quality for a base to be tallied [20]\n");
	printf("  -q INT\tminimum mapping quality for a read to be tallied [20]\n");
	printf("  -l INT\tnum bases to trim from 5' (left) end of each read, after soft clipping [0]\n");
	printf("  -r INT\tnum bases to trim from 3' (right) end of each read, after soft clipping [0]\n");
	printf("  -h print header line with query sites\n");
	printf("\n\n");
}

// SetDebug
// Sets all parameters to their debug state.
void SetDebug(){
	SAM="sam-v2.test";
	QUERY="in-v2.test";
	CHR="4-HA";
	OUTFILE="out-v2.test";
	HEADER=true;
}


//
// StringSplit
// Takes in a string and a character delimiter
// and returns a vector of strings split at that character.
vector<string> StringSplit(string s, char c){
	vector<string> splits;
	string s0;
	unsigned int i=0;

	while(i < s.length()){
		// Skip through delimiter characters at the beginnings of lines.
		while(s[i] == c && i < s.length() - 1){
			i++;
		}
		// Iterate through actual characters until you encounter c.
		while(i < s.length() && s[i] != c){
			s0 += s[i];
			i++;
		}
		// Once c is encountered, stop and save the string, then reset it.
		if(s0.size() > 0){
			splits.push_back(s0);
			s0 = "";
		}
		i++;
	}

	return splits;
}

//
// ExpandCIGAR
// Takes in a CIGAR string as specified by the SAM/BAM file standard
// and expands it into a string the length of the read describing its features.
// For instance, 1S5M4S becomes SMMMMMSSSS
string ExpandCIGAR(string cigar){
	string expanded="";
	string NumBases="";
	// All allowed operations characters.
	string CigarOperations="MIDNSHP=X";
	// Operations characters that sum to the length of the sequence.
	string CigarOperationsSubset="MIS=X";

	// Iterate through the CIGAR string and check for operation characters.
	for(unsigned int i=0; i<cigar.length(); i++){
		// Store non-operation characters that indicate the number of bases.
		if(strchr(CigarOperations.c_str(), cigar[i])==NULL){
			NumBases+=cigar[i];
		}
		// After reaching one of the useful operation characters,
		// add the appropriate number of characters to the expanded string.
		else if(strchr(CigarOperationsSubset.c_str(), cigar[i])!=NULL){
			unsigned int Num=atoi(NumBases.c_str());
			for(unsigned int j=0; j<Num; j++){
				expanded+=cigar[i];
			}
			NumBases="";
		}
		// If the operation character should not be added to the string,
		// reset the number of bases.
		else{
			NumBases="";
		}
	}
	return expanded;
}


//
// ReadSAM
// Given a line of a SAM-format file,
// parse that line and return a SAMRead_t type object.
SAMRead_t ReadSAM(string line){

	// Split the tab-delimited line and load it into BAM format.
	vector<string> fields=StringSplit(line,'\t');

	// Create a new SAMRead_t object.
	SAMRead_t Read;

	// Store the information in the appropriate formats.
	// Fields are hard-coded based on BAM file format.
	// Convert sequence positions from one-indexed to zero-indexed.
	Read.QName=fields[0];
	Read.Flag=atoi(fields[1].c_str());
	Read.Chr=fields[2];
	Read.Pos=atoi(fields[3].c_str())-1;
	Read.MapQ=atoi(fields[4].c_str())-1;
	Read.Cigar=fields[5];
	Read.TLen=atoi(fields[8].c_str());
	Read.Seq=fields[9];
	Read.Quality=fields[10];

	return Read;
}
