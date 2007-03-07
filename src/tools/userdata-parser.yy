/* vim:set ft=yacc: */

/**
 * Grammar for user genetic data. Recognized format looks like this:
 *
 * idno	dncase	rsXXXX	rsXXXX	...
 * 01	0	T:A	G:C	...
 * 02	1		G:C	...
 *
 * Columns are TAB-separated.
 */

%{
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include "getopt.h"
#include "itoa.h"
#include "UserData.h"
#include "UserDataIndividual.h"
#include "UserDataLocus.h"
#include "HapmixmapDataFormatter.h"
#include "userdata-parser.h"

using namespace std;
 
void yyerror(const char *error)
{
	cout << error << endl;
}

extern "C"
{
        int yywrap()
        {
                return 1;
        }
}

/**
 * Global data structures to hold the object data before they get
 * actually used.
 */

vector<UserDataLocus> loci;
unsigned int loci_count = 0;
vector<UserDataIndividual> indivs;
vector<base_genotype_pair_t> gpairs;
base_genotype_pair_t gpair;
UserData ud;
Legend *legend = NULL;
HapmixmapDataFormatter *hdf = NULL;

int yylex(void); 
int main(int, char **);
int yyparse(void);

extern int yydebug;



%}

/**
 * Data structure to temporarily hold values from the input.
 */
%union
{
  char *snp_name;
  char *text;
  int number;
  UserData *ud;
  UserDataIndividual *indiv;
}

/**
 * Recognized tokens.
 */

%token LINE_END COLUMN_SEPARATOR INDIVIDUAL_HEADER OUTCOME_HEADER COLON
%token <number> NUMBER
%token <snp_name> SNP_ID
%token <text> GENOTYPE
%token <text> GENOTYPE_BASE
%token <text> WORD
%type <indiv> data_line;
%type <number> outcome
%type <snp_name> snp_id
%type <text> alphanumeric_individual_id
%type <text> individual_id
%type <text> numeric_individual_id
%type <text> present_genotype
%type <ud> user_data;

%%

user_data:
	 header_line
	 data_lines
	 {
    ud = UserData();
		// Fill the individual with values
    ud.setLoci(loci);
    loci.clear();
		ud.setIndivs(indivs);
		indivs.clear();
	 }
	 ;

header_line:
	INDIVIDUAL_HEADER
	COLUMN_SEPARATOR
	OUTCOME_HEADER
	snp_headers
	end_of_line
	;

snp_headers: /* empty */
	|
	snp_headers
	COLUMN_SEPARATOR
	snp_id
	;

snp_id:
	SNP_ID
	{
		$$ = $1;
    loci.push_back(UserDataLocus(loci_count++, legend->getLocusPointerBySnp($1)));
	}
	;

end_of_line:
	LINE_END
	;

data_lines:
	/* empty */
	|
	data_lines
	data_line
	;

data_line:
	individual_id
	COLUMN_SEPARATOR
	outcome
	genotypes
	end_of_line
	{
		$$ = new UserDataIndividual($1, $3);
		// Insert the genotype data
		$$->setBaseGenotypes(gpairs);
		indivs.push_back(*$$);
		gpairs.clear();
	}
	;

genotypes:
	 /* empty */
	 |
	 genotypes
	 COLUMN_SEPARATOR
	 genotype_field
	 ;

genotype_field:
	missing_genotype
	|
	present_genotype
	;

present_genotype:
	GENOTYPE_BASE
	COLON
	GENOTYPE_BASE
	{
		gpair.present = 1;
		gpair.g[0] = $1[0];
		gpair.g[1] = $3[0];
		gpairs.push_back(gpair);
	}
	;

missing_genotype:
	/* empty */
	{
		gpair.present = 0;
		gpair.g[0] = 0;
		gpair.g[1] = 0;
		gpairs.push_back(gpair);
	}
	;

individual_id:
	numeric_individual_id
	|
	alphanumeric_individual_id
	;

numeric_individual_id:
	NUMBER
	{
		$$ = itoa($1);
	}
;

alphanumeric_individual_id:
	WORD
	{
		$$ = $1;
	}
	;

outcome:
       NUMBER
       {
	       $$ = $1;
       }
       ;
%%

void initParser(string legendFileName)
{
  legend = new Legend(legendFileName);
  hdf = new HapmixmapDataFormatter(&ud, legend);
}

int main(int argc, char **argv)
{
  int c;
  int usage;
  string legendFileName = "";
  string outputNameBase = "";
  while (1) {
    static struct option long_options[] = {
      /* These options set a flag. */
      {"help", no_argument, &usage, 1},
      // {"brief",              no_argument, &verbose_flag, 0},
      /* These options don't set a flag.
         We distinguish them by their indices. */
      // {"add", no_argument, 0, 'a'},
      // {"append", no_argument, 0, 'b'},
      // {"delete", required_argument, 0, 'd'},
      // {"create", required_argument, 0, 'c'},
      // {"file", required_argument, 0, 'f'},
      {"legend",  required_argument, 0, 'l'},
      {"output",  required_argument, 0, 'o'},
      {0, 0, 0, 0}
    };
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long(argc, argv, "l:", long_options,
        &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      /* If this option set a flag, do nothing else now. */
      if (long_options[option_index].flag != 0)
        break;
      /*
      printf("option %s",
             long_options[option_index].name);
      if (optarg)
        printf(" with arg %s", optarg);
      printf("\n");
      break;
      */

    case 'l':
      legendFileName = optarg;
      break;
    case 'o':
      outputNameBase = optarg;
      break;
      
    case '?':
      /* getopt_long already printed an error message. */
      break;

    default:
      abort();
    }
  }

  if ("" == legendFileName) {
    cerr << "Legend file name is required." << endl;
    exit(EXIT_FAILURE);
  }
  /*
  if ("" == outputNameBase) {
    cerr << "Legend file name is required." << endl;
    exit(EXIT_FAILURE);
  }
  */

  // Legend is needed to create loci
  initParser(legendFileName);
  yyparse();

  // TODO: Output file names as parameters
  ofstream genotypes_file;
  genotypes_file.open("user_genotypes.txt");
  hdf->renderFile(genotypes_file, 0);
  genotypes_file.close();

  ofstream locus_file;
  locus_file.open("user_loci.txt");
  hdf->renderFile(locus_file, 1);
  locus_file.close();

  ofstream outcome_file;
  outcome_file.open("user_outcome.txt");
  hdf->renderFile(outcome_file, 2);
  outcome_file.close();

  delete legend;
  delete hdf;
}
