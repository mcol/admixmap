#include "Legend.h"

#define RS_LIMIT_READ 0

/**
 * Constructor, reads a file with data, tab-separated columns
 * rsXXXX <position>  <allele1> <allele2> <chromosome>
 */
Legend::Legend(string fn)
: file_name(fn)
{
  ifstream data_ifstream(file_name.c_str());
  initialize(data_ifstream);
}

Legend::Legend(ifstream& data_ifstream)
{
  if (data_ifstream.is_open()) {
    initialize(data_ifstream);
  } else {
    cerr << "Can't open the data file: '" << file_name << "'." << endl;
    exit(EXIT_FAILURE);
  }
}

Legend::Legend(stringstream& ss)
{
  initialize(ss);
}

void Legend::initialize(istream& data_ifstream)
{
  string line, buf;
  stringstream ss;
  string rs, chromosome;
  unsigned long position, rs_no;
  char a1, a2;
  int counter = 0;
  locus_t *l;
  locus_pointer_t lp;
  locus_pos_t chr_position;
  initialized = false;

  while (!data_ifstream.eof()) {
    counter++;

    // Read data and insert a value into rs map.
    getline(data_ifstream, line);
    if (0 == line.size()) {
      continue;
    }
    stringstream ss(line);
    // Read the identifier
    ss >> rs;
    ss >> position;
    ss >> a1;
    ss >> a2;
    ss >> chromosome;

    l = new locus_t;
    rs_no = atoi(rs.substr(2, rs.size() - 1).c_str());
    l->rs = rs_no;
    l->position = position;
    l->a1 = a1;
    l->a2 = a2;
    l->chromosome = (char)atoi(chromosome.substr(3, chromosome.size() - 2).c_str());
    lp.p = l;
    snp_by_id[rs_no] = lp;
    chr_position.chromosome = l->chromosome;
    chr_position.position = position;
    snp_by_pos[chr_position] = lp;
//      dbg_locus_t(*(snp_by_id[rs_no].p));
    
    if (RS_LIMIT_READ && (0 == counter % (int)1e5)) {
    	cerr << "Rs no. " << counter << endl;
    }
    
    if (RS_LIMIT_READ && (counter > (int)1e5)) {
    	break;
    }

  }

  // Initialize the strand coding table.
  // 
  strand_coding_t forwCoding;
  forwCoding.c['A'] = 'A';
  forwCoding.c['C'] = 'C';
  forwCoding.c['G'] = 'G';
  forwCoding.c['T'] = 'T';
  strandCoding["forward"] = forwCoding;

  // A->T, T->A, C->G, G->C
  strand_coding_t revCoding;
  revCoding.c['A'] = 'T';
  revCoding.c['C'] = 'G';
  revCoding.c['G'] = 'C';
  revCoding.c['T'] = 'A';
  strandCoding["reverse"] = revCoding;
  
  initialized = true;
}

Legend::~Legend()
{
  if (initialized) {
    map<const unsigned long, locus_pointer_t>::iterator i;
    for(i = snp_by_id.begin(); i != snp_by_id.end(); i++) {
      delete i->second.p;
    }
    snp_by_id.clear();
    snp_by_pos.clear();
  }
}

unsigned long Legend::getSnpNumberByString(const string& snp_id)
{
  unsigned long snp_no;
  if (snp_id.size() < 3) {
    throw string("Too short SNP!");
  }
  snp_no = atoi(snp_id.substr(2, snp_id.size() - 1).c_str());
  return snp_no;
}

locus_t *Legend::getLocusPointerBySnp(const string& snp_id)
{
  unsigned long snp_no = getSnpNumberByString(snp_id);
  map<const unsigned long, locus_pointer_t>::iterator i;
  i = snp_by_id.find(snp_no);
  if (snp_by_id.end() == i) {
    cerr << "SNP '" << snp_id << "' doesn't exist!" << endl;
    throw string("Cannae find this SNP!");
    // exit(EXIT_FAILURE);
  } else {
    return i->second.p;
  }
}

locus_t Legend::getLocusBySnp(const string& snp_id)
{
  return *getLocusPointerBySnp(snp_id);
}

string Legend::getChromosomeBySnp(const string& snp_id)
{
  // checkSelf();
  // cerr << "'" << _snp.substr(2, _snp.size() - 1).c_str() << "'" << endl;
  // cerr << "Looking up snp number " << snp_no << endl;
  const locus_t locus = getLocusBySnp(snp_id);
  stringstream ss;
  string ret;
  ss << "chr" << (unsigned int)locus.chromosome;
  ss >> ret;
  return ret;
}

/**
 * Print locus contents
 */
void Legend::dbg_locus_t(locus_t l) {
  cout << l.rs << "\t"
    << l.position << "\t"
    << l.a1 << "\t"
    << l.a2 << "\t"
    << "chr" << (int)l.chromosome << endl;
}

//void Legend::checkSelf(void)
//{
//  if (!initialized) {
//    cerr << "Legend not initialized!" << endl;
//    exit(EXIT_FAILURE);
//  }
//}

short int Legend::getSnpNumberValueByBase(const string& snp_id, char allele)
{
  // checkSelf();
  const unsigned long snp_no = getSnpNumberByString(snp_id);
  map<const unsigned long, locus_pointer_t>::iterator i;
  i = snp_by_id.find(snp_no);
  if (snp_by_id.end() == i) {
    // TODO: SNP not found, raise an exception
    cerr << "SNP '" << snp_id << "' not found in the legend." << endl;
    exit(EXIT_FAILURE);
    return 0;
  }
  // Find out which allele it could be
  // TODO: Implement allele coding consistency check
  short int alleleNo = getSnpNumberValueByBaseAndEncoding(i->second.p, allele, "forward");
  if (0 == alleleNo) {
    alleleNo = getSnpNumberValueByBaseAndEncoding(i->second.p, allele, "reverse");
  }
  if (0 == alleleNo) {
    cerr << "Allele '" << allele << "' not found in any of the forward or reverse strand" << endl;
    cerr << "encodings in locus '" << snp_id << "'." << endl;
    exit(EXIT_FAILURE);
  }
  return alleleNo;
}

short int Legend::getSnpNumberValueByBaseAndEncoding(
    const locus_t *pLocusData,
    char allele,
    const string& encoding)
{
  if (pLocusData->a1 == strandCoding[encoding].c[allele]) {
    return 1;
  } else if (pLocusData->a2 == strandCoding[encoding].c[allele]) {
    return 2;
  } else {
    return 0;
  }
}

bool Legend::isInitialized()
{
  return initialized;
}

///**
// * Find out relation of two loci by rs identifiers. Used for sorting.
// */ 
//bool Legend::isLessThan(const string& snp1, const string& snp2)
//{
//  return (compareLoci(snp1, snp2) < 0);
//}
//
//bool Legend::isGreaterThan(const string& snp1, const string& snp2)
//{
//  return (compareLoci(snp1, snp2) > 0);
//}
  

//int Legend::compareLoci(const string& snp1, const string& snp2)
//{
//  locus_t l1 = getLocusBySnp(snp1);
//  locus_t l2 = getLocusBySnp(snp2);
//  if (l1.chromosome < l2.chromosome) {
//    return -1; 
//  } else if (l1.chromosome > l2.chromosome) {
//    return 1;
//  } else {
//    if (l1.position < l2.position) {
//      return -1;
//    } else if (l1.position > l2.position) {
//      return 1;
//    } else {
//      return 0;
//    }
//  }
//}
