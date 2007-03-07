#include "UserDataLocus.h"


UserDataLocus::UserDataLocus(
    // const string _id,
    unsigned int _idx,
    locus_t *locus)
{
  id = locus->rs;
  idx = _idx;
  coding = FORWARD_STRAND;
  flipped = false;
  monomorphic = false;
  number_alleles = 0;
  faulty = false;
  analyzed = false;
  vpUserData = NULL;
  locusData = locus;
}

void UserDataLocus::setUserData(UserData *pud)
{
  vpUserData = pud;
}

void UserDataLocus::confess(void)
{
  cout << "Locus: "  << id << ", idx: " << idx
    << ", coding "   << coding
    << ", flipped: " << flipped
    << endl;
  analyze();

  set<char>::iterator ai;
  cout << "Alleles (" << (unsigned int)number_alleles << "): ";
  for (ai = alleles.begin(); ai != alleles.end(); ai++) {
    cout << *ai;
  }
  cout << endl;
}

/**
 * Analyze the locus, check if it's monomorphic and whether it's forward
 * or reverse strand.
 */

void UserDataLocus::analyze()
{
  if (analyzed) {
    return;
  }
  // Find all alleles
  vector<UserDataIndividual> indivs = vpUserData->getIndividuals();
  vector<UserDataIndividual>::iterator i;
  base_genotype_pair_t bgp;
  
  for (i = indivs.begin(); i != indivs.end(); i++) {
    // Check if the allele is already in the set
    bgp = i->getBaseGenotype(idx);
    for (int j = 0; j < 2; j++) {
      if (!bgp.present) {
        // Genotype is missing
        continue;
      } else if (alleles.find(bgp.g[j]) != alleles.end()) {
        // Allele is already in the set, do nothing
      } else {
        alleles.insert(bgp.g[j]);
      }
    }
    analyzed = true;
  }
  
  number_alleles = alleles.size();
  switch (number_alleles) {
    case 1:
      monomorphic = true;
      break;
    case 2:
      monomorphic = false;
    default:
      faulty = true;
      break;
  }
}

string& UserDataLocus::getId()
{
  return id;
}

UserDataLocus::~UserDataLocus()
{
}

bool UserDataLocus::isPresent(vector<UserDataIndividual>::iterator& i)
{
  return i->getBaseGenotypes()[idx].present;
}

char UserDataLocus::getBaseAllele(vector<UserDataIndividual>::iterator& i, unsigned char an)
{
  if (an != 0 && an != 1) {
    return 'X';
  }
  return i->getBaseGenotypes()[idx].g[an];  
}

short int UserDataLocus::getNumberAllele(
    vector<UserDataIndividual>::iterator& i,
    unsigned char an,
    Legend& legend)
{
  char g = getBaseAllele(i, an);
  return legend.getSnpNumberValueByBase(id, g);
}

bool UserDataLocus::operator < (const UserDataLocus& other) const {
  return (compare(other) < 0) ? true : false;
}

bool UserDataLocus::operator > (const UserDataLocus& other) const {
  return (compare(other) > 0) ? true : false;
}

bool UserDataLocus::operator == (const UserDataLocus& other) const {
  return (compare(other) == 0) ? true : false;
}

bool UserDataLocus::operator != (const UserDataLocus& other) const {
  return (compare(other) != 0) ? true : false;
}

int UserDataLocus::compare(const UserDataLocus& other) const {
  // Chromosome takes precedence
        
  if (locusData->chromosome < other.locusData->chromosome) {
    return -1;
  } else if (locusData->chromosome > other.locusData->chromosome) {
    return 1;
  } else {
    // The same chromosome, loot at the position
    if (locusData->position < other.locusData->position) {
      return -1;
    } else if (locusData->position > other.locusData->position) {
      return 1;
    } else {
      return 0;
    }
  }
}

unsigned int UserDataLocus::getPosition()
{
  return locusData->position;
}

int UserDataLocus::getChromosome()
{
  return locusData->chromosome;
}

char UserDataLocus::getAllele1()
{
  return locusData->a1;
}

char UserDataLocus::getAllele2()
{
  return locusData->a2;
}
