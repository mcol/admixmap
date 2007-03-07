#include "HapmixmapDataFormatter.h"

HapmixmapDataFormatter::HapmixmapDataFormatter()
{
  vpUserData = 0;
}

HapmixmapDataFormatter::HapmixmapDataFormatter(
    UserData *pUserData,
    Legend *pLegend)
    : DataFormatter(pUserData, pLegend) // super class constructor
{
  // DataFormatter::DataFormatter(pUserData, pLegend);
}

HapmixmapDataFormatter::~HapmixmapDataFormatter()
{
}

/**
 * Hapmixmap data consists of 3 files:
 *
 * 1. genotypes
 * 2. loci
 * 3. outcome
 */
int HapmixmapDataFormatter::getNoFiles()
{
  return 3;
}

void HapmixmapDataFormatter::renderFile(ostream& out, int file_no)
{
  switch (file_no) {
    case 0:
      // Genotypes file
      render_genotypes(out, *vpLegend);
      break;
    case 1:
      // Locus file
      render_loci(out);
      break;
    case 2:
      // Outcome file
      render_outcome(out);
      break;
    default:
      // TODO: Shouldn't happen, perhaps raise an exception?
      break;
  }
}

void HapmixmapDataFormatter::render_genotypes(ostream& out, Legend& legend)
{
  if (NULL == vpUserData) {
    cerr << "vpUserData is empty" << endl;
    exit(EXIT_FAILURE);
  }
  int a1, a2, tmp;
  // TODO: Render genotypes
  out << "idno";
  // All the loci
  vector<UserDataLocus>::iterator i;
  vector<UserDataLocus> loci = vpUserData->getLoci();
  for (i = loci.begin(); i != loci.end(); i++) {
    out << "\t" << i->getId();
  }
  out << endl;
  
  // Output all the individuals
  vector<UserDataIndividual>::iterator j;
  vector<base_genotype_pair_t>::iterator gi;
  vector<base_genotype_pair_t> gpairs;
  vector<UserDataIndividual> indivs = vpUserData->getIndividuals();
  for (j = indivs.begin(); j != indivs.end(); j++) {
    out << j->getId();
    for (i = loci.begin(); i != loci.end(); i++) {
      out << "\t";
      if (i->isPresent(j)) {
        // Map "2,1" to "1,2"
        a1 = i->getNumberAllele(j, (unsigned char)0, legend);
        a2 = i->getNumberAllele(j, (unsigned char)1, legend);
        if (a1 > a2) {
          // swap
          tmp = a1; a1 = a2; a2 = tmp;
        }
        out << "\"" << a1 << "," << a2 << "\""; 
      } else {
        out << "\"0,0\"";
      }
    }
    out << endl;
  }
  return;
}

void HapmixmapDataFormatter::render_loci(ostream& out)
{
  // Render locusfile
  // Header
  out << "\"SNPid\"\t\"NumAlleles\"\t\"DistanceinMb\"" << endl;
  // For every locus
  vector<UserDataLocus> pLoci = vpUserData->getLoci();
  vector<UserDataLocus>::iterator i;
  double prevPosition;
  bool first = true;
  string distance;
  for (i = pLoci.begin(); i != pLoci.end(); i++) {
    // output a line for the locus
    // FIXME: Number of alleles is now hard-coded.
    out << i->getId()
        << "\t" << 2
        << "\t";
    if (first) {
      out << "#";
    } else {
      out << ((i->getPosition() - prevPosition) / 1e+6);
    }
    out << endl;
    prevPosition = i->getPosition();
    first = false;
  }
}

void HapmixmapDataFormatter::render_outcome(ostream& out)
{
  out << "dncase" << endl;
  vector<UserDataIndividual> indivs = vpUserData->getIndividuals();
  vector<UserDataIndividual>::iterator i;
  for (i = indivs.begin(); i != indivs.end(); i++) {
    out << i->getOutcome() << endl;
  }
}

