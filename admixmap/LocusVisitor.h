// *-*-C++-*-*
#ifndef LOCUS_VISITOR_H
#define LOCUS_VISITOR_H 1

/* forward declaration to avoid circular inclusions  */
class CompositeLocus;
class Chromosome;
class Genome;

/**
 * The LocusVisitor Class.
 *
 * This class defines an interface with which to write operations
 * to be performed on the elements in a Genome structure.
 * LocusVisitor lets you define a new operation without changing
 * the CompositeLocus, Chromosome or Genome classes.
 *
 * See LocusVisitorExample for an example of how this interface 
 * can be used.
 */

class LocusVisitor
{
public:
  virtual ~LocusVisitor() {}
  virtual void visitCompositeLocus(CompositeLocus&) = 0;
  virtual void visitChromosome(Chromosome&) = 0;
  virtual void visitGenome(Genome&) = 0;
};

#endif /* !defined LOCUS_VISITOR_H */
