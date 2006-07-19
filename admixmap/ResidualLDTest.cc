/* 
 *   ADMIXMAP
 *   ResidualLDTest.cc 
 *   Class to implement score tests for residual allelic association between adjacent pairs of linked loci
 *   Copyright (c) 2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "ResidualLDTest.h"
#include "gsl/gsl_cdf.h" 
#include "IndividualCollection.h"

ResidualLDTest::ResidualLDTest(){
  Score = 0;
  Info = 0;
  SumScore= 0;
  SumScore2 = 0;
  SumInfo = 0;

  options = 0;
  individuals = 0;
  rank = 0;
  worker_rank = 0;
  NumWorkers = 1;
#ifdef PARALLEL
  Comm = 0;
#endif
}

ResidualLDTest::~ResidualLDTest(){
  if(test){
    for(unsigned j = 0; j < Lociptr->GetNumberOfChromosomes(); ++j){
      unsigned NumberOfLoci = Lociptr->GetSizeOfChromosome(j);
      for(unsigned k = 0; k < NumberOfLoci-1; ++k){
	delete[] Score[j][k];
	delete[] Info[j][k];
	if(rank==0){
	  delete[] SumScore[j][k];
	  delete[] SumScore2[j][k];
	  delete[] SumInfo[j][k];
	}
      }
      delete[] Score[j];
      delete[] Info[j];
      if(rank==0){
	delete[] SumScore[j];
	delete[] SumScore2[j];
	delete[] SumInfo[j];
      }
    }
    if(rank==0){
      delete[] SumScore;
      delete[] SumScore2;
      delete[] SumInfo;

    }
    delete[] Score;
    delete[] Info;
  }
}

#ifdef PARALLEL
void ResidualLDTest::SetComm(const MPI::Intracomm* c, const std::vector<std::string>* locuslabels){
  Comm = c;
  rank = Comm->Get_rank();
  worker_rank = rank - 1;
  NumWorkers = Comm->Get_size()-1;
  LocusLabels = locuslabels;
}
#endif

void ResidualLDTest::Initialise(AdmixOptions* op, const IndividualCollection* const indiv, const Genome* const Loci, 
				LogWriter &Log){
  options = op;
  individuals = indiv;
  chrm = Loci->getChromosomes();
  Lociptr = Loci;
  Log.setDisplayMode(Quiet);
  if(worker_rank==-1) worker_rank =indiv->getSize();//to stop master iterating through individuals
  test = options->getTestForResidualAllelicAssoc();

  if(test){
    //open output file
    if(rank==0)OpenFile(Log, &outputfile, options->getResidualAllelicAssocScoreFilename(), "Tests for residual allelic association", true);
#ifdef PARALLEL
    dimresallelescore = 0;
    dimresalleleinfo = 0;
#endif
    if(rank==0){
      SumScore = new double**[Lociptr->GetNumberOfChromosomes()];
      SumScore2 = new double**[Lociptr->GetNumberOfChromosomes()];
      SumInfo = new double**[Lociptr->GetNumberOfChromosomes()];
    }
    Score = new double**[Lociptr->GetNumberOfChromosomes()];
    Info = new double**[Lociptr->GetNumberOfChromosomes()];

    for(unsigned j = 0; j < Lociptr->GetNumberOfChromosomes(); ++j){
      unsigned NumberOfLoci = Lociptr->GetSizeOfChromosome(j);
      if(rank==0){      
	SumScore[j] = new double*[NumberOfLoci-1];
	SumScore2[j] = new double*[NumberOfLoci-1];
	SumInfo[j] = new double*[NumberOfLoci-1];
      }
      Score[j] = new double*[NumberOfLoci-1];
      Info[j] = new double*[NumberOfLoci-1];

      for(unsigned k = 0; k < NumberOfLoci-1; ++k){
#ifdef PARALLEL
	unsigned dim = 1;
	dimresallelescore += dim;
	dimresalleleinfo += dim*dim;
#else
	int locus = chrm[j]->GetLocus(k);
	unsigned dim = ((*Lociptr)(locus)->GetNumberOfStates()-1) * ((*Lociptr)(locus+1)->GetNumberOfStates()-1);
#endif
	if(rank==0){
	  SumScore[j][k] = new double[dim];
	  fill(SumScore[j][k], SumScore[j][k]+dim, 0.0);
	  SumScore2[j][k] = new double[dim*dim];
	  fill(SumScore2[j][k], SumScore2[j][k]+dim*dim, 0.0);
	  SumInfo[j][k] = new double[dim*dim];
	  fill(SumInfo[j][k], SumInfo[j][k]+dim*dim, 0.0);
	}
	Score[j][k] = new double[dim];
	Info[j][k] = new double[dim*dim];
      }

    }
#ifdef PARALLEL
    sendresallelescore = new double[dimresallelescore];
    sendresalleleinfo = new double[dimresalleleinfo];
    if(rank==0){
      recvresallelescore = new double[dimresallelescore];
      recvresalleleinfo = new double[dimresalleleinfo];
    }
#endif
  }
  
}

void ResidualLDTest::Reset(){
  //resets arrays holding sums of scores and info over individuals to zero; invoked at start of each iteration after burnin.

  if(test){
    for(unsigned j = 0; j < Lociptr->GetNumberOfChromosomes(); ++j){
      for(unsigned k = 0; k < Lociptr->GetSizeOfChromosome(j)-1; ++k){
#ifdef PARALLEL
	unsigned dim = 1;
#else
	int locus = chrm[j]->GetLocus(k);
	unsigned dim = ((*Lociptr)(locus)->GetNumberOfStates()-1) * ((*Lociptr)(locus+1)->GetNumberOfStates()-1);
#endif
	fill(Score[j][k], Score[j][k]+dim, 0.0);
	fill(Info[j][k], Info[j][k]+dim*dim, 0.0);
      }
    }

  }
}

void ResidualLDTest::Update(const array_of_allelefreqs& AlleleFreqs){
  int abslocus = 0;
  if(worker_rank < individuals->getSize())
  for(unsigned c = 0; c < Lociptr->GetNumberOfChromosomes(); ++c){
    for(unsigned j = 0; j < Lociptr->GetSizeOfChromosome(c)-1; ++j){
      UpdateScoresForResidualAllelicAssociation(c, j, AlleleFreqs[abslocus], AlleleFreqs[abslocus+1]);
      ++abslocus;
    }
    ++abslocus;//for last locus on chrm
  }
#ifdef PARALLEL
  //pack score and info into arrays ready to send
  int scoreindex = 0, infoindex = 0;
  for(unsigned j = 0; j < Lociptr->GetNumberOfChromosomes(); ++j){
    for(unsigned k = 0; k < Lociptr->GetSizeOfChromosome(j)-1; ++k){
      unsigned dim = 1;
      copy(Score[j][k], Score[j][k]+dim, sendresallelescore+scoreindex);
      copy(Info[j][k], Info[j][k]+dim, sendresalleleinfo+infoindex);
      scoreindex += dim;
      infoindex += dim*dim;
    }
  }
  //reduce into receive arrays on master process
  Comm->Barrier();
  Comm->Reduce(sendresallelescore, recvresallelescore, dimresallelescore, MPI::DOUBLE, MPI::SUM, 0);
  Comm->Reduce(sendresalleleinfo, recvresalleleinfo, dimresalleleinfo, MPI::DOUBLE, MPI::SUM, 0);

  if(rank==0){
    //accumulate score, square of score and info on master process
    scoreindex = 0; infoindex = 0;
    for(unsigned c = 0; c < Lociptr->GetNumberOfChromosomes(); ++c)
      for(unsigned k = 0; k < Lociptr->GetSizeOfChromosome(c)-1; ++k){
      unsigned dim = 1;

	for(unsigned j = 0; j < dim; ++j){
	  SumScore[c][k][j] += recvresallelescore[scoreindex + j];
	  for(unsigned jj = 0; jj < dim; ++jj){
	    SumScore2[c][k][j*dim +jj] += recvresallelescore[scoreindex + j]*recvresallelescore[scoreindex + jj];
	    SumInfo[c][k][j*dim +jj] = recvresalleleinfo[infoindex + j*dim + jj];
	  }
	}
	scoreindex += dim;
	infoindex+= dim*dim;
      }
  }
#else

  //accumulate score, square of score and info
  for(unsigned c = 0; c < Lociptr->GetNumberOfChromosomes(); ++c)
    for(unsigned k = 0; k < Lociptr->GetSizeOfChromosome(c)-1; ++k){
      int locus = chrm[c]->GetLocus(k);
      unsigned dim = ((*Lociptr)(locus)->GetNumberOfStates()-1) * ((*Lociptr)(locus+1)->GetNumberOfStates()-1);

      for(unsigned j = 0; j < dim; ++j){
	SumScore[c][k][j] += Score[c][k][j];
	for(unsigned jj = 0; jj < dim; ++jj){
	  SumScore2[c][k][j*dim +jj] += Score[c][k][j]*Score[c][k][jj];
	  SumInfo[c][k][j*dim +jj] += Info[c][k][j*dim + jj];
	}
      }
    }
  #endif
}


int delta(int i, int j){
  return (i == j);
}

void ResidualLDTest::UpdateScoresForResidualAllelicAssociation(int c, int locus,  
							   const double* const AlleleFreqsA, const double* const AlleleFreqsB){
  int abslocus = chrm[c]->GetLocus(locus);//number of this locus
#ifdef PARALLEL
  int M = 1, N = 1;
#else
  int M = (*Lociptr)(abslocus)->GetNumberOfStates()-1;
  int N = (*Lociptr)(abslocus+1)->GetNumberOfStates()-1;
#endif
  int dim = M*N;
//   if(dim == 1)UpdateScoresForResidualAllelicAssociation_1D(c, locus, AlleleFreqsA, AlleleFreqsB);
//   else
{
    
    //int Populations = options->getPopulations();

    int count = 0;
    int ancA[2];//ancestry at A
    int ancB[2];//ancestry at B

    for(int i = worker_rank; i < individuals->getSize(); i += NumWorkers){
      Individual* ind = individuals->getIndividual(i);
      if(!ind->GenotypeIsMissing(abslocus)){//skip loci with missing genotypes as hap pairs have not been sampled for these
	ind->GetLocusAncestry(c, locus, ancA);
	ind->GetLocusAncestry(c, locus+1, ancB);
	const int* hA = ind->getSampledHapPair(abslocus);//realized hap pair at locus A
	const int* hB = ind->getSampledHapPair(abslocus+1);//realized hap pair at locus B
	for(int g = 0; g < 2; ++g){
	  //if(ancA[g] == ancB[g]){ 
	  ++count;//count number of gametes with ancestry states the same at both loci
	  
	  for(int i = 0; i < M; ++i)
	    for(int j = 0; j < N; ++j){//i and j index rows
	      //update score
	      Score[c][locus][i*N +j] += ( delta(hA[g], i) - delta(hA[g], M) ) * ( delta(hB[g], j) - delta(hB[g], N) ) //observed
		- ( AlleleFreqsA[i + ancA[g]*(M+1)] - AlleleFreqsA[M + ancA[g]*(M+1)] ) 
		* ( AlleleFreqsB[j + ancB[g]*(N+1)] - AlleleFreqsB[N + ancB[g]*(N+1)] );// - expected
	      
	      for(int m = 0; m < M ; ++m)for(int n = 0; n < N; ++n){//m and n index columns
		//update info
		Info[c][locus][(i*N+j) * dim + (m*N+n)] += 
		  ( delta(i,m)*AlleleFreqsA[i + ancA[g]*(M+1)] + AlleleFreqsA[M + ancA[g]*(M+1)] ) 
		  * ( delta(j,n)*AlleleFreqsB[j + ancB[g]*(N+1)] + AlleleFreqsB[N + ancB[g]*(N+1)] );
	      }
	    }
	  
	}//end gamete loop
	//}//end condition (ancestry equal)
      }
    }//end ind loop
    for(int i = 0; i < dim*dim; ++i)
      Info[c][locus][i] *= (double)count;
  }//end else
}
void ResidualLDTest::UpdateScoresForResidualAllelicAssociation_1D(int c, int locus,  
							      const double* const AlleleFreqsA, const double* const AlleleFreqsB){
  int ancA[2];
  int ancB[2];
  //int Populations = options->getPopulations();
  int count = 0;

  int abslocus = chrm[c]->GetLocus(locus);
  //double sumh = 0.0, sumProbCoupling = 0.0;
  for(int i = worker_rank; i < individuals->getSize(); i += NumWorkers){
    Individual* ind = individuals->getIndividual(i);
    if(!ind->GenotypeIsMissing(abslocus)){//skip loci with missing genotypes as hap pairs have not been sampled for these
      ind->GetLocusAncestry(c, locus, ancA);
      ind->GetLocusAncestry(c, locus+1, ancB);
      const int* hA = ind->getSampledHapPair(abslocus);//realized hap pair at locus A
      const int* hB = ind->getSampledHapPair(abslocus+1);//realized hap pair at locus B
      
      for(int g = 0; g < 2; ++g){
	//if(ancA[g] == ancB[g]){ 
	++count;//count number of gametes with ancestry states the same at both loci
	int h = (hA[g] == hB[g]);//indicator for coupling
	double phiA = AlleleFreqsA[2*ancA[g]];//frequency of first allele in this gamete's ancestry state at locus A
	double phiB = AlleleFreqsB[2*ancB[g]];
	double ProbCoupling = phiA*phiB + (1 - phiA)*(1 - phiB);
	  
	//h - (1-h) evaluates to 1 if allele states equal, -1 if not
	//expected counts factorises to  --\/
	//sumh += h;
	//sumProbCoupling += ProbCoupling;
	Score[c][locus][0] += h - ProbCoupling; //h - (1-h) - (2*phiA - 1.0)*(2*phiB - 1.0);
	Info[c][locus][0] += ProbCoupling *(1 - ProbCoupling);
	//}//end condition on equal ancestry states
      }//end gamete loop
    }
  }//end individual loop

}

void ResidualLDTest::Output(int iterations, bool final){
  string sep = final ? "\t" : ",";//separator
  ofstream* outfile;

  if(test){
    //vector<string> labels;
    if(final){
      string filename(options->getResultsDir());
      filename.append("/ResidualLDTestFinal.txt");
      outfile = new ofstream(filename.c_str(), ios::out);
      *outfile << "Loci\tScore\tCompleteInfo\tObservedInfo\tPercentInfo\tdf\tChiSquared\tPValue\n";
    }else outfile = &outputfile;
    OutputTestsForResidualAllelicAssociation(iterations, outfile, final);
    //     for(unsigned int c = 0; c < Lociptr->GetNumberOfChromosomes(); c++ )
    //       for(unsigned j = 0; j < chrm[c]->GetSize()-1; ++j){
    // 	int abslocus = chrm[c]->GetLocus(j);
    // 	labels.clear();
    // 	labels.push_back("\"" + (*Lociptr)(abslocus)->GetLabel(0) + "/" + (*Lociptr)(abslocus+1)->GetLabel(0) + "\"");
    // 	int M = (*Lociptr)(abslocus)->GetNumberOfStates()-1;
    // 	int N = (*Lociptr)(abslocus+1)->GetNumberOfStates()-1;
    // 	OutputScoreTest(iterations, outfile, M*N, labels, SumScore[c][j], SumScore2[c][j], SumInfo[c][j], final);
    //       }
    if(final)delete outfile;
  }
  
}

void ResidualLDTest::OutputTestsForResidualAllelicAssociation(int iterations, ofstream* outputstream, bool final){
  //TODO: replace with use of function in base class
  double *score = 0, *ObservedInfo = 0;
  string separator = final? "\t" : ",";

  int abslocus = 0;
  for(unsigned int c = 0; c < Lociptr->GetNumberOfChromosomes(); c++ ){
    for(unsigned j = 0; j < Lociptr->GetSizeOfChromosome(c)-1; ++j){
#ifdef PARALLEL
      int M = 1, N = 1;
      const string label1 = (*LocusLabels)[abslocus];
      const string label2 = (*LocusLabels)[abslocus+1];
#else
      int M = (*Lociptr)(abslocus)->GetNumberOfStates()-1;
      int N = (*Lociptr)(abslocus+1)->GetNumberOfStates()-1;
      const string label1 = (*Lociptr)(abslocus)->GetLabel(0);
      const string label2 = (*Lociptr)(abslocus+1)->GetLabel(0);
#endif
      int dim = M*N;
      score = new double[dim];
      ObservedInfo = new double[dim*dim];
      double obsinfo = 0.0;
      double compinfo = 0.0;

      for(int i = 0; i < dim; ++i){
	score[i] = SumScore[c][j][i]/(double) iterations; //score
      }
      for(int i = 0; i < dim; ++i){
	for(int ii = 0; ii < dim; ++ii){
	  ObservedInfo[i*dim +ii] = SumInfo[c][j][i*dim+ii]/ (double) iterations//complete info
	    + score[i]*score[ii] - SumScore2[c][j][i*dim+ii]/(double)iterations;//-missing info
	}
	obsinfo += ObservedInfo[i*dim +i];//trace of Observed Info
	compinfo += SumInfo[c][j][i*dim+i]/ (double) iterations;//trace of Complete Info
      }

      //output labels
      *outputstream << "\"" << label1 << "/" << label2 << "\""<< separator;
      if(final)*outputstream << setiosflags(ios::fixed) << setprecision(3) //output 3 decimal places
			     << score[0] << separator << compinfo << separator <<obsinfo << separator;

      //compute chi-squared statistic
      try{
	double* VinvU = new double[dim];
	HH_solve(dim, ObservedInfo, score, VinvU);
	double chisq = 0.0;
	for(int i = 0; i < dim; ++i){
	  chisq += score[i] * VinvU[i];
	}
	delete[] VinvU;
	
	if(final)*outputstream << double2R((obsinfo*100)/compinfo) << separator<< dim << separator;

	if(chisq < 0.0){
	  *outputstream << "NA" << separator;
	  if(final) *outputstream << "NA" << separator;
	  *outputstream << endl;
	}
	else {
	  //compute p-value
	  double pvalue = gsl_cdf_chisq_Q (chisq, dim);
	  if(final)*outputstream << double2R(chisq) << separator << double2R(pvalue) << separator << endl;
	  else *outputstream << double2R(-log10(pvalue)) << separator << endl;
	}
      }
      catch(...){//in case ObservedInfo is rank deficient
	*outputstream  << "NA" << separator ;
	if(final)*outputstream << dim << separator << "NA" << separator << "NA" << separator;
	*outputstream << endl;
      }
      delete[] score;
      delete[] ObservedInfo;
      ++abslocus;
    }//end loop over loci on chromosome
    ++abslocus;//for last locus on chromosome
  }
}

void ResidualLDTest::ROutput(){
  if(test){
  int numPrintedIterations = (options->getTotalSamples() - options->getBurnIn()) / (options->getSampleEvery() * 10);
    vector<int> dimensions(3,0);
    dimensions[0] = 2;
    dimensions[1] = Lociptr->GetNumberOfCompositeLoci() - Lociptr->GetNumberOfChromosomes();
    dimensions[2] = (int)(numPrintedIterations);
    
    vector<string> labels(dimensions[0],"");
    labels[0] = "Loci";
    labels[1] = "log10Pvalue";

    R_output3DarrayDimensions(&outputfile, dimensions, labels);
  }
}
