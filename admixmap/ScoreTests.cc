#include "ScoreTests.h"

using namespace std;

ScoreTests::ScoreTests(){

  Matrix_d null_Matrix_d(1,1);
  Matrix_i null_Matrix_i(1,1);
  Vector_d null_Vector_d(1);
  Vector_i null_Vector_i(1);
  MatrixArray_d null_MatrixArray_d;
  MatrixArray_i null_MatrixArray_i;

  // these objects are composite-locus specific and could be moved there
  LocusLinkageScore = null_MatrixArray_d;
  LocusLinkageInfo = null_MatrixArray_d;
  SumLocusLinkageScore = null_MatrixArray_d;
  SumLocusLinkageScore2 = null_MatrixArray_d;
  SumLocusLinkageInfo = null_MatrixArray_d;
  AffectedsScore = null_MatrixArray_d;
  AffectedsVarScore = null_MatrixArray_d;
  AffectedsInfo = null_MatrixArray_d;
  LocusLinkageAlleleScore = null_MatrixArray_d;
  LocusLinkageAlleleInfo = null_MatrixArray_d;
  SumAffectedsScore2 = null_MatrixArray_d;
  SumAffectedsScore = null_MatrixArray_d;
  SumAffectedsVarScore = null_MatrixArray_d;
  SumAffectedsInfo = null_MatrixArray_d;
  SumLocusLinkageAlleleScore2 = null_MatrixArray_d;
  SumLocusLinkageAlleleScore = null_MatrixArray_d;
  SumLocusLinkageAlleleInfo = null_MatrixArray_d;
  ScoreWithinHaplotype = 0;
  InfoWithinHaplotype = 0;
  SumScoreWithinHaplotype = null_MatrixArray_d;
  SumScore2WithinHaplotype = null_MatrixArray_d;
  SumInfoWithinHaplotype = null_MatrixArray_d;

  U = null_Matrix_d; 
  I = null_Matrix_d; 
  SumU = null_Matrix_d; 
  SumU2 = null_Matrix_d;
  SumI = null_Matrix_d;

  ancestryAssociationScoreStream = 0;
  SNPsAssociationScoreStream = 0;
  affectedsOnlyScoreStream = 0;

  options = 0;
  individuals = 0;

}//end constructor

ScoreTests::~ScoreTests(){
  if( options->getTestForAllelicAssociation() ){
    delete [] ScoreWithinHaplotype;
    delete [] InfoWithinHaplotype;
  }
}

void ScoreTests::Initialise(AdmixOptions * op, IndividualCollection *indiv, Genome *L, Genome *c,
			    LogWriter *log){
  options = op;
  individuals = indiv;
  chrm = c;
  Lociptr=L;
  Logptr = log;

  pp.SetNumberOfElements( options->getPopulations());
  U.SetNumberOfElements( options->getPopulations(), indiv->getTargetSize() );
  SumU.SetNumberOfElements( options->getPopulations(), indiv->getTargetSize() );
  SumU2.SetNumberOfElements( options->getPopulations(), indiv->getTargetSize() );
  I.SetNumberOfElements( options->getPopulations(), indiv->getTargetSize() );
  SumI.SetNumberOfElements( options->getPopulations(), indiv->getTargetSize() );

  /*----------------------
  | admixture association |
   -----------------------*/
   if( options->getScoreTestIndicator() ){
    if ( strlen( options->getAssocScoreFilename() ) ){
        assocscorestream.open( options->getAssocScoreFilename(), ios::out );
      if( !assocscorestream ){
	Logptr->logmsg(true,"ERROR: Couldn't open assocscorefile\n");
	exit( 1 );}
      else {
	Logptr->logmsg(true,"Admixture score file: ");    
	Logptr->logmsg(true,options->getAssocScoreFilename());
	Logptr->logmsg(true,"\n");
	assocscorestream << setiosflags( ios::fixed );
      }
    }
    else{
      Logptr->logmsg(true,"No assocscorefile given\n");
      exit(1);}
  }

 
  if( options->getTestForAffectedsOnly() ){

    affectedsOnlyScoreStream = new ofstream(options->getAffectedsOnlyScoreFilename());
    if( !affectedsOnlyScoreStream ){
      Logptr->logmsg(true,"ERROR: Couldn't open affectedsonlyscorefile");
      Logptr->logmsg(true,options->getAffectedsOnlyScoreFilename());
      Logptr->logmsg(true,"\n");
      exit( 1 );//remove?
    }
    else{
    Logptr->logmsg(true,"Writing affected only tests for association to ");
    Logptr->logmsg(true,options->getAffectedsOnlyScoreFilename());
    Logptr->logmsg(true,"\n");
    *affectedsOnlyScoreStream << "structure(.Data=c(" << endl;
    }

    AffectedsScore.SetNumberOfElementsWithDimensions(Lociptr->GetNumberOfCompositeLoci(),
						     options->getPopulations(), 1 );
    AffectedsVarScore.SetNumberOfElementsWithDimensions(Lociptr->GetNumberOfCompositeLoci(),
							options->getPopulations(), 1 );
    AffectedsInfo.SetNumberOfElementsWithDimensions(Lociptr->GetNumberOfCompositeLoci(),
						    options->getPopulations(), 1 );
    SumAffectedsScore2.SetNumberOfElementsWithDimensions(Lociptr->GetNumberOfCompositeLoci(),
							 options->getPopulations(), 1 );
    SumAffectedsScore.SetNumberOfElementsWithDimensions(Lociptr->GetNumberOfCompositeLoci(),
							options->getPopulations(), 1 );
    SumAffectedsVarScore.SetNumberOfElementsWithDimensions(Lociptr->GetNumberOfCompositeLoci(),
							   options->getPopulations(), 1 );
    SumAffectedsInfo.SetNumberOfElementsWithDimensions(Lociptr->GetNumberOfCompositeLoci(),
						       options->getPopulations(), 1 );
  }

  /*-----------------------
  | Linkage with ancestry  |
   -----------------------*/
  if( options->getTestForLinkageWithAncestry() ){
    ancestryAssociationScoreStream = new ofstream(options->getAncestryAssociationScoreFilename());
    if( !ancestryAssociationScoreStream ){
      Logptr->logmsg(true,"ERROR: Couldn't open locusscorefile2\n");
      exit( 1 );//remove?
    }
      else{
	Logptr->logmsg(true,"Writing tests for locus linkage to ");
	Logptr->logmsg(true,options->getAncestryAssociationScoreFilename());
	Logptr->logmsg(true,"\n");
	*ancestryAssociationScoreStream << "structure(.Data=c(" << endl;
      }

    LocusLinkageScore.SetNumberOfElementsWithDimensions(Lociptr->GetNumberOfCompositeLoci(),
							2 * options->getPopulations(), 1 );
    LocusLinkageInfo.SetNumberOfElementsWithDimensions(Lociptr->GetNumberOfCompositeLoci(), 2 * options->getPopulations(),
						       2 * options->getPopulations() );
    SumLocusLinkageScore2.SetNumberOfElementsWithDimensions(Lociptr->GetNumberOfCompositeLoci(), options->getPopulations(), 1 );
    SumLocusLinkageScore.SetNumberOfElementsWithDimensions(Lociptr->GetNumberOfCompositeLoci(), options->getPopulations(), 1 );
    SumLocusLinkageInfo.SetNumberOfElementsWithDimensions(Lociptr->GetNumberOfCompositeLoci(), options->getPopulations(), 1 );
  }

  /*----------------------
  | Allelic association  |
   -----------------------*/
  if( options->getTestForAllelicAssociation() ){
    if( options->getAnalysisTypeIndicator() < 2 ){
      Logptr->logmsg(true,"ERROR: To test for allelic association analysistypeindicator must be 2 or 3.\n");
      exit(0);
    }
 
    genescorestream.open( options->getLocusScoreFilename(), ios::out );
    if( !genescorestream ){
      Logptr->logmsg(true,"ERROR: Couldn't open locusscorefile\n");
      exit( 1 );
    }
      else{
	Logptr->logmsg(false,"Test for allelic association written to ");
	Logptr->logmsg(false,options->getLocusScoreFilename());
	Logptr->logmsg(false,"\n");
	genescorestream << "structure(.Data=c(" << endl;
      }

    int kk;
    LocusLinkageAlleleScore.SetNumberOfElements( Lociptr->GetNumberOfCompositeLoci() );
    LocusLinkageAlleleInfo.SetNumberOfElements( Lociptr->GetNumberOfCompositeLoci() );
    SumLocusLinkageAlleleScore2.SetNumberOfElements( Lociptr->GetNumberOfCompositeLoci() );
    SumLocusLinkageAlleleScore.SetNumberOfElements( Lociptr->GetNumberOfCompositeLoci() );
    SumLocusLinkageAlleleInfo.SetNumberOfElements( Lociptr->GetNumberOfCompositeLoci() );
    ScoreWithinHaplotype = new MatrixArray_d[ Lociptr->GetNumberOfCompositeLoci() ];
    InfoWithinHaplotype = new MatrixArray_d[ Lociptr->GetNumberOfCompositeLoci() ];
    SumScoreWithinHaplotype.SetNumberOfElements( Lociptr->GetNumberOfCompositeLoci() );
    SumScore2WithinHaplotype.SetNumberOfElements( Lociptr->GetNumberOfCompositeLoci() );
    SumInfoWithinHaplotype.SetNumberOfElements( Lociptr->GetNumberOfCompositeLoci() );
    for( int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
      
      if((*Lociptr)(j)->GetNumberOfLoci() > 1 )
	kk = 1;
      else if((*Lociptr)(j)->GetNumberOfStates() == 2 )
	kk = 1;
      else
	kk = (*Lociptr)(j)->GetNumberOfStates();
      LocusLinkageAlleleScore(j).SetNumberOfElements( kk + options->getPopulations(), 1 );
      LocusLinkageAlleleInfo(j).SetNumberOfElements( kk + options->getPopulations(),
                                                     kk + options->getPopulations() );
      SumLocusLinkageAlleleScore2(j).SetNumberOfElements( kk, kk );
      SumLocusLinkageAlleleScore(j).SetNumberOfElements( kk, 1 );
      SumLocusLinkageAlleleInfo(j).SetNumberOfElements( kk, kk );
      if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){
	ScoreWithinHaplotype[ j ].SetNumberOfElementsWithDimensions( (*Lociptr)(j)->GetNumberOfLoci(), 1 + options->getPopulations(), 1 );
	InfoWithinHaplotype[ j ].SetNumberOfElementsWithDimensions( (*Lociptr)(j)->GetNumberOfLoci(), 1 + options->getPopulations(), 1 + options->getPopulations() );
	SumScoreWithinHaplotype(j).SetNumberOfElements( (*Lociptr)(j)->GetNumberOfLoci(), 1 );
	SumScore2WithinHaplotype(j).SetNumberOfElements( (*Lociptr)(j)->GetNumberOfLoci(), 1 );
	SumInfoWithinHaplotype(j).SetNumberOfElements( (*Lociptr)(j)->GetNumberOfLoci(), 1 );
      }
    }
}  

  /*---------------------------------
  | Misspecified allele frequencies  |
   ---------------------------------*/
  if( options->getTestForMisspecifiedAlleleFreqs() ){
    allelefreqscorestream.open( options->getAlleleFreqScoreFilename() );
    if( !allelefreqscorestream ){
      Logptr->logmsg(true,"ERROR: Couldn't open allelefreqscorefile\n");
      Logptr->logmsg(true,options->getAlleleFreqScoreFilename());
      Logptr->logmsg(true,"\n");
      exit( 1 );
    }
    else{
    Logptr->logmsg(true,"Writing score tests for mis-specified allele frequencies to ");
    Logptr->logmsg(true,options->getAlleleFreqScoreFilename());
    Logptr->logmsg(true,"\n");
    allelefreqscorestream << "structure(.Data=c(" << endl;
    }
  }

  /*-------------------------------------
  | Misspecified allele frequencies 2 (?) |
   --------------------------------------*/
  if( options->getTestForMisspecifiedAlleleFreqs2() ){
    allelefreqscorestream2.open( options->getAlleleFreqScoreFilename2() );
    if( !allelefreqscorestream2 ){
      Logptr->logmsg(true,"ERROR: Couldn't open allelefreqscorefile\n");
      Logptr->logmsg(true,options->getAlleleFreqScoreFilename2());
      Logptr->logmsg(true,"\n");
      exit( 1 );
    }
    else{
    Logptr->logmsg(true,"Writing score tests for mis-specified allele frequencies to ");
    Logptr->logmsg(true,options->getAlleleFreqScoreFilename2());
    Logptr->logmsg(true,"\n");
    allelefreqscorestream2 << "structure(.Data=c(" << endl;
    }
  }

  /*-------------------
  | SNPs in haplotype  |
   -------------------*/  
 if( strlen( options->getTestsForSNPsInHaplotypeOutputFilename() ) ){
    if( !options->getTestForAllelicAssociation() ){
      Logptr->logmsg(true,"Can't Test For SNPs In Haplotype if options->getTestForAllelicAssociation() is false\n");
      exit(1);
    }
    SNPsAssociationScoreStream = new ofstream( options->getTestsForSNPsInHaplotypeOutputFilename(), ios::out );
    if( !SNPsAssociationScoreStream ){
      Logptr->logmsg(true,"ERROR: Couldn't open testsforSNPsinhaplotypeoutputfilename\n");
      exit( 1 );
    }
    else{
    Logptr->logmsg(true,"Tests for associations of SNPs in haplotypes written to\n");
    Logptr->logmsg(true,options->getTestsForSNPsInHaplotypeOutputFilename());
    Logptr->logmsg(true,"\n");
    *SNPsAssociationScoreStream << "structure(.Data=c(" << endl;
    }
  }
    
}

//Initialise ergodic average score file
void ScoreTests::InitialiseAssocScoreFile(std::string *PLabels){
  if( options->getScoreTestIndicator() ){
    PopLabels = PLabels;
    assocscorestream << "Ergodic averages of score statistic for populations:\n";
    for( int i = 0; i < options->getPopulations(); i++ ){
      if(options->IsPedFile())
	assocscorestream << "\"" << PopLabels[i] << "\"" << " ";
      else
	assocscorestream << PopLabels[i] << " ";
      if( !i )
	assocscorestream << " ";
    }
    assocscorestream << "\ncomplete  missing   statistic  ";
    for( int i = 1; i < options->getPopulations(); i++ )
      assocscorestream << "complete  missing   statistic ";
    assocscorestream << endl;
  }

}

void ScoreTests::Reset(){

  I.SetElements(0);
  U.SetElements(0);
  if( options->getTestForAllelicAssociation() ){
    LocusLinkageAlleleScore.SetElements(0);
    LocusLinkageAlleleInfo.SetElements(0);
    for( int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
      if((* Lociptr)(j)->GetNumberOfLoci() > 1 ){
	ScoreWithinHaplotype[j].SetElements(0);
	InfoWithinHaplotype[j].SetElements(0);
      }
    }
  }
  if( options->getTestForLinkageWithAncestry() ){
    LocusLinkageScore.SetElements(0);
    LocusLinkageInfo.SetElements(0);
  }
  if( options->getTestForAffectedsOnly() ){
    AffectedsScore.SetElements(0);
    AffectedsVarScore.SetElements(0);
    AffectedsInfo.SetElements(0);
  }
}

void ScoreTests::Output2(Vector_d *alpha0, std::ofstream *LogFileStreamPtr){//should be renamed
  //Note: alpha0 = alpha[0]

  for( int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
    if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){

      for( int k = 0; k < options->getPopulations(); k++ )
	pp(k) = (*alpha0)(k) / alpha0->Sum();
      (*Lociptr)(j)->SetDefaultMergeHaplotypes( pp );
      if(options->IsPedFile())
	*LogFileStreamPtr << "\"" << (*Lociptr)(j)->GetLabel(0) << "\"" << endl;
      else
	*LogFileStreamPtr << (*Lociptr)(j)->GetLabel(0) << endl;
      for( int k = 0; k < (*Lociptr)(j)->GetNumberOfStates(); k++ ){
	*LogFileStreamPtr << k << " " << (*Lociptr)(j)->GetMergedHaplotype(k) << endl;
      }
      int kk = (*Lociptr)(j)->GetNumberOfMergedHaplotypes();
      LocusLinkageAlleleScore(j).SetNumberOfElements( kk + options->getPopulations(),
						      1 );
      LocusLinkageAlleleInfo(j).SetNumberOfElements( kk + options->getPopulations(),
						     kk + options->getPopulations() );
      SumLocusLinkageAlleleScore2(j).SetNumberOfElements( kk, kk );
      SumLocusLinkageAlleleScore(j).SetNumberOfElements( kk, 1 );
      SumLocusLinkageAlleleInfo(j).SetNumberOfElements( kk, kk );
    }
  }


}

void
ScoreTests::Update(int TT, Matrix_d* ExpY0,double lambda0)
//Note: TT = TargetType(0), ExpY0 = ExpectedY(0), lambda0=lambda(0)
{
  //the first part was formerly the function UpdateScoreIndLevel
  for( int i = 0; i < individuals->getSize(); i++ ){
    Individual* ind = individuals->getIndividual(i);
    if( options->getTestForMisspecifiedAlleleFreqs() )
      UpdateScoresForMisSpecOfAlleleFreqs( i );
    if( options->getAnalysisTypeIndicator() == 2 )
      if( !options->getScoreTestIndicator() ){
	if( options->getTestForAllelicAssociation() )
	  UpdateScoreForLinkage( ind, i, 0, 0, &LocusLinkageAlleleScore,
				 &LocusLinkageAlleleInfo, (*ExpY0)(i,0),lambda0);
	if( options->getTestForLinkageWithAncestry() ){
	  UpdateScoreForLinkage( ind, i, 0, 1, &LocusLinkageScore,
				 &LocusLinkageInfo, (*ExpY0)(i,0) ,lambda0);
	}
      }
      else
	UpdateScoreForAssociation_Continuous( i, ind->getAncestry(), &U, &I, (*ExpY0)(i,0),lambda0);
    else if( options->getAnalysisTypeIndicator() == 3 || 
	     options->getAnalysisTypeIndicator() == 4 ){
      if( !options->getScoreTestIndicator() ){
	if( options->getTestForAllelicAssociation() ){
	  UpdateScoreForLinkage( ind, i, 1, 0, &LocusLinkageAlleleScore,
				 &LocusLinkageAlleleInfo, (*ExpY0)(i,0),lambda0);
	}
	if( options->getTestForLinkageWithAncestry() ){
	  UpdateScoreForLinkage( ind, i, 1, 1, &LocusLinkageScore,
				 &LocusLinkageInfo, (*ExpY0)(i,0),lambda0);
	}
      }
      else{
	UpdateScoreForAssociation_Binary( i, ind->getAncestry(), &U, &I, (*ExpY0)(i,0));
      }
      if( options->getTestForAffectedsOnly() && (individuals->getTarget(0))(i,0) == 1 ){
	UpdateScoreForLinkageAffectedsOnly( ind); 
      }
    }
    else if( options->getAnalysisTypeIndicator() == 0 ){
      UpdateScoreForLinkageAffectedsOnly( ind); 
    }
    else if( options->getAnalysisTypeIndicator() == 5 ){
      if( options->getTestForAllelicAssociation() ){
	UpdateScoreForLinkage( ind, i, TT, 0,
                               &LocusLinkageAlleleScore,
			       &LocusLinkageAlleleInfo, (*ExpY0)(i,0),lambda0);
      }
      if( options->getTestForLinkageWithAncestry() ){
	UpdateScoreForLinkage( ind, i, TT, 1, &LocusLinkageScore,
			       &LocusLinkageInfo, (*ExpY0)(i,0),lambda0);
      }
    }
  }

  //the rest formerly was the function UpdateScoreStats
  int kk;
  Matrix_d score, info;
  for( int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ )
    if( options->getTestForMisspecifiedAlleleFreqs() )
      (*Lociptr)(j)->SumScoreForMisSpecOfAlleleFreqs();

  if( options->getScoreTestIndicator() ){
    SumU += U;
    SumI += I;
    for( int k = 0; k < options->getPopulations(); k++ )
      for( int kk = 0; kk < individuals->getTargetSize(); kk++ )
	SumU2( k, kk ) += U( k, kk ) * U( k, kk );
  }

  if( options->getAnalysisTypeIndicator() != 1 ){
    // accumulate sums for tests for allelic association
    if( options->getTestForAllelicAssociation() )
      SumScoreForWithinHaplotypeAssociation();
    for( int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
      if( (options->getTestForAllelicAssociation()
	   && (*Lociptr)(j)->GetNumberOfLoci() == 1)
	  || options->getTestForSNPsInHaplotype() ){
	// kk is set to 1 if single diallelic locus, to number of
	// merged haplotypes if compound locus, to number of states
	// if >2 alleles
	if( (*Lociptr)(j)->GetNumberOfStates() == 2 )
	  kk = 1;
	else if( (*Lociptr)(j)->GetNumberOfLoci() > 1 )
	  kk = (*Lociptr)(j)->GetNumberOfMergedHaplotypes();
	else
	  kk = (*Lociptr)(j)->GetNumberOfStates();
	// algorithm corrects for posterior covariance between regression parameters
	TransformScoreStatistics( kk, LocusLinkageAlleleScore(j),
				  LocusLinkageAlleleInfo(j),
				  &score, &info );
	SumLocusLinkageAlleleScore(j) += score;
	SumLocusLinkageAlleleScore2(j) += score * score.Transpose();
	SumLocusLinkageAlleleInfo(j) += info;
      }
        
      // accumulate sums for tests for linkage with ancestry
      if( options->getTestForLinkageWithAncestry() ){
	TransformScoreStatistics( options->getPopulations(),
				  LocusLinkageScore(j), LocusLinkageInfo(j),
				  &score, &info );
	SumLocusLinkageScore(j) += score;
	for( int kk = 0; kk < options->getPopulations(); kk++ )
	  SumLocusLinkageScore2(j)( kk, 0 ) += score( kk, 0 ) * score( kk, 0 );
	SumLocusLinkageInfo(j) += info.GetDiagonal().ColumnMatrix();
      }
      if( options->getTestForAffectedsOnly() ){
	SumAffectedsScore(j) += AffectedsScore(j);
	SumAffectedsVarScore(j) += AffectedsVarScore(j);
	SumAffectedsInfo(j) += AffectedsInfo(j);
	for( int kk = 0; kk < AffectedsScore(j).GetNumberOfRows(); kk++ )
	  SumAffectedsScore2(j)( kk, 0 ) +=
	    AffectedsScore(j)( kk, 0 ) * AffectedsScore(j)( kk, 0 );
      }
    }
  }
}

void ScoreTests::TransformScoreStatistics( int kk, Matrix_d score, Matrix_d info,
					   Matrix_d *newscore, Matrix_d *newinfo )
{
  Matrix_d score1, score2, Vbb, Vab, Vaa;
  score1 = score.SubMatrix( 0, kk - 1, 0, 0 );
  score2 = score.SubMatrix( kk, kk + options->getPopulations() - 1, 0, 0 );
  Vaa = info.SubMatrix( 0, kk - 1, 0, kk - 1 );
  Vbb = info.SubMatrix( kk, kk + options->getPopulations() - 1, kk, kk + options->getPopulations() - 1 );
  Vab = info.SubMatrix( 0, kk - 1, kk, kk + options->getPopulations() - 1 );
  Vbb.InvertUsingLUDecomposition();
  *newscore = score1 - Vab * Vbb * score2;
  *newinfo = Vaa - Vab * Vbb * Vab.Transpose();
}

void
ScoreTests::UpdateScoreForLinkageAffectedsOnly( Individual* ind)
{
  // Different from the notation in McKeigue et  al. (2000). McKeigue
  // uses P0, P1, P2, which relate to individual admixture as follows;
  // P0 = ( 1 - theta0 ) * ( 1 - theta1 )
  // P1 = ( 1 - theta0 ) * theta1 + theta0 * ( 1 - theta1 )
  // P2 = theta0 * theta1

  // Test in McKeigue et al. (2000) was on r defined on the positive
  // real line.  This test is on log r, which should have better
  // assymptotic properties.
  int locus = 0;
  Vector_d M(2), EXX;
  Matrix_d temp;
  MatrixArray_d EX( Lociptr->GetNumberOfCompositeLoci() ),
    Cov( Lociptr->GetNumberOfCompositeLoci(), options->getPopulations(), 1 );
  for( int j = 0; j < chrm->size(); j++ ){
    for( int jj = 0; jj < (*chrm)(j)->GetSize(); jj++ ){
      //         temp = (*chrm)(j)->GetExpectationofAncestry(jj);
      temp = ind->getExpectedAncestry(locus);
      EX(locus) = temp.SubMatrix( 0, options->getPopulations() - 1, 0 , 1 );
      EXX = temp.GetColumn(2);
      for( int k = 0; k < options->getPopulations(); k++ ){
	Cov(locus)(k,0) = EXX(k) - EX(locus)(k,0)*EX(locus)(k,1);
      }
      locus++;
    }
  }
   
  for( int k = 0; k < options->getPopulations(); k++ ){
    locus = 0;
    M(0) = ind->getAncestry()( k, 0 );
    if( options->getModelIndicator() )
      M(1) = ind->getAncestry()( k, 1 );
    else
      M(1) = ind->getAncestry()( k, 0 );
      
    for( int j = 0; j < chrm->size(); j++ ){
      for( int jj = 0; jj < (*chrm)(j)->GetSize(); jj++ ){
	for( int g = 0; g < 2; g++ ){
	  AffectedsScore( locus )( k, 0 ) += 0.5*( EX(locus)(k,g) - M(g) );
	  AffectedsInfo( locus )( k, 0 ) += 0.25 * ( M(g) - M(g)*M(g) );
	  AffectedsVarScore( locus )( k, 0 ) += 0.25*EX(locus)(k,g)*( 1 - EX(locus)(k,g) );
	}
	AffectedsVarScore( locus )( k, 0 ) += 0.5*Cov(locus)(k,0);
	locus++;
	//          cout << M(0) << " " << M(1) << endl << xx << endl;
      }
    }
    //     for( int j = 0; j < Loci.GetNumberOfCompositeLoci(); j++ ){
    //        double x = 0.5 * ( GeneAncestry(i)( k, 2 * j ) + GeneAncestry(i)( k, 2 * j + 1 ) );
    //        AffectedsScore( j )( k, 0 ) += x - 0.5*( M(0) + M(1) );
    //        AffectedsInfo( j )( k, 0 ) += 0.25 * ( M(0) - M(0)*M(0) + M(1) - M(1)*M(1) );
    //     }
  }
}

// p (= alpha) in methods below should be calculated once and stored

void ScoreTests::UpdateScoreForLinkage( Individual* ind, int i, int regressonindicator, int hypothesisindicator, MatrixArray_d *score, MatrixArray_d *info,double EY0, double lambda0)//EY0 = ExpectedY(0)(i,0), lambda0 = lambda(0)
{
  // regressonindicator = 0 => linear regression model
  // regressonindicator = 1 => logistic regression model
  // hypothesisindicator = 0 => test for allelic association
  // hypothesisindicator = 1 => test for association with ancestry

  int dim = 0; //? this removes warning, but isn't necessarily the Right Thing
  int locus = 0;
  Vector_i hap, ancestry, AlleleCounts;
  Vector_d xx;

  double IndividualTargetMinusP = individuals->getTarget(0)( i, 0 ) - EY0;

  for( int j = 0; j < chrm->size(); j++ ){
    for( int jj = 0; jj < (*chrm)(j)->GetSize(); jj++ ){
      if( hypothesisindicator == 0 ){
	if( (*Lociptr)(locus)->GetNumberOfStates() == 2 ){
	  dim = options->getPopulations() + 1;
	} else if( (*Lociptr)(locus)->GetNumberOfLoci() == 1 ){
	  dim = (*Lociptr)(locus)->GetNumberOfStates() + options->getPopulations();
	} else {
	  dim = (*Lociptr)(locus)->GetNumberOfMergedHaplotypes() + options->getPopulations();
	}
      } else if( hypothesisindicator == 1 ) {
	dim = 2 * options->getPopulations();
      }
        
      // Set x-co-ordinates of covariates in model
      Matrix_d cov_x_coord(dim,1);
      cov_x_coord( dim - 1, 0 ) = 1;
      for( int k = 0; k < options->getPopulations() - 1; k++ ){
	cov_x_coord( k + dim - options->getPopulations(), 0 ) = ind->getAncestry()( k, 0 );
      }
     
      // Set x co-ordinate for regression parameter under test
      // Test for allelic assocition
      if( hypothesisindicator == 0 ){
	// special case for SNP
	if( (*Lociptr)(locus)->GetNumberOfStates() == 2 ){
	  //? next line relies on encoding system
	  cov_x_coord( 0, 0 ) = (int)(ind->getGenotype(locus)[0]) + (int)(ind->getGenotype(locus)[1] - 3.0);
	  // general case for simple locus
	} else if((*Lociptr)(locus)->GetNumberOfLoci() == 1 ){
	  for( int k = 0; k < (*Lociptr)(locus)->GetNumberOfStates(); k++ ){
	    if( (int)(ind->getGenotype(locus)[0]) == k + 1 )
	      cov_x_coord( k, 0 )++;
	    if( (int)(ind->getGenotype(locus)[1]) == k + 1 )
	      cov_x_coord( k, 0 )++;
	  }
	  // general case for composite locus
	} else {
	  UpdateScoreForWithinHaplotypeAssociation( i, locus, regressonindicator, EY0,lambda0 );
	  if( options->getTestForSNPsInHaplotype() ){
	    ancestry = ind->GetLocusAncestry( j, jj );
	    hap = (*Lociptr)(locus)->SampleHaplotype( ind->getGenotype(locus), ancestry );
	    for( int k = 0; k <(*Lociptr)(locus)->GetNumberOfStates(); k++ ){
	      if( hap(0) == k )
		cov_x_coord( (*Lociptr)(locus)->GetMergedHaplotype(k), 0 )++;
	      if( hap(1) == k )
		cov_x_coord((*Lociptr)(locus)->GetMergedHaplotype(k), 0 )++;
	    }
	  }
	}
      } else if( hypothesisindicator == 1 ){
	ancestry = ind->GetLocusAncestry( j, jj );
	cov_x_coord( ancestry(0), 0 )++;
	cov_x_coord( ancestry(1), 0 )++;
	//           xx = (*chrm)(j)->GetExpectationofAncestry(jj);
	//           for( int k = 0; k < options->getPopulations(); k++ ){
	//              cov_x_coord(k,0) = xx(k);
	//           }
	//          cout << ancestry << endl;
	//          cout << xx << endl;
      }

      // Check if missing genotype
      //? NB: zero used for missing genotype here
      if( !(hypothesisindicator == 0 && ind->IsMissing(locus)[0] == 0 ) ){
	if( regressonindicator ){
	  (*score)( locus ) += cov_x_coord * IndividualTargetMinusP;
	} else {
	  (*score)( locus ) += (lambda0 * IndividualTargetMinusP) * cov_x_coord;
	}
	if( regressonindicator ){
	  (*info)( locus ) += (cov_x_coord * cov_x_coord.Transpose()) * (EY0 * ( 1 - EY0 ));
	} else {
	  (*info)( locus ) += (cov_x_coord * cov_x_coord.Transpose()) * lambda0;
	}
      }
      locus++;
    }
  }
}

void ScoreTests::UpdateScoreForAssociation_Binary( int i, Matrix_d Theta, Matrix_d  *score, Matrix_d  *info , double EY0)
{
  double x, p;
  p = EY0;
  for( int k = 0; k < options->getPopulations(); k++ ){
    if( options->getModelIndicator() )
      x = 0.5 * ( Theta( k, 0 ) + Theta( k, 1 ) );
    else
      x = Theta( k, 0 );
    (*score)( k, 0 ) += x * ( individuals->getTarget(0)( i, 0 ) - p );
    (*info)( k, 0 ) += x * x * p * ( 1 - p );
  }
}

void ScoreTests::UpdateScoreForAssociation_Continuous( int i, Matrix_d Theta, Matrix_d  *score, Matrix_d  *info, double EY0,double lambda0)
{
  double x, alpha;
  alpha = EY0;
  for( int k = 0; k < options->getPopulations(); k++ ){
    if( options->getModelIndicator() )
      x = 0.5 * ( Theta( k, 0 ) + Theta( k, 1 ) );
    else
      x = Theta( k, 0 );
    (*score)( k, 0 ) += lambda0 * x * ( individuals->getTarget(0)( i, 0 ) - alpha );
    (*info)( k, 0 ) += lambda0 * x * x;
  }
}

void ScoreTests::UpdateScoresForMisSpecOfAlleleFreqs( int i)
{
  Individual* ind = individuals->getIndividual(i);
  Matrix_d phi( options->getPopulations(), options->getPopulations() );

  for( int k = 0; k < options->getPopulations(); k++ )
    for( int kk = 0; kk < options->getPopulations(); kk++ )
      phi( k, kk ) = ind->getAncestry()( k, 0 ) * ind->getAncestry()( kk, 0 );
   
  for( int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ )
    //? relies on 'missing' being encoded as zero
    if( individuals->getIndividual(i)->IsMissing(j)[0] && (*Lociptr)(j)->GetNumberOfLoci() == 1  && !((*Lociptr)(j)->IsRandom()) )
      (*Lociptr)(j)->UpdateScoreForMisSpecOfAlleleFreqs( phi, individuals->getIndividual(i)->getGenotype(j) );
}

// This method calculates score for allelic association at each simple locus within a compound locus
void ScoreTests::UpdateScoreForWithinHaplotypeAssociation( int i, int j, int indicator, double p,double lambda0)
{
  Individual* ind = individuals->getIndividual(i);
  Vector_i AlleleCounts;
  Vector_d x( options->getPopulations() + 1 );
  AlleleCounts = (*Lociptr)(j)->GetAlleleCountsInHaplotype(ind->getGenotype(j));
  Matrix_d info( options->getPopulations() + 1, options->getPopulations() + 1 );

  x( options->getPopulations() ) = 1;
  for( int k = 0; k < options->getPopulations() - 1; k++ )
    x( k + 1 ) = ind->getAncestry()( k, 0 );

  for( int l = 0; l < (*Lociptr)(j)->GetNumberOfLoci(); l++ ){
    x(0) = AlleleCounts( l );
    for( int k = 0; k < options->getPopulations() + 1; k++ ){
      if( indicator && x(0) != 99 ){
	ScoreWithinHaplotype[j](l)( k, 0 ) += x(k) * ( individuals->getTarget(0)( i, 0 ) - p );
      }
      else if( x(0) != 99 ){
	ScoreWithinHaplotype[j](l)( k, 0 ) += lambda0 * x(k) * ( individuals->getTarget(0)( i, 0 ) - p );
      }
      for( int kk = 0; kk < options->getPopulations() + 1; kk++ )
	info( k, kk ) = x( k ) * x( kk );
    }
    if( indicator && x(0) != 99 )
      InfoWithinHaplotype[j](l) += info * p * ( 1 - p );
    else if( x(0) != 99 )
      InfoWithinHaplotype[j](l) += info * lambda0;
  }
}

void ScoreTests::SumScoreForWithinHaplotypeAssociation()
{
  int kk = 1;
  Matrix_d score1, score2, info, Vbb, Vab, Vaa;

  for( int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
    if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){
      for( int l = 0; l < (*Lociptr)(j)->GetNumberOfLoci(); l++ ){
	score1 = ScoreWithinHaplotype[j](l).SubMatrix( 0, kk - 1, 0, 0 );
	score2 = ScoreWithinHaplotype[j](l).SubMatrix( kk, kk + options->getPopulations() - 1, 0, 0 );
	Vaa = InfoWithinHaplotype[j](l).SubMatrix( 0, kk - 1, 0, kk - 1 );
	Vbb = InfoWithinHaplotype[j](l).SubMatrix( kk, kk + options->getPopulations() - 1,
						   kk, kk + options->getPopulations() - 1 );
	Vab = InfoWithinHaplotype[j](l).SubMatrix( 0, kk - 1, kk, kk + options->getPopulations() - 1 );
	Vbb.InvertUsingLUDecomposition();
	score1 = score1 - Vab * Vbb * score2;
	info = Vaa - Vab * Vbb * Vab.Transpose();
	SumScoreWithinHaplotype( j )( l, 0 ) += score1(0,0);
	SumScore2WithinHaplotype( j )( l, 0 ) += score1(0,0) * score1(0,0);
	SumInfoWithinHaplotype( j )( l, 0 ) += info(0,0);
      }
    }
  }
}

void ScoreTests::Output(int iteration,std::string * PLabels){
  PopLabels = PLabels;
  //Allelic association
  if( options->getTestForAllelicAssociation() ){
    OutputTestsForAllelicAssociation( iteration );
  }
  //haplotype association
  if( strlen( options->getTestsForSNPsInHaplotypeOutputFilename() ) ){
    OutputTestsForSNPsInHaplotype( iteration );
  }
  //ancestry association
  if( options->getTestForLinkageWithAncestry() ){
    OutputTestsForLocusLinkage( iteration, ancestryAssociationScoreStream,
				SumLocusLinkageScore,
				SumLocusLinkageScore2, SumLocusLinkageInfo );
  }
  //affectedonly
  if( options->getTestForAffectedsOnly() ){
    OutputTestsForLocusLinkage2( iteration, affectedsOnlyScoreStream,
				 SumAffectedsScore, SumAffectedsVarScore,
				 SumAffectedsScore2, SumAffectedsInfo );
	   
  }
  //misspecified allele freqs
  if( options->getTestForMisspecifiedAlleleFreqs() ){
    OutputTestsForMisSpecifiedAlleleFreqs( iteration);
  }
  if( options->getTestForMisspecifiedAlleleFreqs2() ){
    OutputTestsForMisSpecifiedAlleleFreqs2( iteration);
  }
  //????
  if( options->getScoreTestIndicator() ){
    OutputScoreTest( iteration );
  }

}//end of output definition

// next few methods calculate score tests from the cumulative sums of
// the score, score squared, and information score and info can be
// scalars, or respectively a vector and a matrix should have just one
// method or class to do this
void ScoreTests::OutputTestsForSNPsInHaplotype( int iteration )
{
  int NumberOfMergedHaplotypes;
  Vector_i hap;
  Matrix_d ScoreMatrix, CompleteMatrix, ObservedMatrix;
  for( int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
    if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){
      ScoreMatrix = SumLocusLinkageAlleleScore(j) / ( iteration - options->getBurnIn() );
      CompleteMatrix =  SumLocusLinkageAlleleInfo(j) / ( iteration - options->getBurnIn() );
      ObservedMatrix = CompleteMatrix + ScoreMatrix * ScoreMatrix.Transpose() - SumLocusLinkageAlleleScore2(j) / ( iteration - options->getBurnIn() );
      NumberOfMergedHaplotypes = ScoreMatrix.GetNumberOfRows();
      for( int k = 0; k < NumberOfMergedHaplotypes; k++ ){
	if(options->IsPedFile())
	  *SNPsAssociationScoreStream  << "\"" << (*Lociptr)(j)->GetLabel(0) << "\"" << ",";
	else
	  *SNPsAssociationScoreStream  << (*Lociptr)(j)->GetLabel(0) << ",";
	if( k < NumberOfMergedHaplotypes - 1 ){
	  hap = (*Lociptr)(j)->GetHapLabels(k);
	  *SNPsAssociationScoreStream  << "\"";
	  for( int kk = 0; kk < (*Lociptr)(j)->GetNumberOfLoci() - 1; kk++ ){
	    *SNPsAssociationScoreStream  << hap(kk) << "-";
	  }
	  *SNPsAssociationScoreStream  << hap((*Lociptr)(j)->GetNumberOfLoci() - 1) << "\",";
	}
	else
	  *SNPsAssociationScoreStream  << "\"others\",";
	*SNPsAssociationScoreStream  << double2R(ScoreMatrix(k,0)) << ",";
	*SNPsAssociationScoreStream  << double2R(CompleteMatrix(k,k)) << ",";
	*SNPsAssociationScoreStream  << double2R(ObservedMatrix(k,k)) << ",";
	*SNPsAssociationScoreStream  << double2R(100*ObservedMatrix(k,k) / CompleteMatrix(k,k)) << ",";
	*SNPsAssociationScoreStream  << double2R(ScoreMatrix( k, 0 ) / sqrt( ObservedMatrix( k, k ) ));
	*SNPsAssociationScoreStream  << ",";
	// if not last allele at locus, output NA in chi-square column
	if( k != NumberOfMergedHaplotypes - 1 ){
	  *SNPsAssociationScoreStream  << "\"NA\"," << endl;
	}
      }
      // calculate summary chi-square statistic
      ObservedMatrix = ObservedMatrix.SubMatrix( 0, NumberOfMergedHaplotypes - 2,
						 0, NumberOfMergedHaplotypes - 2 );
      ScoreMatrix = ScoreMatrix.SubMatrix( 0, NumberOfMergedHaplotypes - 2, 0, 0 );
      ObservedMatrix.InvertUsingLUDecomposition();
      *SNPsAssociationScoreStream  << double2R((ScoreMatrix.Transpose() * ObservedMatrix * ScoreMatrix)(0,0)) << "," << endl;
    }
  }
}

void ScoreTests::OutputScoreTest(int iteration)
{
  for( int j = 0; j < options->getPopulations(); j++ ){
    for( int jj = 0; jj < individuals->getTargetSize(); jj++ ){
      double EU = SumU( j, jj ) / ( iteration - options->getBurnIn() );
      double complete = SumI( j, jj ) / ( iteration - options->getBurnIn() );
      double missing = SumU2( j, jj ) / ( iteration - options->getBurnIn() ) - EU * EU;
      assocscorestream.width(9);
      assocscorestream << setprecision(6) << double2R(complete) << " ";
      assocscorestream.width(9);
      assocscorestream << setprecision(6) << double2R(missing) << " ";
      assocscorestream.width(9);
      assocscorestream << setprecision(6) << double2R(EU / sqrt( complete - missing )) << " ";
    }
  }
  assocscorestream << endl;
}

void ScoreTests::OutputTestsForAllelicAssociation( int iteration )
{
  double Score, CompleteInfo, MissingInfo, ObservedInfo;
  for( int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
    for( int l = 0; l < (*Lociptr)(j)->GetNumberOfLoci(); l++ ){
      if((* Lociptr)(j)->GetNumberOfLoci() == 1 ){
	Score = SumLocusLinkageAlleleScore(j)(0,0) / ( iteration - options->getBurnIn() );
	CompleteInfo =  SumLocusLinkageAlleleInfo(j)(0,0) / ( iteration - options->getBurnIn() );
	MissingInfo = SumLocusLinkageAlleleScore2(j)(0,0) / ( iteration - options->getBurnIn() ) - Score * Score;
	ObservedInfo = CompleteInfo - MissingInfo;
      }
      else{
	Score = SumScoreWithinHaplotype( j )( l, 0 ) / ( iteration - options->getBurnIn() );
	CompleteInfo = SumInfoWithinHaplotype( j )( l, 0 ) / ( iteration - options->getBurnIn() );
	MissingInfo = SumScore2WithinHaplotype( j )( l, 0 ) / ( iteration - options->getBurnIn() ) - Score * Score;
	ObservedInfo = CompleteInfo - MissingInfo;
      }
      if(options->IsPedFile())
	genescorestream << "\"" << (*Lociptr)(j)->GetLabel(l) << "\"" << ",";
      else
	genescorestream << (*Lociptr)(j)->GetLabel(l) << ",";
      genescorestream << double2R(Score)        << ",";
      genescorestream << double2R(CompleteInfo) << ",";
      genescorestream << double2R(ObservedInfo) << ",";
      genescorestream << double2R(100*ObservedInfo / CompleteInfo) << ",";
      genescorestream << double2R(Score / sqrt( ObservedInfo ))    << "," << endl;
    }
  }
}

void ScoreTests::OutputTestsForMisSpecifiedAlleleFreqs2( int iteration)
{
  int samples = iteration - options->getBurnIn();
  Matrix_d score, completeinfo, observedinfo;
  for( int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
    CompositeLocus *locus = (CompositeLocus*)(*Lociptr)(j);
    for( int k = 0; k < options->getPopulations(); k++ ){
      score = locus->GetNewScore(k) / samples;
      completeinfo = locus->GetNewInfo(k) / samples;
      observedinfo = completeinfo + score * score.Transpose() - locus->GetNewScoreSq(k) / samples;
      if(options->IsPedFile())
	allelefreqscorestream2 << "\"" << (*Lociptr)(j)->GetLabel(0) << "\"" << ",";
      else
	allelefreqscorestream2 << (*Lociptr)(j)->GetLabel(0) << ",";
      allelefreqscorestream2 << PopLabels[k] << ",";
      allelefreqscorestream2 << double2R(completeinfo.Determinant()) << ",";
      allelefreqscorestream2 << double2R(observedinfo.Determinant()) << ",";
      allelefreqscorestream2 << double2R(100*observedinfo.Determinant() / completeinfo.Determinant()) << ",";
      observedinfo.InvertUsingLUDecomposition();
      allelefreqscorestream2 << double2R((score.Transpose() * observedinfo * score)(0,0)) << "," << endl;
    }
  }
}

void ScoreTests::OutputTestsForMisSpecifiedAlleleFreqs( int iteration)
{
  int samples = iteration - options->getBurnIn();
  Matrix_d ScoreMatrix, CompleteMatrix, ObservedMatrix;
  for( int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
    if( (*Lociptr)(j)->GetNumberOfLoci() == 1 ){
      ScoreMatrix = (*Lociptr)(j)->GetScore() / samples;
      CompleteMatrix = (*Lociptr)(j)->GetInfo() / samples;
      ObservedMatrix = CompleteMatrix + ScoreMatrix * ScoreMatrix.Transpose() - (*Lociptr)(j)->GetScoreSq() / samples;
      for( int k = 0; k < options->getPopulations(); k++ ){
	// Test for mis-specification within each continental-population.
	if(options->IsPedFile())
	  allelefreqscorestream << "\"" << (*Lociptr)(j)->GetLabel(0) << "\"" << ",";
	else
	  allelefreqscorestream << (*Lociptr)(j)->GetLabel(0) << ",";
	allelefreqscorestream << PopLabels[k] << ",";
	allelefreqscorestream << double2R(ScoreMatrix( k, 0 ) ) << ",";
	allelefreqscorestream << double2R(CompleteMatrix( k, k ) ) << ",";
	allelefreqscorestream << double2R(ObservedMatrix( k, k ) ) << ",";
	allelefreqscorestream << double2R(100*ObservedMatrix( k, k ) / CompleteMatrix( k, k ) ) << ",";
	allelefreqscorestream << double2R(ScoreMatrix( k, 0 ) / sqrt( ObservedMatrix( k, k ) ) ) << ",";
	if( k < options->getPopulations() - 1 )
	  allelefreqscorestream  << "\"NA\"," << endl;
      }
      // Summary chi-sq test for mis-specification in all continental-populations.
      ObservedMatrix = ObservedMatrix.SubMatrix( 0, options->getPopulations() - 1, 0, options->getPopulations() - 1 );
      ScoreMatrix = ScoreMatrix.SubMatrix( 0, options->getPopulations() - 1, 0, 0 );
      ObservedMatrix.InvertUsingLUDecomposition();
      allelefreqscorestream << (ScoreMatrix.Transpose() * ObservedMatrix * ScoreMatrix)(0,0) << "," << endl;
    }
  }
}

void
ScoreTests::OutputTestsForLocusLinkage( int iteration, ofstream* outputstream,
					MatrixArray_d Score, MatrixArray_d Score2,
					MatrixArray_d Info )
{
  double EU, missing, complete;
  for( int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
    for( int k = 0; k < options->getPopulations(); k++ ){
      if(options->IsPedFile())
	*outputstream << "\"" << (*Lociptr)(j)->GetLabel(0) << "\"" << ",";
      else
	*outputstream << (*Lociptr)(j)->GetLabel(0) << ",";
      *outputstream << PopLabels[k] << ",";
      
      EU = Score( j )( k, 0 ) / ( iteration - options->getBurnIn() );
      missing = Score2( j )( k, 0 ) / ( iteration - options->getBurnIn() ) - EU * EU;
      complete =  Info( j )( k, 0 ) / ( iteration - options->getBurnIn() );
      
      *outputstream << double2R(EU)                                << ",";
      *outputstream << double2R(complete)                          << ",";
      *outputstream << double2R(complete - missing)                << ",";
      *outputstream << double2R(100*(complete - missing)/complete) << ",";
      *outputstream << double2R(EU / sqrt( complete - missing ))   << "," << endl;
    }
  }
}

void
ScoreTests::OutputTestsForLocusLinkage2( int iteration, ofstream* outputstream,
					 MatrixArray_d Score, MatrixArray_d VarScore,
					 MatrixArray_d Score2, MatrixArray_d Info )
{
  double VU, EU, missing, complete;
  for( int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
    for( int k = 0; k < options->getPopulations(); k++ ){
      if(options->IsPedFile())
	*outputstream << "\"" << (*Lociptr)(j)->GetLabel(0) << "\"" << ",";
      else
	*outputstream << (*Lociptr)(j)->GetLabel(0) << ",";
      *outputstream << PopLabels[k] << ",";
      
      EU = Score( j )( k, 0 ) / ( iteration - options->getBurnIn() );
      VU = VarScore( j )( k, 0 ) / ( iteration - options->getBurnIn() );
      missing = Score2( j )( k, 0 ) / ( iteration - options->getBurnIn() ) - EU * EU + VU;
      complete =  Info( j )( k, 0 ) / ( iteration - options->getBurnIn() );
      
      *outputstream << double2R(EU)                                << ",";
      *outputstream << double2R(complete)                          << ",";
      *outputstream << double2R(complete - missing)                << ",";
      *outputstream << double2R(100*(complete - missing)/complete) << ",";
      *outputstream << double2R(100*(VU/complete))                 << ",";//%missing info attributable to locus ancestry
      *outputstream << double2R(100*(missing-VU)/complete)         << ",";//%missing info attributable to       
      *outputstream << double2R(EU / sqrt( complete - missing ))   << "," << endl;
    }
  }
}
void ScoreTests::ROutput(){
  int numPrintedIterations = options->getTotalSamples()/ options->getSampleEvery() / 10  -  options->getBurnIn() / options->getSampleEvery() / 10;

  /**
   * writes out the dimensions and labels of the 
   * R-matrix previously written to genescorestream
   */
  int count;
  if(options->getTestForAllelicAssociation()){
    count = 0;
    for( int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
      count += (*Lociptr)(j)->GetNumberOfLoci();
    }         
    vector<int> dimensions(3,0);
    dimensions[0] = 6;
    dimensions[1] = count;
    dimensions[2] = (int)(numPrintedIterations);
    vector<string> labels(dimensions[0],"");
    labels[0] = "Locus";
    labels[1] = "Score";
    labels[2] = "CompleteInfo";
    labels[3] = "ObserevedInfo";
    labels[4] = "PercentInfo";
    labels[5] = "StdNormal";
    R_output3DarrayDimensions(&genescorestream,dimensions,labels);
  }

  /** 
   * writes out the dimensions and labels of the        
   * R-matrix for score tests of SNPs in haplotypes.        
   */ 
  if( strlen( options->getTestsForSNPsInHaplotypeOutputFilename() ) ){      
    count = 0;
    for( int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
      if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){
	count += (*Lociptr)(j)->GetNumberOfMergedHaplotypes();
      }
    }         
    vector<int> dimensions(3,0);
    dimensions[0] = 8;
    dimensions[1] = count;
    dimensions[2] = (int)(numPrintedIterations);
    vector<string> labels(dimensions[0],"");
    labels[0] = "Locus";
    labels[1] = "Haplotype";
    labels[2] = "Score";
    labels[3] = "CompleteInfo";
    labels[4] = "ObserevedInfo";
    labels[5] = "PercentInfo";
    labels[6] = "StdNormal";
    labels[7] = "ChiSquare";
    R_output3DarrayDimensions(SNPsAssociationScoreStream,dimensions,labels);
  }
   
  /**
   * writes out the dimensions and labels of the 
   * R-matrix for old test for mis-specified allele frequencies
   */
  if( options->getTestForMisspecifiedAlleleFreqs() ){
    vector<int> dimensions(3,0);
    dimensions[0] = 8;
    dimensions[1] = Lociptr->GetNumberOfCompositeLoci() * options->getPopulations();
    dimensions[2] = (int)(numPrintedIterations);

    vector<string> labels(8,"");
    labels[0] = "Locus";
    labels[1] = "Population";
    labels[2] = "Score";
    labels[3] = "CompleteInfo";
    labels[4] = "ObservedInfo";
    labels[5] = "PercentInfo";
    labels[6] = "StdNormal";
    labels[7] = "ChiSquared";

    R_output3DarrayDimensions(&allelefreqscorestream,dimensions,labels);
  }

  /**
   * writes out the dimensions and labels of the 
   * R-matrix for new test for mis-specified allele frequencies
   */
  if( options->getTestForMisspecifiedAlleleFreqs2() ){
    vector<int> dimensions(3,0);
    dimensions[0] = 6;
    dimensions[1] = Lociptr->GetNumberOfCompositeLoci() * options->getPopulations();
    dimensions[2] = (int)(numPrintedIterations);

    vector<string> labels(6,"");
    labels[0] = "Locus";
    labels[1] = "Population";
    labels[2] = "CompleteInfo";
    labels[3] = "ObservedInfo";
    labels[4] = "PercentInfo";
    labels[5] = "ChiSquared";

    R_output3DarrayDimensions(&allelefreqscorestream2,dimensions,labels);
  }

  /**
   * writes out the dimensions and labels of the 
   * R-matrix previously written to ancestryAssociationScoreStream
   */
  if (options->getTestForLinkageWithAncestry()){
    vector<int> dimensions(3,0);
    dimensions[0] = 7;
    dimensions[1] = Lociptr->GetNumberOfCompositeLoci() * options->getPopulations();
    dimensions[2] = (int)(options->getTotalSamples()/ options->getSampleEvery() / 10  -  options->getBurnIn() / options->getSampleEvery() / 10);

    vector<string> labels(7,"");
    labels[0] = "Locus";
    labels[1] = "Population";
    labels[2] = "Score";
    labels[3] = "CompleteInfo";
    labels[4] = "ObservedInfo";
    labels[5] = "PercentInfo";
    labels[6] = "StdNormal";

    R_output3DarrayDimensions(ancestryAssociationScoreStream,dimensions,labels);
  }

  /**
   * writes out the dimensions and labels of the 
   * R-matrix previously written to affectedsOnlyScoreStream
   */
  if (options->getTestForAffectedsOnly()){
    vector<int> dimensions(3,0);
    dimensions[0] = 9;//was 7
    dimensions[1] = Lociptr->GetNumberOfCompositeLoci() * options->getPopulations();
    dimensions[2] = (int)(numPrintedIterations);
     
    vector<string> labels(9,"");//was 7
    labels[0] = "Locus";
    labels[1] = "Population";
    labels[2] = "Score";
    labels[3] = "CompleteInfo";
    labels[4] = "ObservedInfo";
    labels[5] = "PercentInfo";
    labels[6] = "Missing1";//Percent of missing info due to locus ancestry
    labels[7] = "Missing2";//Percent of missing info due to     
    labels[8] = "StdNormal";

    R_output3DarrayDimensions(affectedsOnlyScoreStream,dimensions,labels);
  }
}

void
ScoreTests::R_output3DarrayDimensions(ofstream* stream,vector<int> dim,vector<string> labels)
{
  *stream << ")," << endl;
  *stream << ".Dim = c(";
  for(unsigned int i=0;i<dim.size();i++){
    *stream << dim[i];
    if(i != dim.size() - 1){
      *stream << ",";
    }
  }
  *stream << ")," << endl;
  *stream << ".Dimnames=list(c(";
  for(unsigned int i=0;i<labels.size();i++){
    *stream << "\"" << labels[i] << "\"";
    if(i != labels.size() - 1){
      *stream << ",";
    }
  }
  *stream << "), character(0), character(0)))" << endl;


}
string
ScoreTests::double2R( double x )
{
  if( isnan(x) )
    return "NaN";
  else{
    stringstream ret;
    ret << x;
    return( ret.str() );
  }
}
