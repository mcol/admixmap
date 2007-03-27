#!/usr/bin/perl -w

#This creates the input files for the ADMIXMAP program using a snp, badsnp, indiv and genotype files
#snpfile: snpcnts file, badsnpfile: list of markers that are to be deleted, genofile: genotype file in ANCESTRYMAP format
#phenofile: phenotype file in ANCESTRYMAP format. Then it outputs the files in ADMIXMAP format
#genooutfile: ADMIXMAP genotype file, indoutfile: ADMIXMAP individual file, historicfreqoutfile: ADMIXMAP historic allele
#frequency file and locicoutfile: ADMIXMAP loci file

$snpfile = $ARGV[0];
$badsnpfile = $ARGV[1];
$genofile = $ARGV[2];
$phenofile = $ARGV[3];
$genooutfile = $ARGV[4];
$indoutfile = $ARGV[5];
$historicfreqoutfile = $ARGV[6];
$locioutfile = $ARGV[7];

unless(open(SNPFILE ,"$snpfile"))           
  {                                              
                                                                      
    print "could not open file $snpfile!\n";    
    exit;    
  }          

unless(open(BADSNPFILE ,"$badsnpfile"))           
  {                                              
                                                                      
    print "could not open file $badsnpfile!\n";    
    exit;    
  }          

unless(open(GENOFILE ,"$genofile"))           
  {                                              
                                                                      
    print "could not open file $genofile!\n";    
    exit;    
  }      
unless(open(PHENOFILE ,"$phenofile"))           
  {                                              
                                                                      
    print "could not open file $phenofile!\n";    
    exit;    
  }              
unless(open(OUTFILE ,">$genooutfile"))           
  {                                              
                                                                      
    print "could not open file $genooutfile!\n";    
    exit;    
  }          
unless(open(INDOUTFILE ,">$indoutfile"))           
  {                                              
                                                                      
    print "could not open file $indoutfile!\n";    
    exit;    
  }          
unless(open(FREQOUTFILE ,">$historicfreqoutfile"))           
  {                                              
                                                                      
    print "could not open file $historicfreqoutfile!\n";    
    exit;    
  }          
unless(open(LOCIOUTFILE ,">$locioutfile"))           
  {                                              
                                                                      
    print "could not open file $locioutfile!\n";    
    exit;    
  }          

#Extact the snp information from badsnp file and remove those from the snp hash
@row = <BADSNPFILE>;
for($i = 0; $i < scalar @row; $i++)
  {
    chomp $row[$i];
    @line = split(' ',$row[$i]);
    $badsnp{$line[0]} = $line[0];
  }

#Extact the snp information from snpcnts file
@row = <SNPFILE>;
$ncnt = 0;
for($i = 0; $i < scalar @row; $i++)
  {
    @line = split(' ',$row[$i]);
    if(!(exists $badsnp{$line[0]}))
      {
	$snp{$line[0]} = $ncnt;
	$chr{$line[0]} = $line[1];
	$genpos{$line[0]} = $line[2];
	$afrref{$line[0]} = $line[5];
	$afrvart{$line[0]} = $line[4];
	$eurref{$line[0]} = $line[7];
	$eurvart{$line[0]} = $line[6];
	$ncnt++;
      }
 
  }

print $ncnt,"\n";

#Get the phenotype information and put the phenotype and indiv_id in hashes
#if the individual is a Case or a Control: eliminate Ignore
@row = <PHENOFILE>;
$nind = 1;
for($i = 0; $i < scalar @row; $i++)
  {
    chomp $row[$i];
    @line = split(' ',$row[$i]);
    if(($line[2] =~ /Case/) || ($line[2] =~ /Control/))
      {
	$indivId{$line[0]} = $nind;
	if($line[2] =~ /Case/)
	  {
	    $pheno{$nind} = 1;
	  }
	elsif($line[2] =~ /Control/)
	  {
	    $pheno{$nind} = 0;
	  }
	$nind++;
      }

  }
print $nind,"\n";

#Extract the genotype information only for the snps and indivs that are in their respective
#hashes

while($row = <GENOFILE>)
  {
    chomp $row;
    @line = split(' ',$row);

    if((exists $snp{$line[0]}) && (exists $indivId{$line[1]}))
      {
	$geno{$snp{$line[0]},$indivId{$line[1]}} = $line[2];
      }
  }




#print OUTFILE "\"","Indiv_".$nind,"\" ";
#print "Indiv_".$nind,"\t",$indivId{$line[0]},"\n";
print INDOUTFILE "\"","MS","\" ","\n";
foreach $keyind(sort {$indivId{$a} <=> $indivId{$b}} keys %indivId)
  {
    print INDOUTFILE $pheno{$indivId{$keyind}},"\n";
  }

#Print the marker line
print OUTFILE "\"","IDs","\" ";
print FREQOUTFILE "\"","Name","\" ","\"","African","\" ","\"","European","\" ","\n"; 
foreach $keysnp(sort {$snp{$a} <=> $snp{$b}} keys %snp)
  {
    print OUTFILE "\"",$keysnp,"\" ";
    print FREQOUTFILE $keysnp,"\t",$afrref{$keysnp},"\t",$eurref{$keysnp},"\n";
    print FREQOUTFILE $keysnp,"\t",$afrvart{$keysnp},"\t",$eurvart{$keysnp},"\n";
    $numsnp++;
  }
print OUTFILE "\n";

@keymarkers = sort {$snp{$a} <=> $snp{$b}} keys %snp;
$i = 0;
$snpdata[$i] = $keymarkers[$i]; 
$snpgen[$i] = 100; 
$i++;             
 
while($i < scalar @keymarkers)                                                                    
  {     

    $snpdata[$i] = $keymarkers[$i];     
    if($chr{$snpdata[$i]} == $chr{$snpdata[$i-1]}) 
      {                       
        $snpgen[$i] = ($genpos{$snpdata[$i]} - $genpos{$snpdata[$i-1]}); 
      }                                                        
    else 
      {                                                      
        $snpgen[$i] = 100; 
      } 
    
    $i++; 
                                          
  }                                                 
                          
print LOCIOUTFILE "\"","SNP_ID","\" ","\"","Num_alleles","\" ","\"","Map_Dist","\" ","\n";                                                                     
for($j = 0; $j < scalar @snpdata; $j++) 
  {                                                  
    print LOCIOUTFILE  "\"",$snpdata[$j],"\" ","2 ",$snpgen[$j],"\n"; 
  } 

@keyindivs = sort {$indivId{$a} <=> $indivId{$b}} keys %indivId;
for($i = 1; $i < $nind; $i++)
  {
    $g = $i-1;
#    print OUTFILE "\"",$keyindivs[$g],"\" ";
    print OUTFILE "\"","Indiv_".$i,"\" ";
    for($k = 0; $k < $ncnt; $k++)
      {
	if (exists $geno{$k,$i})
	    {
	      if($geno{$k,$i} == 0) 
		{
		  print OUTFILE "\"","1,1","\" ";
		}
	      elsif($geno{$k,$i} == 1)
		{
		  print OUTFILE "\"1,2\" ";
		}
	      elsif($geno{$k,$i} == 2)
		{		
		  print OUTFILE "\"2,2\" ";
		}
	      elsif($geno{$k,$i} == -1)
		{		
		  print OUTFILE "\"\" ";
		}
	    }
	else
	  {
	    print OUTFILE "\"\" ";
	  }
      }
    print OUTFILE "\n";
  }

close OUTFILE;
close INDOUTFILE;
