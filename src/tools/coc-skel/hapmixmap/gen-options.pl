# script to generate options files for multiple runs of hapmixmap with different prior or different numbers of states

my $dataprefix="/ichec/HapMap/data";
my @Panels=("CEU", "YRI", "JPTCHB");
my $resultsprefix="/ichec/HapMap/results";
my $resultsdir="Chr22Results";
my $chromosome = 22;

my $arg_hash = {
    deleteoldresults => 0,

#data
    genotypesfile                   => "$datadir/phased_genotypes.txt",
    locusfile                       => "$datadir/phased_loci.txt",

#model
    states     => 6,
    checkdata       => 0,
    # mixturepropsprior => "50, 1",
    fixedmixtureprops => 1,
    fixedmixturepropsprecision =>1,

#main options
    resultsdir      => 'results',
    displaylevel    => 3,
    samples         => $samples,
    burnin          => $burnin,
    every           => $every,
    numannealedruns => 0,
    thermo          => 0,
    hapmixmodel     => 1,

#prior spec
    freqprecisionhiermodel => 0,
    arrivalratesamplerparams => "0.1, 0.00001, 10, 0.9, 20",

		#output files
		logfile =>'logfile.txt',
    paramfile         =>'paramfile.txt',#mean and var of sampled arrival rates
    freqprecisionfile =>'freqprecision.txt', #mean and var of sampled allele freq precision
    arrivalrateposteriormeanfile => "ArrivalRatePosteriorMeans.txt",
		predictgenotypes => 1
};

sub writeOptionsFile{
    my $hash = $_[0];
    my $filename = $_[1];
    open(OPTIONFILE, ">$filename") or die ("Could not open options file");
    foreach my $key (keys %$hash){
      print OPTIONFILE $key . '=' . $hash->{$key} . "\n";
    }
    close OPTIONFILE;
};

for my $pop(@Panels){
## NB '$_' means the current panel code (CEU/YRI/JPTCHB)
print "$pop\n";
#data files
$arg_hash->{locusfile}="$dataprefix/$pop/hapmixmap/chr_"."$chromosome"."_phased_loci.txt";
$arg_hash->{genotypesfile}="$dataprefix/$pop/hapmixmap/chr_"."$chromosome"."_genotypes_train.txt";
$arg_hash->{testgenotypesfile}="$dataprefix/$pop/hapmixmap/chr_"."$chromosome"."_genotypes_masked.txt";

$arg_hash->{arrivalrateprior}="1.2, 0.05, 1, 0.5";
$arg_hash->{residualadprior}="2, 8";
$arg_hash->{resultsdir}="$resultsprefix/$pop/arp-1.2-0.05-1-0.5-afpp-0.25-1/$resultsdir";
writeOptionsFile($arg_hash, "$pop-arp-1.2-0.05-1-0.5-afpp-0.25-1.txt");

}
