mkdir -p -v "Eur/chr22data"
cp -v ~/shared-genepi/genepi/hapmap-source/Eur/chr22data/phased_{genotypes,loci}.txt Eur/chr22data
perl chr22.pl --mask-data --pop=Eur --percent-indivs 17 --percent-loci 30 --limit-loci 0 --genotypes-file mi_genotypes.txt --locus-file phased_loci.txt --maskfile mi_cc_index.txt --case-control-file mi_cc.txt
if [ "$?" != "0" ]; then echo "Perl script returned an error."; exit 1; fi
mkdir -p -v "Afr/chr22data"
cp -v ~/shared-genepi/genepi/hapmap-source/Afr/chr22data/phased_{genotypes,loci}.txt Afr/chr22data
perl chr22.pl --mask-data --pop=Afr --percent-indivs 17 --percent-loci 30 --limit-loci 0 --genotypes-file mi_genotypes.txt --locus-file phased_loci.txt --maskfile mi_cc_index.txt --case-control-file mi_cc.txt
if [ "$?" != "0" ]; then echo "Perl script returned an error."; exit 1; fi
mkdir -p -v "Asian/chr22data"
cp -v ~/shared-genepi/genepi/hapmap-source/Asian/chr22data/phased_{genotypes,loci}.txt Asian/chr22data
perl chr22.pl --mask-data --pop=Asian --percent-indivs 17 --percent-loci 30 --limit-loci 0 --genotypes-file mi_genotypes.txt --locus-file phased_loci.txt --maskfile mi_cc_index.txt --case-control-file mi_cc.txt
if [ "$?" != "0" ]; then echo "Perl script returned an error."; exit 1; fi
