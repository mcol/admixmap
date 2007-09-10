
Step-by-step guide to the bir analyses
-> denotes files produced at each step


 -- Prepare data --

* download HapMap data (chr 22, all 3 panels)
  Use getdata.pl
  We need a wrapper script to get all 3 pops (and all chromosomes?)
  -> 3 files per panel

* download recombination rates (genetic map) from HapMap
  Must be done manually (only one file)

* format for HAPMIXMAP
  (full chromosome) Use FPHD unrestricted
  -> a genotypesfile and a locusfile for each panel
     phased_genotypes.txt, phased_loci.txt
  (5k loci) Use FPHD with locus limit of 5000
  -> phased_5000_loci.txt, phased_5000_genotypes.txt

* Mask genotypes
  Use MaskGenotypes.R
  -> hapmixmap data: train_genotypes.txt, train_loci.txt, test_genotypes.txt; 
     observed allele counts at masked loci among masked individuals: obs_counts.txt;
     observed genotypes of masked individuals at masked loci (coded as 1, 2, 3): obs_masked_genotypes.txt; 
  -> 2 fastPHASE input files
  -> list of masked loci: masked_loci.txt 

* prepare fastphase data
  use functions defined in hapmix2fastphase.R

* prepare IMPUTE data (full chromosome)
  use hapmix2impute
  -> 4 files haplo.txt, geno.txt, legend.txt, map.txt

Above 4 steps may be combined with PrepareData.R. This will also write a list of tasks (ie commands to run fastPHASE and IMPUTE) for a taskfarm.

Data are organised as follows:
HapMap/data/chr22
           /chr22_5kloci
                /rawdata
                /hapmixmap
                /impute
                /fastphase /CEU
                           /YRI
                           /JPTCHB

     /results/hapmixmap/arp-0.5-1-2-radp-0.25-1-s1
             /impute
             /fastphase/ .... /CEU
                              /YRI
                              /JPTCHB
      


 -- comparisons of different numbers of states (5k loci) --
 -- comparisons of different priors (5k loci) --

   - Use gen-options.pl to write options files for HAPMIXMAP 
   - write PBS script and submit
   - run bir.sh, which will read a list of resultsdirs and run BIR.R in each, calculating BIR then collect them into a table

 -- estimates of parameters with weak priors (full chromosome) --

   - write 3 options files with only the panel varying
   - write PBS script and submit
   - run AdmixmapOutput.R on each resultsdir
   - extract parameter estimates from PosteriorQuantiles.txt
   - use hotspotmap.R to plot 


 -- comparisons with fastPHASE and IMPUTE (full chromosome) --

  - write 3 options files for HAPMIXMAP, one for each panel, with the number of states and the priors that give best prediction.
  - write taskfarm list for fastPHASE and IMPUTE. This could be done during data preparation but is also easy to do by hand.
  - write PBS script and submit
  - convert IMPUTE and fastPHASE output to same format as HAPMIXMAP's
  - run bir.sh again with a new list of resultsdirs


 -- Colon Cancer data (2 sets - Hap550 and Hap300) --

 - use FPHD to prepare data, supplying preformatted case/control genotypes
   -> locusfile, genotypesfile, testgenotypesfile
 - write 2 options files, using the optimal number of states and priors as determined above and the score test on
 - write taskfarm list (with 2 jobs !)
 - run AdmixmapOutput on both resultsdirs
 - use hapmixplots.R to produce smoothed plot of info extracted in each array against map position
 - run hetcounts program to produce tables of heterozygosity etc
 - use hapmixplots.R to produce plot of mean info extracted and mean missing info against heterozygosity in each array
