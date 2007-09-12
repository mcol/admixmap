

workdir <- getwd()##"/ichec/work/ndlif006b/genepi/hapmap"##usually run from this dir
Panels <- c("CEU", "YRI", "JPTCHB")
CHR_NUM <- 22
LOCUS_LIMIT <- 5000##loci in small data sets
num.loci <- c(34026,38223, 32713)##loci in full data sets
num.gametes <- c(120, 120, 180)
num.masked.indivs <- c(10, 10, 15)
percent.masked.loci <- 30

taskfarm.list <- "fastphase_impute_tasks.txt"
STATES <- 6## this is a pre-emptive since we have to use some of these data to decide on the number of states

## fastPHASE options: K=number of block states, T=number of starts, C=number of EM iterations
##                    p= print parameter estimates, s=random seed, -M2=fixed mixture props
fastphase.options <- paste("-T20 -C50 -K", STATES, " -s1000 -p -M2", sep="")
Ne.values <- c(11418, 17469, 14269)## for IMPUTE

source("MaskGenotypes.R")
##set random seed, for consistency
set.seed(100)

## loop over populations
for(pop in 1:3){

  ## ###################
  ## data with 5000 loci
  ## ###################
  
  ## format raw data for HAPMIXMAP
  data.dir <- paste(paste(workdir, "data/chr22_5kloci/hapmixmap", sep="/"), Panels[pop], sep="/")
  system(paste("[ ! -d ", data.dir, " ] && mkdir -p ", data.dir, ";", sep=""))

  system(paste("FPHD -p=", workdir, "/data/chr", CHR_NUM, "/rawdata/", Panels[pop], "/chr", CHR_NUM,
               " -c=", CHR_NUM,
               " -l=", data.dir, "/phased_", LOCUS_LIMIT, "_loci.txt",
               " -g=", data.dir, "/phased_", LOCUS_LIMIT, "_genotypes.txt",
               " -n=", LOCUS_LIMIT, sep=""))

  
  ## Mask 5k genotypes and write fastphase data
  fastphase.dir <- paste(paste(workdir, "data/chr22_5kloci/fastphase", sep="/"), Panels[pop], sep="/")
  fastphase.diploid.file <- paste(fastphase.dir, "fastphase.inp", sep="/")
  fastphase.haploid.file <- paste(fastphase.dir, "fphaplotypes.inp", sep="/")

  MaskGenotypes(percent.masked.loci, num.masked.indivs[pop], paste(data.dir, "phased_5000", sep="/"),
                data.dir, 5000, num.gametes[pop], fastphase.diploid.file, fastphase.haploid.file)

  ##convert HAPMIXMAP data to IMPUTE format
  system(paste("hapmix2impute -o=",workdir, "/data/chr22_5kloci/impute/", Panels[pop],  
               " -g=", data.dir, "/train_genotypes.txt",
               " -t=", data.dir, "/test_genotypes.txt",
               " -l=", workdir, "/data/chr22/rawdata/", Panels[pop], "/chr22_legend.txt",
               " -m=", workdir, "/data/chr22/rawdata/genetic_map_chr22.txt"
               , sep=""))

  ## ################
  ## full chromosome
  ## ################

  data.dir <- paste(paste(workdir, "data/chr22/hapmixmap", sep="/"), Panels[pop], sep="/")
  
  ## format raw data for HAPMIXMAP
  system(paste("FPHD -p=", workdir, "/data/chr", CHR_NUM, "/rawdata/", Panels[pop], "/chr", CHR_NUM,
               " -c=", CHR_NUM,
               " -l=", data.dir, "/phased_loci.txt",
               " -g=", data.dir, "/phased_genotypes.txt",
               sep=""))

  ## dir to put fastphase data
  fastphase.dir <- paste(paste(workdir, "data/chr22/fastphase", sep="/"), Panels[pop], sep="/")
  fastphase.diploid.file <- paste(fastphase.dir, "fastphase.inp", sep="/")
  fastphase.haploid.file <- paste(fastphase.dir, "fphaplotypes.inp", sep="/")

  ## Mask full chr22 data and write both HAPMIXMAP and fastphase format
  MaskGenotypes(percent.masked.loci, num.masked.indivs[pop], paste(data.dir, "phased", sep="/"),
                data.dir, num.loci[pop], num.gametes[pop], fastphase.diploid.file, fastphase.haploid.file)

  ##convert HAPMIXMAP data to IMPUTE format
  system(paste("hapmix2impute -o=",workdir, "/data/chr22/impute/", Panels[pop],  
                 " -g=", data.dir, "/train_genotypes.txt",
                 " -t=", data.dir, "/test_genotypes.txt",
                 " -l=", workdir, "/data/chr22/rawdata/", Panels[pop], "/chr22_legend.txt",
                 " -m=", workdir, "/data/chr22/rawdata/genetic_map_chr22.txt"
                 , sep=""))

  ## write taskfarm list

  cat("fastPHASE ", fastphase.options,
      " -b", fastphase.haploid.file,
      " ", fastphase.diploid.file,
      "\n", file=taskfarm.list, sep="", append=T
      )
  
  impute.data.dir <- paste(workdir, "/data/chr22/impute/", Panels[pop], sep="" )
  cat("impute -h ", impute.data.dir, "/haplo.txt",
      " -g ", impute.data.dir, "/geno.txt",
      " -l ", impute.data.dir, "/legend.txt",
      " -m ", impute.data.dir, "/map.txt",
      " -o ", workdir, "/results/impute/CEU_out.txt",
      " -i ", workdir, "/results/impute/CEU_info.txt",
      " -Ne ", Ne.values[pop],
      " -os 2", ## this means output only loci in both genotype file and haplotype file, ie the masked loci
      "\n", file=taskfarm.list, sep="", append=T) 


}
