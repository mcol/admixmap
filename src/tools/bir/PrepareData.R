
workdir <- getwd() ##"/ichec/work/ndlif006b/genepi/hapmap"##usually run from this dir
Panels <- c("CEU", "YRI", "JPTCHB")
CHR.NUM <- 22
LOCUS.LIMIT <- 5000 ##loci in small data sets
num.loci <- c(34026,38223, 32713)# #loci in full data sets
num.gametes <- c(120, 120, 180)
num.masked.indivs <- c(10, 10, 15)
percent.masked.loci <- 20 # 30 gives too much uncertainty

taskfarm.list <- "fastphase_impute_tasks.txt"
STATES <- 6## this is a pre-emptive since we have to use some of these data to decide on the number of states

source("MaskGenotypes.R")
##set random seed, for consistency
set.seed(100)

## loop over populations
for(pop in 1:3){

  ## ################
  ## full chromosome
  ## ################

  data.dir <- paste(paste(workdir, "data/chr22/hapmixmap", sep="/"), Panels[pop], sep="/")
  
  ## format raw data for HAPMIXMAP
  system(paste("FPHD -p=", workdir, "/data/chr", CHR.NUM, "/rawdata/", Panels[pop], "/chr", CHR.NUM,
               " -c=", CHR.NUM,
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
  impute.dir <- paste(paste(workdir, "data/chr22/impute", sep="/"), Panels[pop], sep="/")
  system(paste("[ ! -d ", impute.dir, " ] && mkdir -p ", impute.dir, ";", sep=""))
  system(paste("hapmix2impute -o=",workdir, "/data/chr22/impute/", Panels[pop],  
                 " -g=", data.dir, "/train_genotypes.txt",
                 " -t=", data.dir, "/test_genotypes.txt",
                 " -l=", workdir, "/data/chr22/rawdata/", Panels[pop], "/chr22_legend.txt",
                 " -m=", workdir, "/data/chr22/rawdata/genetic_map_chr22.txt"
                 , sep=""))

  ## ###################
  ## data with 5000 loci
  ## ###################
  
  ## format raw data for HAPMIXMAP
  data.dir <- paste(paste(workdir, "data/chr22_5kloci/hapmixmap", sep="/"), Panels[pop], sep="/")
  system(paste("[ ! -d ", data.dir, " ] && mkdir -p ", data.dir, ";", sep=""))

  system(paste("FPHD -p=", workdir, "/data/chr", CHR.NUM, "/rawdata/", Panels[pop], "/chr", CHR.NUM,
               " -c=", CHR.NUM,
               " -l=", data.dir, "/phased_", LOCUS.LIMIT, "_loci.txt",
               " -g=", data.dir, "/phased_", LOCUS.LIMIT, "_genotypes.txt",
               " -n=", LOCUS.LIMIT, sep=""))

  ## BUG in FPHD: phased_genotypes file is not written correctly with only 5000 loci
  ## TEMP FIX: read full chromosome genotypes file already written, and extract first 5001 cols
  g5000 <- read.table(paste(workdir, "/data/chr22/hapmixmap/", Panels[pop],
                            "/phased_genotypes.txt", sep=""),
                      header=T,
                      nrows=num.gametes[pop],
                      colClasses=c("character", rep("integer", 5000)))
  write.table(g5000[, 1:(1 + LOCUS.LIMIT)],
              file=paste(data.dir, "/phased_", LOCUS.LIMIT, "_genotypes.txt", sep=""),
              row.names=F, col.names=T, quote=F)
  
  ## Mask 5k ge`notypes and write fastphase data
  fastphase.dir <- paste(paste(workdir, "data/chr22_5kloci/fastphase", sep="/"), Panels[pop], sep="/")
  system(paste("[ ! -d ", fastphase.dir, " ] && mkdir -p ", fastphase.dir, ";", sep=""))

  fastphase.diploid.file <- paste(fastphase.dir, "fastphase.inp", sep="/")
  fastphase.haploid.file <- paste(fastphase.dir, "fphaplotypes.inp", sep="/")

  MaskGenotypes(percent.masked.loci, num.masked.indivs[pop], paste(data.dir, "phased_5000", sep="/"),
                data.dir, 5000, num.gametes[pop], fastphase.diploid.file, fastphase.haploid.file)

  ##convert HAPMIXMAP data to IMPUTE format
  impute.dir <- paste(paste(workdir, "data/chr22_5kloci/impute", sep="/"), Panels[pop], sep="/")
  system(paste("[ ! -d ", impute.dir, " ] && mkdir -p ", impute.dir, ";", sep=""))
  system(paste("hapmix2impute -o=", impute.dir,  
               " -g=", data.dir, "/train_genotypes.txt",
               " -t=", data.dir, "/test_genotypes.txt",
               " -l=", workdir, "/data/chr22/rawdata/", Panels[pop], "/chr22_legend.txt",
               " -m=", workdir, "/data/chr22/rawdata/genetic_map_chr22.txt"
               , sep=""))

}
