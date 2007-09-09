

workdir <- "."##"/ichec/work/ndlif006b/genepi/hapmap"##usually run from this dir
Panels <- c("CEU", "YRI", "JPTCHB")
num.loci <- c(34026,38224, 32713)
num.gametes <- c(120, 120, 180)
num.masked.indivs <- c(10, 10, 15)
percent.masked.loci <- 30

taskfarm.list <- "fastphase_impute_tasks.txt"
STATES <- 6## this is a pre-emptive since we have to use some of these data to decide on the number of states
fastphase.options <- paste("-T20 -C50 -K", STATES, " -s1000 -p -M2", sep="")
Ne.values <- c(11418, 17469, 14269)

source("MaskGenotypes.R")

## loop over populations
for(pop in 1:3){

  ## Mask 5k genotypes, no fastphase data
  data.dir <- paste(paste(workdir, "data/chr22_5kloci/hapmixmap", sep="/"), Panels[pop], sep="/")
 
  MaskGenotypes(percent.masked.loci, num.masked.indivs[pop], paste(data.dir, "phased_5000", sep="/"),
                data.dir, 5000, num.gametes[pop], "", "")

  ## Mask full chr22 data and write both HAPMIXMAP and fastphase format
  data.dir <- paste(paste(workdir, "data/chr22/hapmixmap", sep="/"), Panels[pop], sep="/")
  
  ## dir to put fastphase data
  fastphase.dir <- paste(paste(workdir, "data/chr22/fastphase", sep="/"), Panels[pop], sep="/")
  fastphase.diploid.file <- paste(fastphase.dir, "fastphase.inp", sep="/")
  fastphase.haploid.file <- paste(fastphase.dir, "fphaplotypes.inp", sep="/")

  MaskGenotypes(percent.masked.loci, num.masked.indivs[pop], paste(data.dir, "phased", sep="/"),
                data.dir, num.loci[pop], num.gametes[pop], "", "")

  ##convert HAPMIXMAP data to IMPUTE format
  command <- paste("hapmix2impute -o=",workdir, "/data/chr22/impute/", Panels[pop],  
                 " -g=", data.dir, "/train_genotypes.txt",
                 " -t=", data.dir, "/test_genotypes.txt",
                 " -l=", workdir, "/data/chr22/rawdata/", Panels[pop], "/chr22_legend.txt",
                 " -m=", workdir, "/data/chr22/genetic_map_chr22.txt"
                 , sep="")
  system(command)

  ## write taskfarm list

  cat("fastPHASE ", fastphase.options,
      " -b", fastphase.haploid.file,
      " ", fastphase.diploid.file,
      "\n", file=taskfarm.list, sep=""
      )
  
  impute.data.dir <- paste(workdir, "/data/chr22/impute/", Panels[pop], sep="" )
  cat("impute -h ", impute.data.dir, "/haplo.txt",
      " -g ", impute.data.dir, "/geno.txt",
      " -l ", impute.data.dir, "/legend.txt",
      " -m ", impute.data.dir, "/map.txt",
      " -o ", workdir, "/results/impute/CEU_out.txt",
      " -i ", workdir, "/results/impute/CEU_info.txt",
      " -Ne ", Ne.values[pop],
      "\n", file=taskfarm.list, sep="") 


}
