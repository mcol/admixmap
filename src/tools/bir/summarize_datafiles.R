## script to read all data files and print their dimensions

system("find data/chr22*/hapmixmap -name '*.txt' > filenames.txt")

filenames <- read.table(file="filenames.txt", header=F, as.is=T)[, 1]


for(file in 1:length(filenames)) {
  y <- try(read.table(filenames[file], header=T, fill=T))
  
  cat(filenames[file], dim(y), "\n")
  cat(filenames[file], dim(y), "\n", file="datatables.txt", append=T)
  rm(y)
}

