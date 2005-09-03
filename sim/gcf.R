


Fdr.Func <- function(p,alpha)
{

   m <- length(p)
   p <- sort(p)
   p <- c(0,p)
   p[p==1] <- 1-10e-10
   G.hat <- (0:m)/m
   G.hat[G.hat<p] <- p[G.hat < p]
   k <- floor(m/2) + 1
   a.hat <- max( (G.hat[1:k] - p[1:k])/(1-p[1:k]) )
   thresh <- max(p[ (p <= (0:m)*alpha/(m*(1-a.hat))) ])
   list(thresh=thresh,a.hat=a.hat)
}


Error.Ch <- function(string, out)
{
   cat("\n\n--------------------------------------------------------------------------", file=out, append=T)
   cat("\n\n ERROR: ", string, file=out, append=T)  
   return(1)
}


Main.Func <- function()
{

   errorflag     <- 0
   driver        <- NULL
   driverfile    <- "gcdriver.txt"
#   cat("\nPlease enter name of driver file (no quotes) below then hit Enter twice:\n")
#   driverfile    <- scan(what='character')
#   
#   warntemp      <- options()$warn
#   options(warn=-1)
#   tryreaddriver <- try(read.table(file=driverfile),silent=T)
#
#   if (attr(tryreaddriver,"class")=="try-error")
#   {
#      cat("\n\n\nERROR: Driver file does not exist or can not be read. Aborting.\n\n\n")
#      return()
#   }
      
   driver     <- read.table(file=driverfile)
   drivertemp <- strsplit(as.vector(driver$V1), "=")

   datafile   <- NA
   covfile    <- NA
   nullloci   <- NA
   model      <- NA
   modeltype  <- NA
   alpha      <- NA
   cov        <- NA
   outmodel   <- NA
   outnull    <- NA

   for (i in 1: length(drivertemp))
   {

      if (drivertemp[[i]][1] == "DATAFILE")  datafile  <- drivertemp[[i]][2]
      if (drivertemp[[i]][1] == "COVFILE")   covfile   <- drivertemp[[i]][2]
      if (drivertemp[[i]][1] == "NULLLOCI")  nullloci  <- drivertemp[[i]][2]
      if (drivertemp[[i]][1] == "MODEL")     model     <- drivertemp[[i]][2]
      if (drivertemp[[i]][1] == "MODELTYPE") modeltype <- drivertemp[[i]][2]
      if (drivertemp[[i]][1] == "ALPHA")     alpha     <- drivertemp[[i]][2]
      if (drivertemp[[i]][1] == "COVLABEL")  cov       <- drivertemp[[i]][2]
      if (drivertemp[[i]][1] == "OUTMODEL")  outmodel  <- drivertemp[[i]][2]
      if (drivertemp[[i]][1] == "OUTNULL")   outnull   <- drivertemp[[i]][2]
   }

   if (is.na(outmodel))
   {
      cat("\n\n\nERROR: Output file not provided in driver file. Aborting.\n\n\n")
      return()
   }

   cat("\nStatus: Error checking")
   cat("\n--------------------------------------------------------------------------", file=outmodel, append=F)
   cat("\nAssociation Studies for Quantitative Traits in Structured Populations:", file=outmodel, append=T)
   cat("\nModels and Summary Statistics", file=outmodel, append=T)
   cat("\n\nDate:",date(), file=outmodel, append=T)
   cat("\n--------------------------------------------------------------------------", file=outmodel, append=T)

   dataflag    <- 1
   tryreaddata <- try(read.table(file=datafile,sep=""),silent=T)

   if ((attr(tryreaddata,"class")=="try-error") | (is.na(datafile)))
   {
      dataflag <- 0
      errorflag <- Error.Ch("Data file (DATAFILE=) not specified or does not exist.", outmodel)
   }

   readdata    <- NULL
   locimattemp <- NULL
   outvect     <- NULL

   if (dataflag == 1)
   {
      readdata    <- read.table(file=datafile,sep="")
      locimattemp <- cbind(readdata[,3:(dim(readdata)[2])])
      outvect     <- readdata[,2]
      if ((dim(locimattemp)[2] %% 2) == 1)
      {
         dataflag  <- 0
         errorflag <- Error.Ch("Must have even number of loci columns in data file.", outmodel)
      }
   }
      
   if (dataflag == 1)
   {
      locimat  <- matrix(rep(0, dim(locimattemp)[1]*(dim(locimattemp)[2]/2)),
                         ncol=dim(locimattemp)[2]/2)


      for (i in 1:(dim(locimattemp)[2]/2))
      {
         locimat[,i]  <- locimattemp[,(2*i-1)] + locimattemp[,(2*i)] 
      }

      locimatch <- locimat
      locimatch[(locimatch==4) | (locimatch==1) | (locimatch==2) | (locimatch==3) | (locimatch==0)] <- 0

      if (sum(locimatch) > 0)
      {
         dataflag  <- 0
         errorflag <- Error.Ch("Unexpected values observed in data file (DATAFILE=).", outmodel)
      }
         
      locimat[locimat == 0] <- NA
      locimat[locimat == 2] <- 0
      locimat[locimat == 3] <- 1
      locimat[locimat == 4] <- 2
                                        
      locimatc <- locimat - (matrix(rep(colMeans(locimat,na.rm=T),each=dim(locimat)[1]),
                                    ncol=dim(locimat)[2]))
   }

   covflag <- ifelse(is.na(covfile), 0, 1)
   tryreadcov <- try(read.table(file=covfile,sep=""),silent=T)

   if ((attr(tryreadcov,"class")=="try-error") & (covflag==1))
   {
      covflag <- 0
      errorflag <- Error.Ch("Covariate file (COVFILE=) does not exist.", outmodel)
   }

   readcov   <- NULL
   covmat    <- NULL
   covlabels <- NA

   if (covflag == 1)
   {
      readcov       <- read.table(file=covfile,sep="")
      covmatu       <- cbind(readcov[,2:(dim(readcov)[2])])
      covmat        <- as.data.frame(covmatu - (matrix(rep(colMeans(covmatu,na.rm=T),each=dim(covmatu)[1]),
                                                ncol=dim(covmatu)[2])))
      covlabels     <- c(unlist(strsplit(cov, ",")))

      if (dim(covmat)[2] != length(covlabels))
         errorflag <- Error.Ch("Incorrect number of covariate labels (COVLABEL=).", outmodel)
      else if ((!is.null(covmat)) & (!is.na(cov))) names(covmat) <- covlabels
   }
   
   if ((dataflag == 1) & (covflag == 1))
   {
      if (dim(locimatc)[1] != dim(covmat)[1])
         errorflag <- Error.Ch("Data file and Covariate file have unequal number of observations.", outmodel)
   }

   if (covflag == 0) ncov <- 0
   if (covflag == 1) ncov <- dim(covmat)[2]

   modelflag  <- 1
   tempmodel1 <- NA

   if (length(strsplit(model, "~")[[1]]) > 1)
   {
      tempvec = c(unlist(strsplit(model, ",")))
      newvec = NA
      for (i in 1:length(tempvec))
      {
         if (length(strsplit(tempvec[i], "~")[[1]]) > 1)
         {
            tempset = as.numeric(unlist(strsplit(tempvec[i], "~")))
            for (j in min(tempset):(max(tempset)-1))
            {
               for (k in (j+1):(max(tempset)))
               {
                  newvec = c(newvec,paste(j,k,sep=":"))
               }
            }
         }
         else newvec = c(newvec,tempvec[i])
      }
      model = paste(newvec[2:length(newvec)],collapse=",")
   }

   tempmodel1 <- strsplit(c(unlist(strsplit(model, ","))), ":")

   if (is.na(model)) tempmodel1 <- NA
                    
   if (any(is.na(unlist(tempmodel1))))
   {
      errorflag <- Error.Ch("Need to specify model loci.", outmodel)
      modelflag <- 0
   }

   margint <- NULL
   if (length(unlist(tempmodel1)) == length(tempmodel1) )   margint <- 1
   if (length(unlist(tempmodel1)) == 2*length(tempmodel1) ) margint <- 4


   if (is.null(margint))
   {
      errorflag <- Error.Ch("Need to specify only marginal or interaction models.", outmodel)
      modelflag <- 0
   }

   if (modelflag == 1)
   {
      if (margint == 1)
         tempmodel2 <- strsplit(c(unlist(tempmodel1)), "-")

      if (margint == 4)
         tempmodel2 <- tempmodel1
   }

   if ( (is.na(as.numeric(paste(unlist(tempmodel2),collapse="")))) |
        ( any(is.na(as.numeric(unlist(tempmodel2)))) ) )
   {
      errorflag <- Error.Ch("Unexpected characters found in specification of model loci.", outmodel)
      modelflag <- 0
   }
   
   if ((dataflag == 1) & (modelflag == 1))
   {
      if ((max(as.numeric(unlist(tempmodel2))) > dim(locimatc)[2]) | 
          (min(as.numeric(unlist(tempmodel2))) < 1))
         errorflag <- Error.Ch("Model loci as specified in driver file do not exist.", outmodel)
   }

   nullflag <- 1
   tempnull <- NA
   if (is.na(nullloci))
      nullloci <- paste(1,dim(locimatc)[2],sep="-")
   tempnull <- strsplit(c(unlist(strsplit(nullloci, ","))), "-")

   if (any(is.na(unlist(tempnull))))   nullflag  <- 0

   if ((is.na(as.numeric(paste(unlist(tempnull),collapse="")))) & (nullflag == 1)) 
   {
      errorflag <- Error.Ch("Unexpected characters found in specification of null loci.", outmodel)
      nullflag <- 0
   }

   if ((dataflag == 1) & (nullflag == 1))
   {
      if ((max(as.numeric(unlist(tempnull))) > dim(locimatc)[2]) | 
          (min(as.numeric(unlist(tempnull))) < 1))
         errorflag <- Error.Ch("Null loci as specified in driver file do not exist.", outmodel)
   }
   
   fam   <- as.numeric(modeltype)

   if ((fam == 0) & (length(attr(as.factor(outvect),"levels")) != 2) )
      errorflag <- Error.Ch("Need 2 levels for dichotomous response variable", outmodel)   

   alpha <- as.numeric(alpha)


   if ((is.na(fam)) | ((fam != 1) & (fam != 0)))
      errorflag <- Error.Ch("Need to specify type of response variable (MODELTYPE=).", outmodel)


   if ((is.na(alpha)) | (alpha <= 0) | (alpha >= 1))
      errorflag <- Error.Ch("Alpha invalid", outmodel)


   if (errorflag == 1)
   {
      cat("\n\n\nError(s) occured - see output file" , outmodel, "\n\n\n")
      return()
   }


   cat("\nStatus: Computing estimate(s) of lambda - 0%")
   seqnull      <- NULL
   nullmat      <- NULL
   tnullmarg    <- NULL
   tnullint     <- NULL
   lamhatmarg   <- 1
   lamhatint    <- 1
   lamhatomn    <- 1

   if (nullflag == 1)
   {
      for (i in 1:length(tempnull))
      {
         tempmat   <- matrix(as.numeric(tempnull[[i]]),ncol=2)
         seqnull   <- c(seqnull,tempmat[1,1]:tempmat[1,2])
      }
      nullmat      <- locimatc[,seqnull]
      tnullmarg    <- as.data.frame(matrix(nrow=dim(nullmat)[2],ncol=3))
      names(tnullmarg) <- c("nullmodel", "stat", "pval")

      if (margint == 1) nnullops <- dim(nullmat)[2]
      if (margint == 4) nnullops <- dim(nullmat)[2] + choose(dim(nullmat)[2], 2)
      nnullops <- floor(nnullops) + (4-(nnullops %% 4))
      cnullops <- 0

      for (i in 1:dim(tnullmarg)[1])
      {
         if (covflag == 1)
         {
            fmlatemp <- paste("covmat$",paste(names(covmat),collapse=" + covmat$"), sep="")
            fmla <- as.formula(paste("outvect ~ nullmat[,i]", fmlatemp, sep=" + "))
         }
         else if (covflag == 0)
         {
            fmla <- as.formula("outvect ~ nullmat[,i]")     
         }
         if (fam == 0)
            lmtemp <- glm(fmla, na.action=na.omit,family=binomial(link="logit"))
         else if (fam == 1)
            lmtemp <- glm(fmla, na.action=na.omit,family=gaussian(link ="identity"))
         tnullmarg[i,1] <- seqnull[i]
         ## summary(lmtemp)$coef) is a matrix with only one row if all genotypes same
         if(dim(summary(lmtemp)$coef)[1] == 1) {
           tnullmarg[i,2:3] <- NA
         } else {
           tnullmarg[i,2] <- summary(lmtemp)$coef[2,1] / summary(lmtemp)$coef[2,2]
           tnullmarg[i,3] <- summary(lmtemp)$coef[2,4]
         }
           cnullops <- cnullops + 1
         if (cnullops/nnullops == 0.25)
            cat("\nStatus: Computing estimate(s) of lambda - 25%")
         if (cnullops/nnullops == 0.5)
            cat("\nStatus: Computing estimate(s) of lambda - 50%")
         if (cnullops/nnullops == 0.75)
            cat("\nStatus: Computing estimate(s) of lambda - 75%")
      }

      if (margint == 4)
      {
         tnullint     <- as.data.frame(matrix(nrow=choose(dim(nullmat)[2],2),ncol=6))
         names(tnullint) <- c("nullmodel", "flag", "inttermstat", "omnibusstat", "inttermpval", "omnibuspval")

         for (i in 1:(length(seqnull)-1))
         {
            for (j in (i+1):length(seqnull))
            {
               if (covflag == 1)
               {
                  fmlatemp <- paste("covmat$",paste(names(covmat),collapse=" + covmat$"), sep="")
                  fmla <- as.formula(paste("outvect ~ nullmat[,i] + nullmat[,j] + nullmat[,i]:nullmat[,j]",
                                        fmlatemp, sep=" + "))
                  fmlanull <- as.formula(paste("outvect ~ ", fmlatemp, sep=" + "))
               }
               else if (covflag == 0)
               {
                  fmla <- as.formula("outvect ~ nullmat[,i] + nullmat[,j] + nullmat[,i]:nullmat[,j]")     
                  fmlanull <- as.formula("outvect ~ 1")
               }
               if (fam == 0)
               {
                  lmtemp     <- glm(fmla,     na.action=na.omit,family=binomial(link="logit"))
                  lmtempnull <- glm(fmlanull, na.action=na.omit,family=binomial(link="logit"),
                     subset=((is.na(nullmat[,i])==F) & (is.na(nullmat[,j])==F)))
               }
               else if (fam == 1)
               {
                  lmtemp     <- glm(fmla,     na.action=na.omit,family=gaussian(link ="identity"))
                  lmtempnull <- glm(fmlanull, na.action=na.omit,family=gaussian(link ="identity"),
                     subset=((is.na(nullmat[,i])==F) & (is.na(nullmat[,j])==F)))
               }
               row <- (i-1)*length(seqnull) - ((i-1)*i/2) + (j-i)
               tnullint[row,1] <- paste(c(seqnull[i],seqnull[j]),collapse=":")

               clmtemp <- summary(lmtemp)$coef

               if (dim(clmtemp)[1] == (4+ncov))
               {
                  tnullint[row,2] <- 1
                  tnullint[row,3] <- clmtemp[dim(clmtemp)[1],1] / clmtemp[dim(clmtemp)[1],2]
                  tnullint[row,5] <- clmtemp[dim(clmtemp)[1],4]

                  almtemp     <- anova(lmtemp)
                  almtempnull <- anova(lmtempnull)

                  numss   <- (almtempnull[dim(almtempnull)[1],dim(almtempnull)[2]]   -
                              almtemp[dim(almtemp)[1],dim(almtemp)[2]])
                  numdf   <- (almtempnull[dim(almtempnull)[1],dim(almtempnull)[2]-1] -
                              almtemp[dim(almtemp)[1],dim(almtemp)[2]-1])
                  num     <- numss/numdf
                  denss   <- almtemp[dim(almtemp)[1],dim(almtemp)[2]]
                  dendf   <- almtemp[dim(almtemp)[1],dim(almtemp)[2]-1]
                  den     <- denss/dendf
                  tnullint[row,4] <- num / den
                  tnullint[row,6] <- 1-pf(tnullint[row,4],numdf,dendf)
               }
               else tnullint[row,2] <- -9
               cnullops <- cnullops + 1
               if (cnullops/nnullops == 0.25)
                  cat("\nStatus: Computing estimate(s) of lambda - 25%")
               if (cnullops/nnullops == 0.5)
                  cat("\nStatus: Computing estimate(s) of lambda - 50%")
               if (cnullops/nnullops == 0.75)
                  cat("\nStatus: Computing estimate(s) of lambda - 75%")
            }
         }
      }   

      lamhatmarg <- mean(tnullmarg[,2]^2,na.rm=T)
      lamhatint <- NULL
      if (margint == 4)
      {
         lamhatint  <- mean(tnullint[,3]^2,na.rm=T)
         lamhatomn  <- mean(tnullint[,4],  na.rm=T)
      }

   }
   
   tnullmargfs <- formatC(tnullmarg[,2], digits=3,format="f")
   tnullmargfp <- formatC(tnullmarg[,3], digits=4,format="f")
   tnullintfis <- formatC(tnullint[,3],  digits=3,format="f")
   tnullintfos <- formatC(tnullint[,4],  digits=3,format="f")
   tnullintfip <- formatC(tnullint[,5],  digits=4,format="f")
   tnullintfop <- formatC(tnullint[,6],  digits=4,format="f")

   cat("\nStatus: Computing estimate(s) of lambda - 100%")
   cat("\nStatus: Running models")

   tempmodel <- NULL
   seqmodel  <- NULL
   if (margint == 1)
   {
      for (i in 1:length(tempmodel2))
      {
         tempmat    <- matrix(as.numeric(tempmodel2[[i]]),ncol=2)
         seqmodel   <- c(seqmodel,tempmat[1,1]:tempmat[1,2])
      }
      tempmodel <- as.list(seqmodel)
   }
   else tempmodel <- tempmodel2 
   
   stats  <- matrix(NA, ncol=(1+margint+ncov), nrow=length(tempmodel))
   pvals  <- matrix(NA, ncol=(1+margint+ncov), nrow=length(tempmodel))


   for (i in 1:length(tempmodel))
   {
      if (margint == 1)
      {
         fmlatemp1 <- "outvect ~ locimatc[,as.numeric(tempmodel[[i]][1])]"
         if (covflag == 1)
         {
            fmlatemp2 <- paste("covmat$",paste(names(covmat),collapse=" + covmat$"), sep="") 
            fmla <- as.formula(paste(fmlatemp1, fmlatemp2, sep=" + "))
         }
         else if (covflag == 0)
         {
            fmla <- as.formula(fmlatemp1)     
         }
         if (fam == 0)
            lmtemp <- glm(fmla, na.action=na.omit, family=binomial(link="logit"))
         else if (fam == 1)
            lmtemp <- glm(fmla, na.action=na.omit, family=gaussian(link="identity"))
         if (dim(summary(lmtemp)$coef)[1] == (2+ncov))
         {
            stats[i,1] <- 1
            stats[i,2] <- summary(lmtemp)$coef[2,1] / (sqrt(lamhatmarg) * summary(lmtemp)$coef[2,2])
            pvals[i,2]  <- 1-pchisq((stats[i,2]^2),df=1)
            if (ncov > 0)
            {
               for (j in 3:(3+ncov-1))
               {
                  stats[i,j]   <- summary(lmtemp)$coef[j,1] / (sqrt(lamhatmarg) * summary(lmtemp)$coef[j,2])  
                  pvals[i,j]   <- 1-pchisq((stats[i,j]^2),df=1)
               }  
            }
         }
         else stats[i,1] <- -9
      }

      if (margint == 4)
      {      
         fmlatemp1 <- paste(
            "outvect ~ locimatc[,as.numeric(tempmodel[[i]][1])]",
            "locimatc[,as.numeric(tempmodel[[i]][2])]",
            "locimatc[,as.numeric(tempmodel[[i]][1])]:locimatc[,as.numeric(tempmodel[[i]][2])]",
            sep=" + ")
         if (covflag == 1)
         {
            fmlatemp2 <- paste("covmat$",paste(names(covmat),collapse=" + covmat$"), sep="") 
            fmla <- as.formula(paste(fmlatemp1, fmlatemp2, sep=" + "))
            fmlanull <- as.formula(paste("outvect ~ ", fmlatemp2, sep=" + "))
         }
         else if (covflag == 0)
         {
            fmla     <- as.formula(fmlatemp1)        
            fmlanull <- as.formula("outvect ~ 1") 
         }
         if (fam == 0)
         {
            lmtemp <- glm(fmla, na.action=na.omit, family=binomial(link="logit"))
            lmtempnull <- glm(fmlanull, na.action=na.omit,family=binomial(link="logit"),
               subset=((is.na(locimatc[,as.numeric(tempmodel[[i]][1])] )==F) &
                       (is.na(locimatc[,as.numeric(tempmodel[[i]][2])] )==F)))
         }
         else if (fam == 1)
         {
            lmtemp <- glm(fmla, na.action=na.omit, family=gaussian(link ="identity")) 
            lmtempnull <- glm(fmlanull, na.action=na.omit,family=gaussian(link ="identity"),
               subset=((is.na(locimatc[,as.numeric(tempmodel[[i]][1])] )==F) &
                       (is.na(locimatc[,as.numeric(tempmodel[[i]][2])] )==F)))
         }

         clmtemp <- summary(lmtemp)$coef
         if (dim(clmtemp)[1] == (4+ncov))
         {
            stats[i,1]        <- 1
            stats[i,2]        <- clmtemp[2,1] / (sqrt(lamhatmarg) * clmtemp[2,2]) 
            stats[i,3]        <- clmtemp[3,1] / (sqrt(lamhatmarg) * clmtemp[3,2]) 
            stats[i,4]        <- clmtemp[(4+ncov),1] / (sqrt(lamhatint) * clmtemp[(4+ncov),2]) 

            pvals[i,2]         <- 1-pf((stats[i,2]^2),1,dim(tnullmarg)[1])
            pvals[i,3]         <- 1-pf((stats[i,3]^2),1,dim(tnullmarg)[1])
            pvals[i,4]         <- 1-pf((stats[i,4]^2),1,dim(tnullmarg)[1])
         
            if (ncov > 0)
            {
               for (j in 6:(6+ncov-1))
               {
                  stats[i,j]  <- clmtemp[(j-1),1] / (sqrt(lamhatmarg) * clmtemp[(j-1),2]) 
                  pvals[i,j]   <- 1-pf((stats[i,j]^2),1,dim(tnullmarg)[1])
               }  
            }
            almtemp <- anova(lmtemp)
            almtempnull <- anova(lmtempnull)
            numss   <- (almtempnull[dim(almtempnull)[1],dim(almtempnull)[2]]   -
                        almtemp[dim(almtemp)[1],dim(almtemp)[2]])
            numdf   <- (almtempnull[dim(almtempnull)[1],dim(almtempnull)[2]-1] -
                        almtemp[dim(almtemp)[1],dim(almtemp)[2]-1])
            num     <- numss/numdf
            denss   <- almtemp[dim(almtemp)[1],dim(almtemp)[2]]
            dendf   <- almtemp[dim(almtemp)[1],dim(almtemp)[2]-1]
            den     <- denss/dendf
            stats[i,5]  <- num / den         
            pvals[i,5]  <- 1-pf((stats[i,5]/lamhatomn),numdf,dim(tnullmarg)[1])
         }
         else stats[i,1] <- -9
      }
   }

   pvalsf <- formatC(pvals,digits=4,format="f")
   statsf <- formatC(stats,digits=3,format="f")

   if (margint == 1)
      pvalvec      <- as.vector(na.exclude(as.vector(pvals[,2])))
   if (margint == 4)
      pvalvec      <- as.vector(na.exclude(as.vector(pvals[,5])))


   bonadj       <- length(pvalvec)
   boncutoff    <- alpha/bonadj
   bonmat       <- matrix(" ", nrow=dim(pvals)[1],ncol=dim(pvals)[2])
   bonmat[pvals <= boncutoff] <- "Yes"
   bonmat

   fdr          <- Fdr.Func(pvalvec,alpha)
   fdrcutoff    <- fdr$thresh
   if (fdrcutoff == 0)
      fdrcutoff <- NA
   fdrmat       <- matrix(" ", nrow=dim(pvals)[1],ncol=dim(pvals)[2])
   fdrmat[pvals <= fdrcutoff] <- "Yes"
   fdrmat

   cat("\nStatus: Printing output")
   if (fam==0)
   famstring <- "Binomial"

   if (fam==1)
   famstring <- "Gaussian"

   cat("\n\n--------------------------------------------------------------------------", file=outmodel, append=T)
   cat("\nProgram: GCF", file=outmodel, append=T)
   cat("\nResponse Type (Y):", famstring, file=outmodel, append=T)
   if (covflag == 0)
      cat("\nCovariates: No", file=outmodel, append=T)
   if (covflag == 1)
      cat("\nCovariates: Yes", file=outmodel, append=T)
   cat("\nNumber of loci in data file:", dim(locimattemp)[2]/2, file=outmodel, append=T)
   cat("\nNumber of subjects:", dim(locimatc)[1], file=outmodel, append=T)
   if (nullflag == 1)
   {
      cat("\nNull Loci:", nullloci, file=outmodel, append=T)
      cat("\nEstimate(s) of lambda:", file=outmodel, append=T)
      cat("\n   Marginal:   ", round(lamhatmarg,digits=5), file=outmodel, append=T)
      if (margint == 4)
      {
         cat("\n   Interaction:", round(lamhatint,digits=5), file=outmodel, append=T)
         cat("\n   Omnibus:    ", round(lamhatomn,digits=5), file=outmodel, append=T)
      }
   }
   
   if (nullflag == 0)
   {
      cat("\nNo null loci specified: Estimate(s) for lambda set to 1.", file=outmodel, append=T)
   }
   if (is.na(outnull) == T)
   cat("\nNOTE: No output file specified for null models.", file=outmodel, append=T)
   cat("\nAlpha:", alpha, file=outmodel, append=T)
   cat("\nBonferroni p-value cutoff:", formatC(boncutoff,digits=4,format="f"), file=outmodel, append=T)
   if (is.na(fdrcutoff))
      cat("\nFDR p-value cutoff: NA", file=outmodel, append=T)
   else
      cat("\nFDR p-value cutoff:", formatC(fdrcutoff,digits=4,format="f"), file=outmodel, append=T)
   cat("\n--------------------------------------------------------------------------", file=outmodel, append=T)


   for (i in 1:dim(stats)[1])
   {
      cat("\n\n\n----------", file=outmodel, append=T)
      cat("\nModel", i, file=outmodel, append=T)   
      cat("\n----------", file=outmodel, append=T)

      if (margint == 1)
      {
         cat("\nModel: Y = X", tempmodel[[i]][1], file=outmodel, append=T, sep="")
 
         if (stats[i,1] == -9)
         {
            cat("\n\nNOTE: A problem occured when fitting this model", file=outmodel, append=T)
         }

         else
         {
            cat("\n\n                   Test                              Bonferroni     FDR", file=outmodel,append=T)
            cat("\n  Term            statistic           p-value          cutoff      cutoff", file=outmodel,append=T) 
            cat("\n--------------------------------------------------------------------------", file=outmodel,append=T) 
            cat("\n  X", tempmodel[[i]][1],
                paste(as.character(rep(" ",22-nchar(tempmodel[[i]][1]) - nchar(statsf[i,2]) )),collapse=""),
                statsf[i,2], "             ",
                pvalsf[i,2], "           ",
                bonmat[i,2],
                paste(as.character(rep(" ",13-nchar(bonmat[i,2]) )),collapse=""), 
                fdrmat[i,2],
                file=outmodel,append=T,sep="")

            if (ncov > 0)
            {
               for (j in 3:(3+ncov-1))
               {
                  cat("\n  ", names(covmat)[(j-2)],
                      paste(as.character(rep(" ",23-nchar(names(covmat)[(j-2)]) - nchar(statsf[i,j]) )),collapse=""),
                      statsf[i,j], "             ",
                      pvalsf[i,j], "           ",
                      paste(as.character(rep(" ",13-nchar(bonmat[i,j]) )),collapse=""), 
                      file=outmodel,append=T,sep="")
               }  
            }
         }
      }

      if (margint == 4)
      {
         cat("\nModel: Y = X", tempmodel[[i]][1], " + X", tempmodel[[i]][2], " + X", tempmodel[[i]][1],
             "*X", tempmodel[[i]][2],
             file=outmodel, append=T, sep="")
 
         if (stats[i,1] == -9)
         {
            cat("\n\nNOTE: A problem occured when fitting this model", file=outmodel, append=T)
         }

         else
         {
            cat("\n\n                   Test                              Bonferroni     FDR", file=outmodel,append=T)
            cat("\n  Term            statistic           p-value          cutoff      cutoff", file=outmodel,append=T) 
            cat("\n--------------------------------------------------------------------------", file=outmodel,append=T) 
            cat("\n  X", tempmodel[[i]][1],
                paste(as.character(rep(" ",22-nchar(tempmodel[[i]][1]) - nchar(statsf[i,2]) )),collapse=""),
                statsf[i,2], "             ",
                pvalsf[i,2], "           ",
                paste(as.character(rep(" ",13-nchar(bonmat[i,2]) )),collapse=""), 
                file=outmodel,append=T,sep="")
            cat("\n  X", tempmodel[[i]][2],
                paste(as.character(rep(" ",22-nchar(tempmodel[[i]][2]) - nchar(statsf[i,3]) )),collapse=""),
                statsf[i,3], "             ",
                pvalsf[i,3], "           ",
                paste(as.character(rep(" ",13-nchar(bonmat[i,3]) )),collapse=""), 
                file=outmodel,append=T,sep="")
            cat("\n  X", tempmodel[[i]][1], "*X", tempmodel[[i]][2],
                paste(as.character(
                   rep(" ",20-nchar(tempmodel[[i]][1])-nchar(tempmodel[[i]][2]) - nchar(statsf[i,4]) )),collapse=""),
                statsf[i,4], "             ",
                pvalsf[i,4], "           ",
                paste(as.character(rep(" ",13-nchar(bonmat[i,4]) )),collapse=""), 
                file=outmodel,append=T,sep="")
            cat("\n  omnibus           ",
                statsf[i,5], "             ",
                pvalsf[i,5], "           ",
                bonmat[i,5],
                paste(as.character(rep(" ",13-nchar(bonmat[i,5]) )),collapse=""), 
                fdrmat[i,5],
                file=outmodel,append=T,sep="")

            if (ncov > 0)
            {
               for (j in 6:(6+ncov-1))
               {
                  cat("\n  ", names(covmat)[(j-5)],
                      paste(as.character(rep(" ",23-nchar(names(covmat)[(j-5)]) - nchar(statsf[i,j]) )),collapse=""),
                      statsf[i,j], "             ",
                      pvalsf[i,j], "           ",
                      paste(as.character(rep(" ",13-nchar(bonmat[i,j]) )),collapse=""),
                      file=outmodel,append=T,sep="")
               }  
            }
         }
      }
   }

   
   cat("\n\n\n\n--------------------------------------", file=outmodel, append=T)
   cat("\nSummary Statistics: Response Variable", file=outmodel, append=T)   
   cat("\n--------------------------------------", file=outmodel, append=T)
   
   if (fam == 0)
   {
      outlevels  <- as.numeric(attr(as.factor(outvect),"levels"))
      sumstatout <- c(sum(outvect==outlevels[1]), sum(outvect==outlevels[2]))

      cat("\n\nLevel : Count\n-------------",  file=outmodel, append=T)
      cat("\n   ", outlevels[1], ": ", sumstatout[1], file=outmodel, append=T)
      cat("\n   ", outlevels[2], ": ", sumstatout[2], file=outmodel, append=T)
   }   


   if (fam == 1)
   {
      outtemp  <- c(mean(outvect, na.rm=T),
                    sd(outvect, na.rm=T),
                    min(outvect, na.rm=T),
                    quantile(outvect, c(0.25,0.5,0.75), na.rm=T),
                    max(outvect, na.rm=T))
      outtempf <- formatC(outtemp, digits=3, format="f")

      cat("\n\n  Mean:", outtempf[1],
            "\n  SD:  ", outtempf[2],
            "\n  Min: ", outtempf[3],
            "\n  Q1:  ", outtempf[4],
            "\n  Med: ", outtempf[5],
            "\n  Q3:  ", outtempf[6],
            "\n  Max: ", outtempf[7],
          file=outmodel, append=T)
   }

   cat("\n\n\n\n--------------------------------------", file=outmodel, append=T)
   cat("\nSummary Statistics: Loci", file=outmodel, append=T)   
   cat("\n--------------------------------------", file=outmodel, append=T)
   
   sumstatloci <- matrix(0, ncol=5, nrow=dim(locimat)[2])

   sumstatloci[,1] <- seq(from=1, to=dim(locimat)[2], by=1)
   sumstatloci[,2] <- colSums(locimat==0,na.rm=T)
   sumstatloci[,3] <- colSums(locimat==1,na.rm=T)
   sumstatloci[,4] <- colSums(locimat==2,na.rm=T)
   sumstatloci[,5] <- dim(locimat)[1] - rowSums(sumstatloci[,2:4])


   cat("\n\n  Loci        11        12        22     Missing", file=outmodel, append=T)
   cat("\n--------------------------------------------------", file=outmodel, append=T)
   for (i in 1:dim(sumstatloci)[1])
   {
      cat("\n  X", sumstatloci[i,1],
          paste(as.character(rep(" ",13-nchar(sumstatloci[i,1]) - nchar(sumstatloci[i,2]) )),collapse=""),
          sumstatloci[i,2],
          paste(as.character(rep(" ",10-nchar(sumstatloci[i,3]))),collapse=""),
          sumstatloci[i,3],
          paste(as.character(rep(" ",10-nchar(sumstatloci[i,4]))),collapse=""),
          sumstatloci[i,4],
          paste(as.character(rep(" ",9-nchar(sumstatloci[i,5]))),collapse=""),
          sumstatloci[i,5],
          file=outmodel, append=T, sep="")
   }
   
   if (covflag == 1)
   {
      cat("\n\n\n\n--------------------------------------", file=outmodel, append=T)
      cat("\nSummary Statistics: Covariates", file=outmodel, append=T)   
      cat("\n--------------------------------------", file=outmodel, append=T)
   
      sumstatcov <- matrix(NA, ncol=7, nrow=ncov)

      for (i in 1:ncov)
      {
         sumstatcov[i,1]        <- mean(covmatu[,i], na.rm=T)
         sumstatcov[i,2]        <- sd(covmatu[,i], na.rm=T)
         sumstatcov[i,3]        <- min(covmatu[,i], na.rm=T)
         sumstatcov[i,c(4,5,6)] <- quantile(covmatu[,i], c(0.25,0.5,0.75), na.rm=T)
         sumstatcov[i,7]        <- max(covmatu[,i], na.rm=T)
      }

      sumstatcov <- formatC(sumstatcov, digits=3, format="f")

      for (i in 1:dim(sumstatcov)[1])
      {
         cat("\n\n", names(covmat)[i], "\n---------------", file=outmodel, append=T, sep="")
         cat("\n\n  Mean:", sumstatcov[i,1],
               "\n  SD:  ", sumstatcov[i,2],
               "\n  Min: ", sumstatcov[i,3],
               "\n  Q1:  ", sumstatcov[i,4],
               "\n  Med: ", sumstatcov[i,5],
               "\n  Q3:  ", sumstatcov[i,6],
               "\n  Max: ", sumstatcov[i,7],
             file=outmodel, append=T)
      }
   }

   
   if ((nullflag == 1) & (is.na(outnull) == F))
   {
      # tempnullint  <- strsplit(c(unlist(strsplit(tnullint[,1], ","))), ":")

      cat("\n--------------------------------------------------------------------------", file=outnull, append=F)
      cat("\nAssociation Studies for Quantitative Traits in Structured Populations:", file=outnull, append=T)
      cat("\nResults from Models Using Null Loci", file=outnull, append=T)
      cat("\n\nDate:",date(), file=outnull, append=T)
      cat("\n--------------------------------------------------------------------------", file=outnull, append=T)
      cat("\n\n--------------------------------------------------------------------------", file=outnull, append=T)
      cat("\nProgram: GCF", file=outnull, append=T)
      cat("\nResponse Type (Y):", famstring, file=outnull, append=T)
      if (covflag == 0)
         cat("\nCovariates: No", file=outnull, append=T)
      if (covflag == 1)
         cat("\nCovariates: Yes", file=outnull, append=T)
      cat("\nNumber of loci in data file:", dim(locimattemp)[2]/2, file=outnull, append=T)
      cat("\nNumber of subjects:", dim(locimatc)[1], file=outnull, append=T)
      cat("\nNull Loci:", nullloci, file=outnull, append=T)
      cat("\n--------------------------------------------------------------------------", file=outnull, append=T)

      cat("\n\n---------------\nMarginal models\n---------------", file=outnull, append=T)
      cat("\n\n  Model         Statistic       p-value", file=outnull, append=T)
      cat("\n--------------------------------------------------------------------------", file=outnull, append=T)

      for (i in 1:dim(tnullmarg)[1])
      {
         cat("\n  X", tnullmarg[i,1],
             paste(as.character(rep(" ",20-nchar(tnullmarg[i,1]) - nchar( tnullmargfs[i]) )),collapse=""),
             tnullmargfs[i], "         ",
             tnullmargfp[i],
             file=outnull, append=T, sep="")
      }
      
      if (margint == 4)
      {
         cat("\n\n\n------------------\nInteraction models\n------------------", file=outnull, append=T)
         cat("\n\n                      Interaction                     Omnibus", file=outnull, append=T)
         cat("\n                -----------------------       -----------------------", file=outnull, append=T)
         cat("\n  Model         Statistic       p-value       Statistic       p-value", file=outnull, append=T)
         cat("\n--------------------------------------------------------------------------", file=outnull, append=T)

         for (i in 1:dim(tnullint)[1])
         {
            temp <- strsplit(tnullint[i,1],":")[[1]]
            if (tnullint[i,2] == -9)
               cat("\n  X", temp[1], "*X", temp[2],
                   "     A problem occured when fitting this null model.",
                   file=outnull, append=T, sep="")
            else
               cat("\n  X", temp[1], "*X", temp[2],
                   paste(as.character(rep(" ",18-nchar(temp[1])-nchar(temp[2]) - nchar(tnullintfis[i]) )),collapse=""),
                   tnullintfis[i], "         ",
                   tnullintfip[i],
                   paste(as.character(rep(" ",15-nchar(tnullintfos[i] ) )),collapse=""),
                   tnullintfos[i], "         ",
                   tnullintfop[i],
                   file=outnull, append=T, sep="")
         }
      }
   }
   
   write.table(data.frame(statsf[,2], pvalsf[,2]), file="GCTests.txt", sep="\t", quote=F, row.names=F)
   cat("\nStatus: Done\n\n")
   ## options(warn=warntemp)
   
}

Main.Func()



##############################################################



