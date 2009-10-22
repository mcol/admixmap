#place headers in bclib dir
includedir=$(prefix)/include/bclib

#headers to be installed
top_headers = bclib.h common.h

sampler_headers = AdaptiveRejection.h DirichletParamSampler.h DispersionSampler.h Gaussian.h GaussianProposalMH.h HamiltonianMonteCarlo.h ModeFinder.h MuSampler.h rand.h StepSizeTuner.h

regression_headers = Regression.h LogisticRegression.h LinearRegression.h CoxRegression.h

utils_headers = ColumnIter.h DataMatrix.h DataReader.h dist.h DelimitedBuf.h Exceptions.h Delimitedostream.h DelimitedFileWriter.h GSLExceptions.h GSLErrorHandler.h linalg.h misc.h LogWriter.h Array.hh OptionReader.h ProbVector.h pvector.h ranker.h RObjectWriter.h StringSplitter.h StringConvertor.h exceptions.h estr.h NRand.h

#headers to be installed
include_HEADERS = $(top_headers) $(sampler_headers) $(regression_headers) $(utils_headers)
