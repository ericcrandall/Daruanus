## This is free software written by a free software
## covered by the CeCILL licence (GPL compatible).
## © F. Rousset & R. Leblois 2007- ; contributor: C.R.Beeravolu
## Migraine version 0.6 (Built on Apr 14 2021 at 21:34:18).
##  rosglobal$latt2Ns2 : estimate of 2{Nh}m condS2 where Nh is haploid number per lattice step and condS2 is is lattice units
## ...

if (interactive()) {options(error=recover)} else {
  options(echo = FALSE)
  options(error = quote(dump.frames(paste("dump",blackbox.getOption("jobSampleNbr"),sep=""), TRUE)))
}
rm(list=ls())
chk <- try(library(blackbox))
if (class(chk)=="try-error") stop("(!) Failed to load the 'blackbox' library")
chk1 <- try(library(doSNOW))
if (class(chk1)=="try-error") message("(!) Failed to load the 'doSNOW' library, parallelization will use 'pbapply()' instead ")
GP<-list()
GP$usedBy <- "Migraine"
GP$jobSampleNbr<-1
GP$coreNbr<-4
GP$dataFile <- paste("pointls_",GP$jobSampleNbr,".txt",sep="")
GP$estimOutf <- file(paste("output_",GP$jobSampleNbr,".txt",sep=""), "w")
GP$nextPointsf <- paste("nextpoints_",GP$jobSampleNbr,".txt",sep="")
GP$cleanResu <- file(paste("results_",GP$jobSampleNbr,".txt",sep=""), "w")
GP$plotFiles<-list() ## sinon le premier appel de 'plotFiles[[filename]]<-' cree un vecteur de lgr 1 et le 2 appel plante
GP$ParameterNames<-c("twoNmu","twoNm","g") ## Canonical order of elements in parameter vector
GP$paramnbr <- length(GP$ParameterNames) ## fitted or not
GP$timeScale <- "popSize" 
GP$graphicsFormat<-"pdf"
GP$miscOptions<-c("optimizeKriging","ignoreSmallEigen")
GP$verbosity <- 1 ## display of info about GCV
GP$DemographicModel <- "IBD" 
GP$subsetRows <- NULL ## no selection of rows in pointls
GP$FONKgNames <- c("twoNmu","twoNm","condS2")
GP$FONKgScale <- c("Canonical","Canonical","logscale")
names(GP$FONKgScale) <- GP$FONKgNames
GP$samplingSpace<-c("twoNmu","twoNm","condS2") ## ONLY for next points
GP$samplingScale<-c("Canonical","Canonical","logscale") ## ONLY for next points
names(GP$samplingScale)<-GP$samplingSpace
GP$extraScale<-c(latt2Ns2="logscale") ## 
GP$profile3passesBool<-TRUE ##  should become user option ?
GP$ICthresh<-1e8 ## default threshold for ignoreSmallEigen
## (!) plotRange can contain Nb rather than latt2Ns2
GP$plotRange<-list()
GP$graphicPars<-list()
GP$plotOptions<-c("All1DProfiles","NoSurfaces")
GP$oneDimCIvars<-c("twoNmu","twoNm","g") ## (an empty c() value will be meaningful for nextBounts)
GP$spec1DProfiles<-c()
GP$spec2DProfiles<-c()
GP$S2factor<-1 ##  other values add more confusion that anything else
GP$D2IBDbool<-FALSE;GP$D1IBDbool<-TRUE 
GP$Nbfactor<-58 ##  correction from 2NmS2 to 2Ds2=Nb
GP$GeoUnit<-"ind.km" ## string that gives the Nb unit in plot
GP$mincondS2<-1
if ( ! is.null(GP$plotRange[["Nb"]])) GP$plotRange[["latt2Ns2"]]<-GP$plotRange[["Nb"]]/GP$Nbfactor
GP$fittedLoci <- NULL
## Booleans: 
GP$LikStatistic<-"PAC" ## PAC/IS...
GP$designRetain<-1 ## fraction of points in the Top retained for Kriging, removing occurrences of close points 
GP$GCVdesignRetain<-1 ## fraction of points in the Top retained for GCV, removing occurrences of close points 
GP$GCVptnbr<--1 ## nbr of optimization steps for GCV
GP$GCVoptimSteps<-100 ## nbr of optimization steps for GCV
GP$covFamily<-"Matern" ## !!!!! covariance function for kriging
GP$krigmax<--1
GP$kriglength<-2000## will be overwritten if krigmax<0
GP$krigoverlap<-500## will be overwritten if krigmax<0
GP$minKrigPtNbr<--1 ## nbr of optimization steps for GCV
GP$maxKrigPtNbr<-FALSE ## max pts for Kriging
## larger bounds below => fewer boundary gcv problems but not necessarily meaningful
GP$GCVlowerFactor<-20 ## 2 *w*as suggested by Martin & Simpson Am. Inst. Aeron. Astron. 43: 853- (2005)
GP$GCVupperFactor<-5 ## as suggested by ibidem
if (GP$covFamily %in% c("Matern", "GneitingMatern")) {
	GP$minSmoothness<-2
	GP$maxSmoothness<-4
}
GP$CovFnParam <- numeric(0);GP$CovFnParamInSettingsBool<-FALSE
GP$metarange<-3 ## controls covariance range, except if CovFnParam vector is estimated by GCV 
GP$hullExpandFactor<-2
GP$profileEdges<-c() ## will accumulate masking points encountered during profile computations
GP$gridstepsNbr<-21
GP$LRTlist<-list()
GP$CIlevel<-0.05
GP$NextBoundsLevel<-0.001
GP$nextBounds<-"pointsFromR"
GP$nextPointNumber<-250
GP$ptSamplingSeed<-67144630 ## random seed for sampling new points
GP$parDigits<-7 ## digits of parameter coordinates in files
GP$scalefactor<-1 ## for (fnscale,parscale in optim())

write("Migraine 0.6 (Built on Apr 14 2021 at 21:34:18)",file=GP$cleanResu)

write(paste(paste(packageDescription("blackbox",fields = c("Package", "Version")), collapse=", version "),"loaded"),file=GP$cleanResu)
write(paste("R code run on ",date()),file=GP$cleanResu)
write("",file=GP$cleanResu)
write("Data file: ../../Daruanus_FP.txt",file=GP$cleanResu)
write("Settings file: migraine.txt",file=GP$cleanResu)
write("",file=GP$cleanResu)
GP$testPointList <- NULL
# GP$interactiveGraphics <- FALSE # Uncomment this to switch to file graphic output in interactive session
GP <- preprocessbboptions(optionList=GP)
do.call(blackbox.options, GP)
if (!blackbox.getOption("interactiveGraphics"))  {providePlotFile(blackbox.getOption("basicRplotsfile"))}
## *************************
## **** Processing data ****
## *************************

pointls <- buildPointls(dataFile=blackbox.getOption("dataFile"), 
                        respCols=NULL,
                        subsetRows=NULL,
                        ycolname="-ln(L)", 
                        cleanResu=blackbox.getOption("cleanResu"))

FONKgpointls <- buildFONKgpointls(pointls=pointls) # includes blackbox.options(FONKgpointls = FONKgpointls)

GCVblob <- calcGCV(sorted_data=FONKgpointls, 
                   CovFnParam=blackbox.getOption("CovFnParam"),
                   GCVptnbr=blackbox.getOption("GCVptnbr"), 
                   topmode="dLnL",
                   verbose=blackbox.getOption("verbosity"), 
                   cleanResu=blackbox.getOption("cleanResu"), 
                   optimizers=blackbox.getOption("optimizers") 
                   )

provideDevice(bbdefaultPars=TRUE) ## one window for _each_ graphic

fitobject <- calcPredictorOK(FONKgpointls=FONKgpointls, 
                   minKrigPtNbr=blackbox.getOption("minKrigPtNbr"), 
                   krigmax=blackbox.getOption("krigmax"), 
                   topmode="dLnL", 
                   rawPlots=TRUE, 
                   cleanResu=blackbox.getOption("cleanResu")) ## includes rawProfiles


primarymax <- maximizeOK(fitobject=blackbox.getOption("fitobject"), cleanResu=blackbox.getOption("cleanResu"))
calc1Dprofiles()


calc2D3Dplots() ## **** plain or slice plots of predicted likelihood ****

## **** More graphics : 2D profiles ****
calcProfileLR(varNames=blackbox.getOption("fittedNames"), cleanResu=blackbox.getOption("cleanResu")) 

## regression observed/ kriging prediction
provideDevice(bbdefaultPars=TRUE) ## one window for _each_ graphic
plot(blackbox.getOption("fitobject"), which=1) ## Diagnostic plot lnL Estimates vs lnL predicted values

## **** Likelihood-ratio tests ****
LRTs <- calcLRTs(testPointList=blackbox.getOption("testPointList"),
         cleanResu=blackbox.getOption("cleanResu"))

## **** 1D CIs ****
if (length(blackbox.getOption("oneDimCIvars"))>0 
           || blackbox.getOption("nextBounds") %innc% c("from1dci", "frompoints")) {
  CIpointsList <- calc1DCIs(oneDimCIvars=blackbox.getOption("oneDimCIvars"), 
                     FONKgNames=blackbox.getOption("FONKgNames"), 
                     fittedNames=blackbox.getOption("fittedNames"), 
                     CIlevel=blackbox.getOption("CIlevel"), 
                     nextBounds=blackbox.getOption("nextBounds"), 
                     boundsOutfile=paste("bounds_", blackbox.getOption("jobSampleNbr"), ".txt", sep=""), 
                     dataString=blackbox.getOption("dataName"), 
                     cleanResu=blackbox.getOption("cleanResu")
  ) 
  blackbox.options(CIpointsList=CIpointsList) ## important for not yet transparent code
} 

## **** Sampling design for next iteration ****
if (blackbox.getOption("nextBounds") %==nc% "pointsfromr") {
  goodpoints <- sampleByResp(size=blackbox.getOption("nextPointNumber"), 
                  outfile = blackbox.getOption("nextPointsf"), 
                  useEI=blackbox.getOption("useEI"), 
                  rnd.seed=blackbox.getOption("ptSamplingSeed")
  )
}  

writeFinalInfo(cleanResu=blackbox.getOption("cleanResu"))
# close(blackbox.getOption("estimOutf")) ## ! leaves connection open if interactive run !

