#To make a genind object, I need to combine every two columns together for a standard msat table input.
#More difficult than it looks! Used code from here:
#https://stackoverflow.com/questions/28686848/speedy-elegant-way-to-unite-many-pairs-of-columns
mash_alleles <- function(df,first_col) {
  mapply(function(x, y) {
    paste(x, y, sep = "/")},
    df[ ,seq(first_col, ncol(df), by = 2)],
    df[ ,seq(first_col+1, ncol(df), by = 2)])
}



# two functions following the two-trait equations in Waples et al. 2014 Table 3
nbadj2 <- function(Nb, Adult.Lifespan = Adult.Lifespan , Alpha.maturity = Alpha.maturity){
  Nbadj <- Nb / (1.103 - (0.245 * log10(Adult.Lifespan/Alpha.maturity)))
  return(Nbadj)
}

neadj2 <- function(Nbadj, Adult.Lifespan = Adult.Lifespan , Alpha.maturity = Alpha.maturity){
  Neadj <- Nbadj / (0.485 + (0.758 * log10(Adult.Lifespan/Alpha.maturity)))
  return(Neadj)
}


neadj_chr <- function(Neadj2, oneNchromosome_number){
  # a third correction for chromosome linkage based on Equation 1a of Waples et al. 2016.
  #Returns the bias estimate and the correction in a vector
  bias <- 1 - (0.098 + 0.219*log(oneNchromosome_number))
  neadjchr <- Neadj2 + bias*Neadj2
  return(c(bias=bias,correctedNe=neadjchr))
}

WDFilter <- function(Nb_estimates, minSampSize){
  # Filter Ne_Estimator output following the "rule of thumb" of Waples and Do 2010, pg 251
  # Extended to include sample sizes greater than 100, as described by the second paragraph on that page.
  # Ne Estimates should be run for a range of critical values e.g. 0.01 - 0.05 by 0.01
  # Columns of the input dataframe need to be titled SampSize and Pcrit
  # minSampSize sets a minimum sample size for results from that population to be included in the output
  Nb_estimates %>% filter(SampSize >= 100 & Pcrit == 0.01| SampSize < 100 & SampSize >= 25  & Pcrit == 0.02 | 
                          SampSize < 25 & Pcrit >= 1/(2*SampSize) ) %>%
  group_by(Population) %>% filter(Pcrit == min(Pcrit)) %>% filter(SampSize >= minSampSize)
}



# Equations 3,4,5 from Supplemental Methods of Pinsky et al. 2017
eq3Pinsky <- function(Adult.Lifespan,Alpha.maturity){
  nb_over_ne <- 0.485 + 0.758*(log10(Adult.Lifespan/Alpha.maturity))
}

eq4Pinsky <- function(Nb, nb_over_ne){
  Nbadj <- Nb/(1.26+0.323*nb_over_ne)
}

eq5Pinsky <- function(Nbadj,nb_over_ne){
  Neadj <- Nbadj/nb_over_ne
}


g_to_sigma2 <- function(g, twoD = F){
# This converts estimates of g to estimates of sigma^2 according to equations 1 & 2 in
# Rousset and Leblois 2012. twoD is a boolean (T/F) for whether to convert g for
# a 1 or 2 dimensional IBD model (default is 1 dimension).
# This was originally an internal function from Rousset and Leblois (2007) blackbox R package.
  if(twoD) {
    return((1 + g)/((2 - g) * (1 - g)^2))
  }
  else return((1 + g)/((1 - g)^2))
}

#Literature Cited
#Pinsky, M. L., Saenz-Agudelo, P., Salles, O. C., Almany, G. R., Bode, M., Berumen, M. L., … Planes, S. (2017). Marine Dispersal Scales Are Congruent over Evolutionary and Ecological Time. Current Biology : CB, 1–16. doi: 10.1016/j.cub.2016.10.053
#Rousset, F., & Leblois, R. (2007). Likelihood and Approximate Likelihood Analyses of Genetic Structure in a Linear Habitat: Performance and Robustness to Model Mis-Specification. Molecular Biology And Evolution, 24(12), 2730–2745. doi: 10/d3vn37
#Rousset, Francois, & Leblois, R. (2012). Likelihood-based inferences under isolation by distance: Two-dimensional habitats and confidence intervals. Molecular Biology And Evolution, 29(3), 957–973. doi: 10/fgs233
#Waples, R. S., & Do, C. (2010). Linkage disequilibrium estimates of contemporary N-e using highly variable genetic markers: A largely untapped resource for applied conservation and evolution. Evolutionary Applications, 3(3), 244–262. doi: 10.1111/j.1752-4571.2009.00104.x
#Waples, R. S., Antao, T., & Luikart, G. (2014). Effects of overlapping generations on linkage disequilibrium estimates of effective population size. Genetics, 197(2), 769–780. doi: 10.1534/genetics.114.164822
#Waples, R. K., Larson, W. A., & Waples, R. S. (2016). Estimating contemporary effective population size in non-model species using linkage disequilibrium across thousands of loci. Heredity, 117(4), 233–240. doi: 10/f9d8kn


read_migraine <- function(runDir,outfileName="results_1.txt"){
  outfile <- scan(file=file.path(runDir,outfileName),what="character",sep="\n") 
  lattice2geog <- as.numeric(strsplit(grep("bin width",outfile,value=T),
                                      "=", perl = T)[[1]][2])
  NS <- as.numeric(strsplit(grep("Neighborhood",outfile,value=T)," ")[[1]][8])
  NSCI <- as.numeric(strsplit(grep("interval for Nb",outfile,value=T)," ")[[1]][c(8,10)])
  pointEstimatesline <- grep("Point estimates", outfile) +2
  Nmu <- as.numeric(strsplit(outfile[pointEstimatesline]," +")[[1]][2])
  Nm <- as.numeric(strsplit(outfile[pointEstimatesline]," +")[[1]][3])
  g <- as.numeric(strsplit(outfile[pointEstimatesline]," +")[[1]][4])
  result <- c("lattice2geog"=lattice2geog, "NS"=NS, "NSCI" = NSCI, "Nmu"=Nmu, "Nm" = Nm, "g"= g)
}


read_NeEstimator <- function(file) {
# read in table output from the NeEstimator program

  col.names <- c("Population","SampSize","Pcrit","WeightedMean","IndAlleles",
                 "r2","ExpR2","Ne","ParametricLow","ParametricHigh",
                 "JackknifeLow","JackknifeHigh","EffDF")
  col.types = "cdddddddddddd"
  output <- read_tsv(file = file, skip = 14, 
                 col_names = col.names, col_types = col.types) %>% 
                fill(Population) %>% fill(SampSize) %>% filter(!is.na(Ne))
  
}

r2p <- function(r2, expr2){
  #r2 is r^2 and expr2 is E(r^2_sample)
  # see footnote to Table 2 in Waples 2006 10.1007/s10592-005-9100-y
  r2p = r2 - expr2
  return(r2p)
}
WaplesMonoNe <- function(r2p){
  #based on the "Monogamy" entry of Table 2 in Waples 2006 10.1007/s10592-005-9100-y
  #r2p is r^2 prime as given in the footnote to that table
  Ne = ((2/3)+sqrt((4/9) - (7.2*r2p))) / (2*r2p)
  return(Ne)
  }
