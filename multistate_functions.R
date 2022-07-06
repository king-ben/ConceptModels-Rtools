library(tidyverse)
library(ape)
library(XML)
library(devtools)
source_url("https://raw.githubusercontent.com/king-ben/R-functions/master/asr_functions.R")


multistate_xml <- function(file){
  wp <- find_partitions(file)
  wn <- find_site_names(file)
  d <- read.nexus.data(file)
  n <- length(d)
  #get mutation rates
  mr <- wp[,3]-wp[,2] + 1
  mr <- mr/mean(mr)
  
  likelihood <- newXMLNode("distribution", attrs = c(id="likelihood", spec="util.CompoundDistribution", useThreads="true"))
  mutationratesoperator <- newXMLNode("operator", attrs = c(id="mutationRatesOperator", spec="MutationRateReweight", parameter="@mutationrate_spread", weight="3.0"))
  logger <- newXMLNode("logger")
  codemaplist <- list()
  for(j in 1:nrow(wp)){
    concept <- wp[j,1]
    v <- rep("?", n)
    wn2<- wn[wp[j,2]:wp[j,3]]
    namb <- 0
    ambiguities <- list()
    for(i in 1:n){
      x <- wn2[which(d[[i]][wp[j,2]:wp[j,3]]=="1")]
      if(length(x)==1){
        v[i]<- x[1]
      }
      if(length(x)>1){
        namb = namb+1
        v[i] <- paste0("ambiguity", namb)
        ambiguities[[namb]] <- x
      }
    }
    
    
    wn3 <- wn2[2:length(wn2)]
    treelikelihood <- newXMLNode("distribution", parent=likelihood, attrs = c(id=paste0("treeLikelihood.", concept), spec="TreeLikelihood", tree="@Tree.t:tree", branchRateModel="@clock", useAmbiguities="true"))
    data <- newXMLNode("data", parent=treelikelihood, attrs = c(id=paste0("data.", concept), spec="AlignmentFromTrait"))
    traitset <- newXMLNode("traitSet", parent=data, attrs = c(id=paste0("traitSet.", concept), spec="beast.evolution.tree.TraitSet", taxa="@TaxonSet", traitname="discrete"))
    z <- paste0(names(d), "=", v, ",")
    z[n] <- gsub(",", "", z[n])
    xmlValue(traitset) <- paste0(z, collapse=" ")
    m <- length(wn3)
    codemap <- paste0(wn3, "=", 0:(m-1))
    codemaplist[[j]] <- data.frame(code=0:(m-1), cognate=wn3)
    unknown <- paste("?", "=", paste(0:(m), collapse=" "))
    unobserved <- paste0("unobserved", "=", m)
    codemap2 <- paste(paste(codemap, collapse=","), unobserved, unknown, sep=",")
    if(length(ambiguities)>0){
      for(i in 1:length(ambiguities)){
        codemap2 <- paste(codemap2, paste0("ambiguity", i, "=", paste(match(ambiguities[[i]], wn3)-1, collapse=" ")), sep=",")
      }
    }
    
    datatype <- newXMLNode("userDataType", parent=data, attrs=c(id=paste0("traitDataType.", concept), spec="beast.evolution.datatype.UserDataType", codeMap=codemap2, codelength="-1", states=m+1))
    sitemodel <- newXMLNode("siteModel", parent=treelikelihood, attrs = c(id=paste0("siteModel.", concept), spec="SiteModel", gammaCategoryCount="1"))
    mutationrate <- newXMLNode("parameter", parent=sitemodel, attrs = c(id=paste0("mutationRate.", concept), spec="parameter.RealParameter", estimate="false", name="mutationRate"))
    xmlValue(mutationrate) <- mr[j]
    substmodel <- newXMLNode("substModel", parent=sitemodel, attrs = c(id=paste0("substModel.", concept), spec="LewisMK", datatype=paste0("@traitDataType.", concept)))
    frequencies <- newXMLNode("frequencies", parent=substmodel, attrs=c(spec="Frequencies"))
    # freqs <- newXMLNode("frequencies", parent=frequencies, attrs=c(id=paste0("frequencies.", concept), spec="parameter.RealParameter", value=paste(paste0(rep((1-unobservedfreq)/(m), m), collapse=" "), unobservedfreq)))
    freqs <- newXMLNode("frequencies", parent=frequencies, attrs=c(id=paste0("frequencies.", concept), spec="parameter.RealParameter", value=paste0(rep(1/(m+1), m+1), collapse=" ")))
    
    #make mutationrate operator
    mutrate <- newXMLNode("mutationrate", parent=mutationratesoperator, attrs=c(idref=paste0("mutationRate.", concept)))
    # nstates <- append(nstates, m+1)
    #mutationrate loggers
    log <- newXMLNode("log", parent=logger, attrs = c(idref=paste0("mutationRate.", concept)))
  }

  names(codemaplist) <- wp[,1]
  nstates <- wp[,3]-wp[,2]
  nstates_xml <- newXMLNode("nstates", parent=mutationratesoperator, attrs = c(id="nstates", spec="parameter.IntegerParameter", dimension=nrow(wp), estimate="false"))
  xmlValue(nstates_xml) <- paste(nstates, collapse=" ")

  operators <- newXMLNode("operators")
  freqlogger <- newXMLNode("logger", attrs=c(id="tracelog.frequencies", spec="Logger", fileName="trace.frequencies.log", logEvery="10000", model="@posterior", sort="smart"))
  for(j in 1:nrow(wp)){
    concept <- wp[j,1]
    freqop <- newXMLNode("operator", parent=operators, attrs=c(id=paste0("freqop.", concept), spec="UnobservedFrequencyExchanger", windowSize="0.05", parameter=paste0("@frequencies.", concept), weight="0.3"))
    freqlog <- newXMLNode("log", parent=freqlogger, attrs=c(id=paste0("unobs.", concept), spec="beast.core.util.UnobservedFrequencyLogger", parameter=paste0("@frequencies.", concept)))
  }
  
  state <- newXMLNode("state")
  for(j in 1:nrow(wp)){
    concept <- wp[j,1]
    mutationratenode <- newXMLNode("stateNode", parent=state, attrs=c(idref=paste0("mutationRate.", concept)))
    frequencynode <- newXMLNode("stateNode", parent=state, attrs=c(idref=paste0("frequencies.", concept)))
  }
  
  out <- list(likelihood, mutationratesoperator, logger, operators, freqlogger, state, codemaplist)
  names(out) <- c("likelihood", "mutationrates_operator", "mutationrates_logger", "frequency_operators", "frequency_logger", "state", "codemap")
  return(out)
}

saveall_multistate <- function(xml){
  saveXML(xml$likelihood, "likelihood.xml")
  saveXML(xml$mutationrates_operator, "mutationratesoperator.xml")
  saveXML(xml$mutationrates_logger, "mutationrateslogger.xml")
  saveXML(xml$state, "state.xml")
  saveXML(xml$frequency_operators, "frequency_operators.xml")
  saveXML(xml$frequency_logger, "frequency_logger.xml")
  codemap <- xml$codemap
  save(codemap, file="codempa.Rdata")
}
