setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(ape)
d <- read.nexus.data("Exp0187_2021-07-18_IE-CoR_Lgs161_Mgs170_BEAUti_sorted.nex")
wn <- read.table("word_names.txt")
wn <- wn[,2]
wp <- read.table("words_positions.txt")
n <- length(d)

nstates <- vector()
for(j in 1:nrow(wp)){
  # if(j %in% keep){
  concept <- wp[j,1]
  v <- rep("?", n)
  for(i in 1:n){
    wn2<- wn[wp[j,2]:wp[j,3]]
    x <- wn2[which(d[[i]][wp[j,2]:wp[j,3]]=="1")]
    if(length(x)>0){
      v[i]<- x[1]
    }
  }
  m <- length(unique(v)[! unique(v) %in% "?"])
  nstates <- append(nstates, m+1)
}

nstates <- nstates[order(nstates)]
frequencies <- 1/nstates

power_mrates <- function(nstates, power){
  nstatesp <- nstates^power
  mratesp <- nstatesp/mean(nstatesp)
  return(mratesp)
}

library(viridis)

pdf(file="mutationrates.pdf")
plot(x=nstates, y=power_mrates(nstates, 2), pch=20, col="#440154FF", ylab="Mutation rate", xlab="Number of observed cognate sets", main="Mutation rates")
points(x=nstates, y=power_mrates(nstates, 1.8), pch=20, col="#443A83FF")
points(x=nstates, y=power_mrates(nstates, 1.6), pch=20, col="#31688EFF")
points(x=nstates, y=power_mrates(nstates, 1.4), pch=20, col="#21908CFF")
points(x=nstates, y=power_mrates(nstates, 1.2), pch=20, col="#35B779FF")
points(x=nstates, y=power_mrates(nstates, 1.0), pch=20, col="#8FD744FF")
points(x=nstates, y=power_mrates(nstates, 0.0), pch=20, col="#FDE725FF")
legend(x = "topleft",  
       legend = c("2.0", "1.8", "1.6", "1.4", "1.2", "1.0", "0.0"), 
       fill = viridis(7),
       title="spread") 
dev.off()

pdf(file="Instantaneousrates.pdf")
plot(x=nstates, y=power_mrates(nstates, 2)*frequencies, pch=20, col="#440154FF", ylab="Instantaneous transition rate", xlab="Number of observed cognate sets", main="Instantaneous transition rates")
points(x=nstates, y=power_mrates(nstates, 1.8)*frequencies, pch=20, col="#443A83FF")
points(x=nstates, y=power_mrates(nstates, 1.6)*frequencies, pch=20, col="#31688EFF")
points(x=nstates, y=power_mrates(nstates, 1.4)*frequencies, pch=20, col="#21908CFF")
points(x=nstates, y=power_mrates(nstates, 1.2)*frequencies, pch=20, col="#35B779FF")
points(x=nstates, y=power_mrates(nstates, 1.0)*frequencies, pch=20, col="#8FD744FF")
points(x=nstates, y=power_mrates(nstates, 0.0)*frequencies, pch=20, col="#FDE725FF")
legend(x = "topleft",  
       legend = c("2.0", "1.8", "1.6", "1.4", "1.2", "1.0", "0.0"), 
       fill = viridis(7),
       title="spread") 
dev.off()

pdf(file="Instantaneousrates_fullscale.pdf")
plot(x=nstates, y=power_mrates(nstates, 2)*frequencies, pch=20, col="#440154FF", ylab="Instantaneous transition rate", xlab="Number of observed cognate sets", main="Instantaneous transition rates", ylim=c(0, 0.5))
points(x=nstates, y=power_mrates(nstates, 1.8)*frequencies, pch=20, col="#443A83FF")
points(x=nstates, y=power_mrates(nstates, 1.6)*frequencies, pch=20, col="#31688EFF")
points(x=nstates, y=power_mrates(nstates, 1.4)*frequencies, pch=20, col="#21908CFF")
points(x=nstates, y=power_mrates(nstates, 1.2)*frequencies, pch=20, col="#35B779FF")
points(x=nstates, y=power_mrates(nstates, 1.0)*frequencies, pch=20, col="#8FD744FF")
points(x=nstates, y=power_mrates(nstates, 0.0)*frequencies, pch=20, col="#FDE725FF")
legend(x = "topright",  
       legend = c("2.0", "1.8", "1.6", "1.4", "1.2", "1.0", "0.0"), 
       fill = viridis(7),
       title="spread") 
dev.off()

pdf(file="Instantaneousrates_estimation.pdf")
plot(x=nstates, y=power_mrates(nstates, 2)*frequencies, pch=20, col="#44015410", ylab="Instantaneous transition rate", xlab="Number of observed cognate sets", main="Instantaneous transition rates")
points(x=nstates, y=power_mrates(nstates, 1.8)*frequencies, pch=20, col="#443A8310")
points(x=nstates, y=power_mrates(nstates, 1.6)*frequencies, pch=20, col="#31688E10")
points(x=nstates, y=power_mrates(nstates, 1.4)*frequencies, pch=20, col="#21908C10")
points(x=nstates, y=power_mrates(nstates, 1.2)*frequencies, pch=20, col="#35B77910")
points(x=nstates, y=power_mrates(nstates, 1.0)*frequencies, pch=20, col="#8FD74410")
points(x=nstates, y=power_mrates(nstates, 0.0)*frequencies, pch=20, col="#FDE72510")
points(x=nstates, y=power_mrates(nstates, 1.367)*frequencies, pch=18, col="#E3311D")

legend(x = "topleft",  
       legend = c("2.0", "1.8", "1.6", "1.4", "1.2", "1.0", "0.0"), 
       fill = viridis(7),
       title="spread") 
dev.off()


