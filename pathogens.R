setwd("C:/Users/USUARIO/Documents/cietrop/analisis_2.2_250724")
dir()
source("metagenomic_commands_2.R")

##################
### BACTERIOMA ###
##################

## hisopado ##

baclist <- read.csv("bacteria.txt", header=T)
head(baclist)
v <- sort(unique(baclist[,1]))

bact <- read.csv("hisopado_bactSTND/Bacterioma_total_species_0.csv", header=T)
head(bact)

a1 <- 0 ; a2 <- 0 ; 
for (i in 1:length(v)){
a1 <- grep(v[i], bact$species, value=T)
a2 <- append(a2,a1)
}

a3 <- sort(unique(a2)[2:length(a2)])

bact_pt <- bact[bact$species %in% a3, ]
head(bact_pt)

metaheat(bact_pt, 100,"Bacterias patógenas (hisopado)","Bacteria")

## suero ##

bact <- read.csv("suero_bactSTND/Bacterioma_total_suero.csv", header=T)
head(bact)

a1 <- 0 ; a2 <- 0 ; 
for (i in 1:length(v)){
a1 <- grep(v[i], bact$species, value=T)
a2 <- append(a2,a1)
}

a3 <- sort(unique(a2)[2:length(a2)])

bact_pt <- bact[bact$species %in% a3, ]
head(bact_pt)

metaheat(bact_pt, 100,"Bacterias patógenas (suero)","Bacteria")

##############
### VIROMA ###
##############

## hisopado ##

virlist <- read.csv("virus.txt", header=T)
head(virlist)
v <- sort(unique(virlist[,1]))

vir <- read.csv("hisopado_vir/Viroma_total_hispoado.csv", header=T)
head(vir)

a1 <- 0 ; a2 <- 0 ; 
for (i in 1:length(v)){
a1 <- grep(v[i], vir$species, value=T)
a2 <- append(a2,a1)
}

a3 <- sort(unique(a2)[2:length(a2)])

vir_pt <- vir[vir$species %in% a3, ]
head(vir_pt)

metaheat(vir_pt, 100,"Virus patógenos (hisopado)","Virus")

## suero ##

vir <- read.csv("suero_vir/Viroma_total_suero.csv", header=T)
head(vir)

a1 <- 0 ; a2 <- 0 ; 
for (i in 1:length(v)){
a1 <- grep(v[i], vir$species, value=T)
a2 <- append(a2,a1)
}

a3 <- sort(unique(a2)[2:length(a2)])

vir_pt <- vir[vir$species %in% a3, ]
head(vir_pt)

metaheat(vir_pt, 100,"Virus patógenos (suero)","Virus")