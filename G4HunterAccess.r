#### R script for the G4Access Project
#### For JC Andrau's lab

#### map G4 with G4Hunter for SacCer3 and DM6

source("G4HunterAccess_function.r")

### DM6
library("BSgenome.Dmelanogaster.UCSC.dm6")
genome <- BSgenome.Dmelanogaster.UCSC.dm6
seqinf <- seqinfo(genome)
seqlevels(seqinf) <- seqlevels(genome)[1:8]
w=25

G4hk=mclapply((8:1),function(i) {x=mG4huntlistref(i,k=w,hl=c(1.2,1.5,2),gen=genome,with.seq=F);print(i);return(x)},mc.cores=8L)
G4H_w25=unlist(do.call(c,G4hk))
seqlevels(G4H_w25, pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(G4H_w25) <- seqinf

save(G4H_w25,file="DM6_G4H_w25.RData")
lapply(c(1.2,1.5,2), function(i) {
	export(G4H_w25[G4H_w25$hl==i], con=paste0("DM6_G4H",i,".bed"))
})

### Yeast data
library("BSgenome.Scerevisiae.UCSC.sacCer3")
genome <- BSgenome.Scerevisiae.UCSC.sacCer3
seqinf <- seqinfo(genome)
seqlevels(seqinf) <- seqlevels(genome)[1:17]
w=25

G4hk=mclapply((1:17),function(i) {x=mG4huntlistref(i,k=w,hl=c(1.2,1.5,2),gen=genome,with.seq=F);print(i);return(x)},mc.cores=8L)
G4H_w25=unlist(do.call(c,G4hk))
seqlevels(G4H_w25, pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(G4H_w25) <- seqinf

save(G4H_w25,file="SC3_G4H_w25.RData")
lapply(c(1.2,1.5,2), function(i) {
	export(G4H_w25[G4H_w25$hl==i], con=paste0("SC3_G4H",i,".bed"))
})

#### Script to affect G4 to GenomicRanges (from bed file)

### first G4 calling at low threshold with the different windows
### to be done once and reuse for different bed

library("BSgenome.Hsapiens.UCSC.hg19")
genome <- BSgenome.Hsapiens.UCSC.hg19
source("G4HunterAccess_function.r")
seqinf <- seqinfo(genome)
seqlevels(seqinf) <- seqlevels(genome)[1:24]

w=100
newG4hk=mclapply((24:1),function(i) {x=mG4huntref(i,k=w,hl=0.5,gen=genome,with.seq=F);print(i);return(x)},mc.cores=8L)
mG4H_hg19_hl05_w100=do.call(c,newG4hk)
seqlevels(mG4H_hg19_hl05_w100, pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(mG4H_hg19_hl05_w100) <- seqinf
save(mG4H_hg19_hl05_w100,file='mG4H_hg19_hl05_w100.RData')

w=50
newG4hk=mclapply((24:1),function(i) {x=mG4huntref(i,k=w,hl=0.5,gen=genome,with.seq=F);print(i);return(x)},mc.cores=8L)
mG4H_hg19_hl05_w50=do.call(c,newG4hk)
seqlevels(mG4H_hg19_hl05_w50, force=T) <- seqlevels(seqinf)
seqinfo(mG4H_hg19_hl05_w50) <- seqinf
save(mG4H_hg19_hl05_w50,file='mG4H_hg19_hl05_w50.RData')

w=25
newG4hk=mclapply((24:1),function(i) {x=mG4huntref(i,k=w,hl=0.5,gen=genome,with.seq=F);print(i);return(x)},mc.cores=8L)
mG4H_hg19_hl05_w25=do.call(c,newG4hk)
seqlevels(mG4H_hg19_hl05_w25, force=T) <- seqlevels(seqinf)
seqinfo(mG4H_hg19_hl05_w25) <- seqinf
save(mG4H_hg19_hl05_w25,file='mG4H_hg19_hl05_w25.RData')

#### then I just load the proper RData to do the attribution
### example for Hg19
load("mG4H_hg19_hl05_w50.RData")
load("mG4H_hg19_hl05_w100.RData")
load("mG4H_hg19_hl05_w25.RData")

addG4H <- function(G4in, nan=NA)
{
	G4data <- mG4H_hg19_hl05_w25
	fov.test <- findOverlaps(G4in,G4data)
	G4in$G4H25 <- sapply(1:length(G4in), function(i)
	{
		if (i%in%queryHits(fov.test))
		{ x=max(abs(G4data[subjectHits(fov.test)[queryHits(fov.test)==i]]$max_score))}
		else
		{ x=nan}
		return(x)
	}
	)
	G4data <- mG4H_hg19_hl05_w50
	fov.test <- findOverlaps(G4in,G4data)
	G4in$G4H50 <- sapply(1:length(G4in), function(i)
	{
		if (i%in%queryHits(fov.test))
		{ x=max(abs(G4data[subjectHits(fov.test)[queryHits(fov.test)==i]]$max_score))}
		else
		{ x=nan}
		return(x)
	}
	)
	G4data <- mG4H_hg19_hl05_w100
	fov.test <- findOverlaps(G4in,G4data)
	G4in$G4H100 <- sapply(1:length(G4in), function(i)
	{
		if (i%in%queryHits(fov.test))
		{ x=max(abs(G4data[subjectHits(fov.test)[queryHits(fov.test)==i]]$max_score))}
		else
		{ x=nan}
		return(x)
	}
	)

	return(G4in)
}

### to use on a folder with bed files

### path for the input bed (created before)
path2 <- "bedFile_to_score/"
newbed <- dir(path2)
newbed2 <- tools::file_path_sans_ext(tools::file_path_sans_ext(newbed))
### path for the output bed (created before)
path2out <- "bedFile_scored/"

for (n in 1:length(newbed))
{
	G4in <- import(paste0(path2,newbed[n]),format="bed")
	seqlevels(G4in, pruning.mode="coarse") <- seqlevels(seqinf)
	seqinfo(G4in) <- seqinf
	G4in2 <- addG4H(G4in,nan=0)
	save(G4in2,file=paste0(path2out,newbed[n],".RData"))
	sc25 <- G4in2
	sc25$score <- sc25$G4H25
	export(sc25,con=paste0(path2out,newbed2[n],"_G4H25.bed"))
	sc50 <- G4in2
	sc50$score <- sc50$G4H50
	export(sc50,con=paste0(path2out,newbed2[n],"_G4H50.bed"))
	sc100 <- G4in2
	sc100$score <- sc100$G4H100
	sc100[is.na(sc100$score)]$score <- 0
	export(sc100,con=paste0(path2,newbed2[n],"_G4H100.bed"))
}

#### addG4H can be modified according to the Genome used as the G4data are hardcoded. The G4H search should also be performed again
