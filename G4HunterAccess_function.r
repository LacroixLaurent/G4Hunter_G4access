##### G4Hunter function for G4access
##### Function modified from the G4Hunter paper #########
##### Laurent Lqcroix 20210125 (laurent.lacroix@inserm.fr)

################################################################################
###### G4translate change the DNA code into G4Hunter code.
###### Only G or C are taken into account. non G/C bases are translated in 0
###### It is OK if N or U are in the sequence, but might be a problem if other letters or numbers are present
###### lowercase ARE not welcome
##########################################################
G4translate <- function(y,v1=1,v2=2,v3=3,v4=4)
	# x a DNAString or a DNAStringSet (or just a string of char)
{
	require(XVector)
	x= toupper(Rle(strsplit(as.character(y),NULL)[[1]]))
	xres=x
	runValue(xres)[runValue(x)=='C' & runLength(x)>3] <- -v4
	runValue(xres)[runValue(x)=='C' & runLength(x)==3] <- -v3
	runValue(xres)[runValue(x)=='C' & runLength(x)==2] <- -v2
	runValue(xres)[runValue(x)=='C' & runLength(x)==1] <- -v1
	runValue(xres)[runValue(x)=='G' & runLength(x)>3] <- v4
	runValue(xres)[runValue(x)=='G' & runLength(x)==3] <- v3
	runValue(xres)[runValue(x)=='G' & runLength(x)==2] <- v2
	runValue(xres)[runValue(x)=='G' & runLength(x)==1] <- v1
	runValue(xres)[runValue(x)!='C' & runValue(x)!='G'] <- 0	# N or U are not a problem
	Rle(as.numeric(xres))
}
################################################################################

################################################################################
##### return the G4Hscore. y is a string of char (DNA sequence)
G4Hscore <- function(y)
{
	y3=G4translate(y)
	mean(y3)
}
################################################################################

################################################################################
##### function to refine G4hunt results

G4startrun=function(y,chrom=chr,letter='C')	#y is a START
{
	if (letter(chrom,y)==letter)
	{
		while (letter(chrom,y)==letter & y!=1)
		{
			if (letter(chrom,y-1)==letter)
			{
				y <- y-1
			}else{
				break
			}
		}
	}else{
		y=y+1
		while (letter(chrom,y)!=letter) {y=y+1}
	}
	return(y)
}

G4endrun=function(y,chrom=chr,letter='C')	#y is a END
{
	if (letter(chrom,y)==letter)
	{
		while (letter(chrom,y)==letter & y!=length(chrom))
		{
			if (letter(chrom,y+1)==letter)
			{
				y <- y+1
			}else{
				break
			}
		}
	}else{
		y=y-1
		while (letter(chrom,y)!=letter) {y=y-1}
	}
	return(y)
}
################################################################################

#####################################################
## modified G4hunt to add the refining procedure and the max_score

mG4huntref <- function(i,k=25,hl=1.5,gen=genome,masked=5,with.seq=T,Gseq.only=T)
{
	require(GenomicRanges,Biostrings)
	#### k=RUNMEAN WINDOW SIZE, hl=threshold, i=chromosom number in the genome

	chr=gen[[i]]
	if (masked==2) {chr=injectHardMask(chr)}
	if (masked==3) {active(masks(chr))['RM']=T;chr=injectHardMask(chr)}
	if (masked==0) {active(masks(chr))=F;chr=injectHardMask(chr)}
	if (masked==4) {active(masks(chr))=T;chr=injectHardMask(chr)}

	tchr <- G4translate(chr)
	chr_G4hk <- runmean(tchr,k)
	if (class(gen)=="DNAStringSet")
	{seqname=names(gen)[i]
	}else{
		seqname=seqnames(gen)[i]
	}

	j <- hl
	chrCh <- Views(chr_G4hk, chr_G4hk<=(-j))
	chrGh <- Views(chr_G4hk, chr_G4hk>=j)

	IRC <- reduce(IRanges(start=start(chrCh),end=(end(chrCh)+k-1)))
	if (length(IRC)==0)
	{
		nxC <- GRanges()
	}else{
		nnIRC=IRC
		start(nnIRC)=sapply(start(IRC),G4startrun,letter='C',chrom=chr)
		end(nnIRC)=sapply(end(IRC),G4endrun,letter='C',chrom=chr)
		seqC=as.character(Views(chr,nnIRC))
		if (Gseq.only)
		{
			nnseqC=as.character(reverseComplement(Views(chr,nnIRC)))
		}else{
			nnseqC=as.character(Views(chr,nnIRC))
		}
		nG4scoreC=sapply(seqC,function(x) signif(G4Hscore(x),3))
		mscoreC <- signif(min(Views(chr_G4hk,IRC)),3)
		straC <- Rle(rep('-',length(IRC)))
		hlC <- Rle(rep(j,length(IRC)))
		kC <- Rle(rep(k,length(IRC)))
		maskC <- Rle(rep(masked,length(IRC)))
		nxC <- GRanges(seqnames=Rle(seqname),
									 ranges=nnIRC,
									 strand=straC,
									 score=nG4scoreC,
									 max_score=mscoreC,
									 hl=hlC,
									 k=kC,
									 mask=maskC,
									 sequence=nnseqC)
	}

	IRG <- reduce(IRanges(start=start(chrGh),end=(end(chrGh)+k-1)))
	if (length(IRG)==0)
	{
		nxG <- GRanges()
	}else{
		nnIRG=IRG
		start(nnIRG)=sapply(start(IRG),G4startrun,letter='G',chrom=chr)
		end(nnIRG)=sapply(end(IRG),G4endrun,letter='G',chrom=chr)
		nnseqG=as.character(Views(chr,nnIRG))
		nG4scoreG=sapply(nnseqG,function(x) signif(G4Hscore(x),3))
		mscoreG <- signif(max(Views(chr_G4hk,IRG)),3)
		straG <- Rle(rep('+',length(IRG)))
		hlG <- Rle(rep(j,length(IRG)))
		kG <- Rle(rep(k,length(IRG)))
		maskG <- Rle(rep(masked,length(IRG)))
		nxG <- GRanges(seqnames=Rle(seqname),
									 ranges=nnIRG,
									 strand=straG,
									 score=nG4scoreG,
									 max_score=mscoreG,
									 hl=hlG,
									 k=kG,
									 mask=maskG,
									 sequence=nnseqG)
	}

	nx <- sort(c(nxC,nxG),ignore.strand=T)
	names(nx) <- NULL
	if (with.seq==F) {nx$sequence=NULL}
	return(nx)
}


###############################################################
# function working directly with a list of threshold

mG4huntlistref <- function(i,k=25,hl=c(1,1.2,1.5,1.75,2),gen=genome,masked=5,with.seq=T,Gseq.only=T)
{
	require(GenomicRanges,Biostrings)
	#### k=RUNMEAN WINDOW SIZE, hl=threshold, i=chromosom number in the genome

	chr=gen[[i]]
	if (masked==2) {chr=injectHardMask(chr)}
	if (masked==3) {active(masks(chr))['RM']=T;chr=injectHardMask(chr)}
	if (masked==0) {active(masks(chr))=F;chr=injectHardMask(chr)}
	if (masked==4) {active(masks(chr))=T;chr=injectHardMask(chr)}

	tchr <- G4translate(chr)
	chr_G4hk <- runmean(tchr,k)
	if (class(gen)=="DNAStringSet")
	{seqname=names(gen)[i]
	}else{
		seqname=seqnames(gen)[i]
	}

	nx=list()

	hl=round(hl,2)		# there is a rounding error in the seq funtcion

	for (j in hl)
	{
		chrCh <- Views(chr_G4hk, chr_G4hk<=(-j))
		chrGh <- Views(chr_G4hk, chr_G4hk>=j)

		IRC <- reduce(IRanges(start=start(chrCh),end=(end(chrCh)+k-1)))
		if (length(IRC)==0)
		{
			nxC <- GRanges()
		}else{
			nnIRC=IRC
			start(nnIRC)=sapply(start(IRC),G4startrun,letter='C',chrom=chr)
			end(nnIRC)=sapply(end(IRC),G4endrun,letter='C',chrom=chr)
			seqC=as.character(Views(chr,nnIRC))
			if (Gseq.only)
			{
				nnseqC=as.character(reverseComplement(Views(chr,nnIRC)))
			}else{
				nnseqC=as.character(Views(chr,nnIRC))
			}
			nG4scoreC=sapply(seqC,function(x) signif(G4Hscore(x),3))
			mscoreC <- signif(min(Views(chr_G4hk,IRC)),3)
			straC <- Rle(rep('-',length(IRC)))
			hlC <- Rle(rep(j,length(IRC)))
			kC <- Rle(rep(k,length(IRC)))
			maskC <- Rle(rep(masked,length(IRC)))
			nxC <- GRanges(seqnames=Rle(seqname),
										 ranges=nnIRC,
										 strand=straC,
										 score=nG4scoreC,
										 max_score=mscoreC,
										 hl=hlC,
										 k=kC,
										 mask=maskC,
										 sequence=nnseqC)
		}

		IRG <- reduce(IRanges(start=start(chrGh),end=(end(chrGh)+k-1)))
		if (length(IRG)==0)
		{
			nxG <- GRanges()
		}else{
			nnIRG=IRG
			start(nnIRG)=sapply(start(IRG),G4startrun,letter='G',chrom=chr)
			end(nnIRG)=sapply(end(IRG),G4endrun,letter='G',chrom=chr)
			nnseqG=as.character(Views(chr,nnIRG))
			nG4scoreG=sapply(nnseqG,function(x) signif(G4Hscore(x),3))
			mscoreG <- signif(max(Views(chr_G4hk,IRG)),3)
			straG <- Rle(rep('+',length(IRG)))
			hlG <- Rle(rep(j,length(IRG)))
			kG <- Rle(rep(k,length(IRG)))
			maskG <- Rle(rep(masked,length(IRG)))
			nxG <- GRanges(seqnames=Rle(seqname),
										 ranges=nnIRG,
										 strand=straG,
										 score=nG4scoreG,
										 max_score=mscoreG,
										 hl=hlG,
										 k=kG,
										 mask=maskG,
										 sequence=nnseqG)
		}
		nx[[which(hl==j)]]=sort(c(nxC,nxG),ignore.strand=T)
		if (with.seq==F) {nx[[which(hl==j)]]$sequence=NULL}
		names(nx[[which(hl==j)]])=NULL
	}

	names(nx) <- as.character(hl)
	nx <- GRangesList(nx)
	return(nx)
}

