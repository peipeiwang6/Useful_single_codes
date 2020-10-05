args=commandArgs(TRUE)
your_data = args[1]
res <- read.table(your_data,head=T,sep='\t')
res <- res[res[,4]!=0,]
Enrichment <- function(k,n,C,G){
	return((k / C) / (n/ G))
	}
enrichment <- c()
for(i in 1:nrow(res)){
	numbers <- matrix(as.numeric(res[i,2:5]),nrow = 2)
	p <- fisher.test(numbers, alternative = "two.sided")[[1]][1]
	a = as.numeric(res[i,2])
	b = as.numeric(res[i,3])
	cc = as.numeric(res[i,4])
	d = as.numeric(res[i,5])
	if(Enrichment(a, a+b, a+cc, a+b+cc+d) >= 1) direction = '+' else direction = '-'
	enrichment <- rbind(enrichment, c(res[i,],direction,p))
	}
write.table(enrichment,paste(your_data,'.fisher.pvalue',sep=''),row.names=F,col.names=F,quote=F,sep='\t')

dat <- read.table(paste(your_data,'.fisher.pvalue',sep=''),head=F,sep='\t',stringsAsFactors=F)
dat <- cbind(dat,p.adjust(dat[,7], method = "BH"))
dat <- dat[order(dat[,8]),]
dat <- cbind(dat,'BP'='')
dat <- cbind(dat,'CC'='')
dat <- cbind(dat,'MF'='')
write.table(dat,paste(your_data,'.fisher.qvalue',sep=''),row.names=F,col.names=F,quote=F,sep='\t')

library('GOstats')
library('GSEABase')
dat <- read.table(paste(your_data,'.fisher.qvalue',sep=''),head=F,sep='\t',stringsAsFactors=F)
for(i in 1:nrow(dat)){
	tryCatch( {
		if(!is.null(getGOTerm(dat[i,1])$BP[1])) dat[i,9] <- getGOTerm(dat[i,1])$BP[1]
		if(!is.null(getGOTerm(dat[i,1])$CC[1])) dat[i,10] <- getGOTerm(dat[i,1])$CC[1]
		if(!is.null(getGOTerm(dat[i,1])$MF[1])) dat[i,11] <- getGOTerm(dat[i,1])$MF[1]
		},
		error = function(e) {print(paste("no GO for ",dat[i,1]));NaN},
		finally = {})
	}
write.table(dat,paste(your_data,'.fisher.qvalue_GO_term',sep=''),row.names=F,col.names=F,quote=F,sep='\t')
