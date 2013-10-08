args <- commandArgs(TRUE);
frst <- as.integer(args[1]);
last <- as.integer(args[2]);

#######################################################################################################
# load functions
source('readms.output.msms.sel.R', chdir = TRUE); # Version of Dick Hudson's script from msdir modified for msms.
source('calc_wcFstats.R', chdir = TRUE); # From Eva Chan.

#########################################################################
# parameters
u=10^-8 
Nsub = 25000 #Ne for each subpop
len <- 2000;
theta = 4*Nsub*u*len;
rho = 4*Nsub*u*10*len;
splits <- c(500,2000,5000,10000,15000,20000,25000,30000,40000,50000,60000,80000,100000);
ns <- c(rep(3000,length(splits)));
nsam <- 100;
subpop <- c(rep(1,(nsam/4)),rep(2,(nsam/4)));
idx <- 0;
res.max <- matrix(ncol=24,nrow=(1000*length(splits)));
cut <- 5;

cat('Updated at 12.28 on 5.3 \n');
#########################################################################
for(t in frst:last){
	tsplit = splits[t]/(4*Nsub);#subs split this many gens before admixture
	tsel = tsplit*0.99;
	
	#################################
	## simulate set of selected loci, 1 per iteration
	system(paste('/progs/msms/bin/msms -N ',Nsub,' -ms 100 1000 -t ',theta,' -r ',rho,' ',as.integer(len),' -I 2 50 50 -n 1 1 -n 2 1 -ej ',tsplit,' 2 1 -SI ',tsel,' 2 ',(1/Nsub),' 0 -Sc 0 1 ',ns[t],' ',(ns[t]-1000),' 0 -Sp 0.5 -SFC -Smark > Admix.sel.N25K.',splits[t],'.out',sep=''));
	sel <- read.ms.output(file.ms.output=paste('Admix.sel.N25k.',splits[t],'.out',sep=''));
	selset <- which(sel[[1]]>1);
	nsam <- sel[[6]];
	
	#################################
	# iterate on this split grid point.  Each iteration is a new exome
	for(i in 1:100){
		if(i%%5==0){
			cat('Now searching splits ',splits[t],' on iteration ',i,' \n');
		}
		idx <- idx+1;
		fst.temp <- fst.temp2 <- cls <- clse <- chi <- che <- NULL;
		res.max[idx,1] <- splits[t];
		res.max[idx,2] <- nsam/2;
		
		## start with stats for selected locus	
		segs <- sel[[1]][[selset[i]]];
		posits <- sel[[5]][[selset[i]]];	
		haps <- sel[[2]][[selset[i]]];
		geno <- matrix(ncol=nsam/2,nrow=segs)
		for(g in 1:(nsam/2)){
			geno[,g] <- haps[(g*2)-1,]+haps[(g*2),];
		}
		x <- calc_wcFstats(geno,subpop);
		res.max[idx,3] <- x[[1]][which(posits==0.5),5];
		res.max[idx,4] <- x[[2]][2];
		fst.temp <- c(fst.temp,x[[2]][2]); # fst for whole exon
		fst.temp2 <- c(fst.temp2,x[[1]][,5]); # fst by snp
		snpcls <- rep(2,segs);
		snpcls[which(sel[[5]][[i]]==0.5)] <- 3;
		cls <- c(cls,snpcls);
		sel.snp <- which(cls==3);  	

		##############################
		## simulate neutral loci
		system(paste('/progs/msms/bin/msms -N ',Nsub,' -ms 100 550 -t ',theta,' -r ',rho,' ',as.integer(len),' -I 2 50 50 -n 1 1 -n 2 1 -ej ',tsplit,' 2 1 > Admix.neut.simplesplit.N25k.',splits[t],'.',i,'.out',sep=''));
		chroms <- sample(c(1:1000),25);
		for(c in chroms){
			out <- read.ms.output(file.ms.output=paste('Admix.neut.simplesplit.N25k.',splits[t],'.',i,'.out',sep=''));
			set2 <- which(out[[1]]>1);
		
			#############################
		
			for(f in set2[1:length(set2)]){
				segs <- out[[1]][f];
				haps <- out[[2]][[f]];
				geno <- matrix(ncol=nsam/2,nrow=segs)
				for(g in 1:(nsam/2)){
					geno[,g] <- haps[(g*2)-1,]+haps[(g*2),];
				}
			
				x <- calc_wcFstats(geno,subpop);
				fst.temp <- c(fst.temp,x[[2]][2]);
				fst.temp2 <- c(fst.temp2,x[[1]][,5]);
				cls <- c(cls,rep(1,segs));	
			}
		}
		tots1 <- rowSums(geno[,1:(nsam/4)]);tots2 <- rowSums(geno[,((nsam/4)+1):(nsam/2)]);
		
		#############################
		## collect results
		# Fst by whole exon ----------------
		res.max[idx,5] <- length(set2);	
		res.max[idx,6] <- mean(fst.temp);    ## Fst for whole exon
		res.max[idx,7] <- median(fst.temp);
		res.max[idx,8] <- max(fst.temp);
		res.max[idx,9] <- length(which(fst.temp>=res.max[idx,4]));
		res.max[idx,10] <- sort(fst.temp,decreasing=TRUE)[round(length(fst.temp)*0.05)];
		res.max[idx,11] <- sort(fst.temp,decreasing=TRUE)[round(length(fst.temp)*0.01)];
		res.max[idx,12] <- sort(fst.temp,decreasing=TRUE)[2];
		# Fst by SNP ------------------------
		res.max[idx,13] <- length(which(fst.temp2>=fst.temp2[which(cls==3)]));   
		res.max[idx,14] <- length(which(is.element(fst.temp2[which(cls==2)],sort(fst.temp2[which(cls!=3)],decreasing=TRUE)[1:10])==TRUE));
		#res.max[idx,15] <- sort(fst.temp2,decreasing=TRUE)[2];
		# Chi-sq method ---------------------
		res.max[idx,18] <- mean(fst.temp2);
		res.max[idx,20] <- length(which(tots1==(nsam/4) & tots2==0))+length(which(tots2==(nsam/4) & tots1==0));
		res.max[idx,21] <- var(fst.temp2);
		res.max[idx,23] <- length(which(fst.temp2[which(cls==2)]==fst.temp2[which(cls==3)]));	
	}			
	#colnames(res.max) <- c('SplitGens','ChromosPop','SelSite.Fst','Num.Neut.Exons','MeanFst','MedFst','MaxFst','SelSiteRank','Fst.5Percent','Fst.1Percent','Fst.Num2');
	write.table(res.max,file=paste('Fst.grid.results.GEN.',frst,'.txt',sep=''),row.names=FALSE,quote=FALSE);	
}

###########################################################################################################################
