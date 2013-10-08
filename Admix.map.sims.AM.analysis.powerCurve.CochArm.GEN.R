args <- commandArgs(TRUE);
mix <- as.integer(args[1]);
strt <- as.integer(args[2]);
end <- as.integer(args[3]);

timestrt=proc.time()

cat('Updated at 228 on 5/10 \n');
#########################################################################################################

splits <- c(500,5000,10000,20000,30000,50000,60000,80000,100000,150000,450000);
mixers <- c(2,4,6,8,seq(from=10,to=50,by=10));
nsams <- c(100);
idx <- 0;
iters <- 100;
dir= '/Admix.sims.new/AM/admixed.pops'
dir2='/Admix.sims.new'

g=mix;
#########################################################################################################
res.max <- matrix(ncol=25,nrow=(100*length(splits)*length(nsams)));
for(t in end:strt){
	selsite <- read.table(file=paste(dir,'Admixed.pop.n500.',splits[t],'.selsites',sep=''),header=F)[,1];
		p1 <- p2 <- NULL;
		for(i in 1:200){ 																									
			phens <- read.table(file=paste(dir,'Admixed.pop.n500.',splits[t],'.',g,'.',i,'.phenos',sep=''),header=TRUE);
			p1[i] <- length(which(phens[,1]==1))
			p2[i] <- length(which(phens[,1]==2))
		}

	hits <- which(p1>=50 & p2>=50);
	for(i in hits[1:100]){
		# Pull in 
		phens <- read.table(file=paste(dir,'Admixed.pop.n500.',splits[t],'.',g,'.',i,'.phenos',sep=''),header=TRUE);
		temp <- read.table(file=paste(dir,'Admixed.pop.n500.',splits[t],'.',g,'.',i,'.pop',sep=''),header=F);
		segs <- read.table(file=paste(dir,'Admixed.pop.n500.',splits[t],'.',g,'.',i,'.sites',sep=''),header=F)[,1];
		breaks <- which(segs==0);	
		segs.cln <- segs[-breaks];
				
		exon <- NULL;breaks2 <- c(0,breaks,(length(segs)+1));cls.big <- NULL;
		for(e in 1:(length(breaks2)-1)){
			if(e==250){
				cls.big <- c(cls.big,rep(2,(breaks2[e+1]-breaks2[e]-1)));
			}else{
				cls.big <- c(cls.big,rep(1,(breaks2[e+1]-breaks2[e]-1)));
			}
			exon <- c(exon,rep(e,(breaks2[e+1]-breaks2[e]-1)))
		}
		cls.big[which(segs.cln==selsite[i])] <- 3;
							
		idx <- idx+1;
		nsam <- nsams[1];
		cut <- (0.2*nsam);
		cat('Now on sam mixers ',g,' for iteration ',which(hits==i),' and splitgen ',splits[t],' \n');
					
		# Sub-sample down to desired sample size and assemble hap list
		set <- c(sample(which(phens[,1]==1),nsam/2),sample(which(phens[,1]==2),nsam/2));
		set2 <- NULL;
		for(d in 1:length(set)){
			set2 <- c(set2,((set[d]*2)-1),((set[d]*2)-0));
		}
		mends <- phens[set,1];
		haps.temp <- temp[set2,];
					
		# Convert coding of chromosome SNPs back to 0/1 in proper matrix structure
		haps <- matrix(0,ncol=ncol(haps.temp),nrow=nrow(haps.temp));
		for(d in 1:nrow(haps.temp)){
			haps[d,which(is.element(haps.temp[d,],c(1,4))==TRUE)] <- 1;												
		}
					
		# Convert coding of chromosome SNPs back to 0/1 in proper matrix structure
		geno <- matrix(0,ncol=(nsam),nrow=ncol(haps));
		for(d in 1:(nsam)){
			geno[,d] <- haps[(d*2)-1,]+haps[(d*2),];						
		}
				 	
		## Run association tests
		arm <- NULL;tots <- rowSums(geno);idx3 <-0;
		keeps=which(tots>cut & tots<((nsam*2)-cut));
		for(d in which(tots>cut & tots<((nsam*2)-cut))){
			#idx3 <-idx3+1;
					
			# Coch Arm test on Mend-Dom trait
			data <- data.frame(
				ph=factor(mends),
				genotype=ordered(geno[d,],
				levels=c('0','1','2')))						
			table.1 <- table(data);
			x <- prop.trend.test(
				x=table.1['2',],
				n=margin.table(table.1,margin=2),
				score=c(0,1,1))
			if(length(which(x[[3]][1]>0))==1){
				#idx3<-idx3+1;
				arm[d] <- c(x[[3]][1]);
			}else{
				keeps=keeps[-(which(keeps==d))];
			}	
		}
					
		## Adjust indexing to accomodate SNP frequency selection
		exon2 <- exon[keeps];
		cls <- cls.big[keeps];
		arm=arm[keeps];	
					
		sel.snp <- which(cls==3);
		sel.linked <- which(cls==2);
							
		## Assemble exon-wise p-values
		arm.p<- NULL;
		for(e in 1:(length(breaks)+1)){
			if(length(which(exon2==e))>0){											
				chi.temp4 <- arm[which(exon2==e)];				
				arm.p[e] <- pchisq((-2*sum(log(chi.temp4),na.rm=TRUE)),(2*length(which(chi.temp4>0))),ncp=0,lower.tail=FALSE,log.p=FALSE); 
			}
		}
		
		## Store results in result matrix
		res.max[idx,1] <- splits[t];
		res.max[idx,2] <- mix;
		res.max[idx,3] <- nsam;
		
		# By SNP
		res.max[idx,13] <- length(which(arm[sel.linked]==arm[sel.snp]));
		res.max[idx,14] <- length(which(arm<=arm[sel.snp]));
		res.max[idx,15] <- arm[sel.snp];
		res.max[idx,16] <- mean(log(arm),na.rm=TRUE);
				
		# By exon					
		res.max[idx,19] <- length(which(arm.p<=arm.p[250])); # by exon with mend trait and arm test
		res.max[idx,22] <- arm.p[250];										
		
		write.table(res.max[1:idx,],file=paste(dir2,'GAtest.SamSize.n',nsam,'.g',mix,'.',strt,'.txt',sep=''),row.names=F,quote=F);

			
		##############################################################	
		chroms <- sample(c(1:200)[-i],25);arm.big <- arm;
		for(c in 1:25){
			if(c%%5==0){
				cat('Now on sam mixers ',g,' for iteration ',which(hits==i),' and splitgen ',splits[t],' and chrom ',c,' \n');
			}	
			temp <- read.table(file=paste(dir,'Admixed.pop.n500.',splits[t],'.',g,'.',chroms[c],'.pop',sep=''),header=F);
			segs <- read.table(file=paste(dir,'Admixed.pop.n500.',splits[t],'.',g,'.',chroms[c],'.sites',sep=''),header=F)[,1];
			breaks <- which(segs==0);	
			segs.cln <- segs[-breaks];
			
			exon <- NULL;breaks2 <- c(0,breaks,(length(segs)+1));cls.big <- NULL;
			for(e in 1:(length(breaks2)-1)){
				if(e==250){
					cls.big <- c(cls.big,rep(2,(breaks2[e+1]-breaks2[e]-1)));
				}else{
					cls.big <- c(cls.big,rep(1,(breaks2[e+1]-breaks2[e]-1)));
				}
				exon <- c(exon,rep(e,(breaks2[e+1]-breaks2[e]-1)))
			}
						
			# Sub-sample down to desired sample size and assemble hap list
			haps.temp <- temp[set2,];
					
			# Convert coding of chromosome SNPs back to 0/1 in proper matrix structure
			haps <- matrix(0,ncol=ncol(haps.temp),nrow=nrow(haps.temp));
			for(d in 1:nrow(haps.temp)){
				haps[d,which(is.element(haps.temp[d,],c(1,4))==TRUE)] <- 1;												
			}
				
			# Convert coding of chromosome SNPs back to 0/1 in proper matrix structure
			geno <- matrix(0,ncol=(nsam),nrow=ncol(haps));
			for(d in 1:(nsam)){
				geno[,d] <- haps[(d*2)-1,]+haps[(d*2),];						
			}
			
			geno <- geno[which(cls.big==1),]
			cls.big=cls.big[which(cls.big==1)]
			
			## Run association tests
			arm <- NULL;tots <- rowSums(geno);idx3 <-0;
			keeps=which(tots>cut & tots<((nsam*2)-cut));
			for(d in which(tots>cut & tots<((nsam*2)-cut))){
				#idx3 <-idx3+1;
						
				# Coch Arm test on Mend-Dom trait
				data <- data.frame(
					ph=factor(mends),
					genotype=ordered(geno[d,],
					levels=c('0','1','2')))						
				table.1 <- table(data);
				x <- prop.trend.test(
					x=table.1['2',],
					n=margin.table(table.1,margin=2),
					score=c(0,1,1))
				if(length(which(x[[3]][1]>0))==1){
					#idx3<-idx3+1;
					arm[d] <- c(x[[3]][1]);
				}else{
					keeps=keeps[-(which(keeps==d))];
				}	
			}
			arm.big <- c(arm.big,arm[keeps]);		
			
			## Adjust indexing to accomodate SNP frequency selection
			exon2 <- exon[keeps];
			cls <- cls.big[keeps];
			arm=arm[keeps];					
					
			## Assemble exon-wise p-values
			for(e in c(1:(length(breaks)+1))){
				if(length(which(exon2==e))>0){											
					chi.temp4 <- arm[which(exon2==e)];
					arm.p <- c(arm.p,pchisq((-2*sum(log(chi.temp4),na.rm=TRUE)),(2*length(which(chi.temp4>0))),ncp=0,lower.tail=FALSE,log.p=FALSE)); 
				}
			}
			cat('After chrom ',c,' the rank is ',length(which(arm.p<=arm.p[250])),'\n');
			
		}# end of neut chroms loop
				
		## Store results in result matrix
		res.max[idx,1] <- splits[t];
		res.max[idx,2] <- mix;
		res.max[idx,3] <- nsam;
		
		# By SNP
		res.max[idx,13] <- length(which(arm.big[sel.linked]==arm.big[sel.snp]));
		res.max[idx,14] <- length(which(arm.big<=arm.big[sel.snp]));
		res.max[idx,15] <- arm.big[sel.snp];
		res.max[idx,16] <- mean(log(arm.big),na.rm=TRUE);
				
		# By exon					
		res.max[idx,19] <- length(which(arm.p<=arm.p[250])); # by exon with mend trait and arm test
		res.max[idx,22] <- arm.p[250];										
		
		write.table(res.max[1:idx,],file=paste(dir2,'GAtest.SamSize.GEN.n',nsam,'.g',mix,'.',strt,'.txt',sep=''),row.names=F,quote=F);
		
	}# end of iters loop
	
}# end of splits loop
timend<-proc.time();
cat(timestrt,'\n');
cat(timend,'\n');