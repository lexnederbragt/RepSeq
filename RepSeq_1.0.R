##################################################################################################################################################
# Estimation of Repeated Sequences from read Depths, Version 1.0 (Nov. 2009).
# Publication:
#	Nederbragt, A.J., Rounge, T.B., Kausrud, K and Jakobsen, K.S. 2009.
#	Identification and quantification of genomic repeats and sample contamination in assemblies of 454 pyrosequencing reads,
#	Sequencing
# Kyrre Linné Kausrud kyrreka@bio.uio.no
##################################################################################################################################################
# This program for R (http://www.r-project.org/) estimates the number of copies of each contig depending on the observed distributions of read depths for a sequenced genome.
# Being experimental, it comes with absolutely no warranty or liability whatsoever for accuracy, injury or damages.
# It generates some plots as it goes along, saving the results as a table (in .txt or .csv) and a figure (.pdf).

#Set parameters:

# Specify the path (in "quotes")of the folder containing the data file. The output figure and table will also be placed here. NB: Replace Windows single backslash (\)with double backslash (\\) in the file path.
	Folder.path<-"~/Documents/Coverage/Kyrre poisson/e_coli/"

# Specify (in "quotes") the name of the data file in the folder named above:
	Data.name<-"454AlignmentInfo.tsv"

# Remove "noise" contigs with read Depth below a certain minimum threshold. Recommended: 10
	min.lambda<-10

# Contigs shorter than this number will be excluded from the analysis:
	Length.minimum<-500

# Contigs with an estimated number of repeats larger than this number will be labeled in the output figure.
	show.numbers<-2	

# Either "csv" (comma separated values) or "txt" (tab separated values) for writing the results in a table of the appropriate format for your system.
	tabletype<-"txt"  	

#Choose between parametric and non-parametric test. (Parametric recommended.)
	Parametric<-F


## Optional adjustments:

#If Parametric=F, choose the quantiles of read Depth extracted from each contig's read Depth distributions. Determine the CI of the results.
# Recommended: 0.05 and 0.95
QL<-0.05
QU<-0.95

#If parametric=T, choose the number of standard errors of the mean to include in the final CI:
	N.SE<-3

#To adjust the CI for non-repeated read Depths, determine how far "down" the relative probability peak the CI is included. Smaller "Prob" values gives a larger CI. Recommended: 0.5
	Prob<-0.5

#############################################################################
##################################################################################
#Reads the innput data file, extracts the contig names and reformats the data so that they are usable for R:

setwd(Folder.path)
#A function for a "progress bar":
ticker<-function(x,x1,txt=NULL){plot(0:1,0:1,type="n",axes=F,xlab="Proportion completed",ylab="",main=txt);lines(c(0,x/x1),c(0.5,0.5),lwd=10);axis(1,at=c(0:10)/10);abline(v=c(0:10)/10,col="grey")}
ticker(1,1,txt="Reading Data file (be patient)")

#Reading raw data:
dat.gen<-read.delim2(Data.name,comment.char="#")
attach(dat.gen)
names(dat.gen)
#Find the lines where the raw innput file have put breaks for contig names:
Contig.names<-substring(dat.gen[which(is.finite(dat.gen[,3])==F),1],2)
brpts<-c(which(is.finite(dat.gen[,3])==F),length(dat.gen[,3])+1)
Contig.lengths<-brpts[2:length(brpts)]-brpts[1:(length(brpts)-1)]-1
Contig<-NULL;Contig.length<-NULL

#Make vectors giving contig number and length for each base pair:
for(i in 1:length(Contig.names)){
Contig<-append(Contig,rep(i,brpts[i+1]-brpts[i]-1))
Contig.length<-append(Contig.length,rep(Contig.lengths[i],brpts[i+1]-brpts[i]-1))
ticker(i,length(Contig.names),txt="Reading Contig lengths")
}
#Removing raw data to avoid memory overflow:
detach(dat.gen)
rm(dat.gen)

#Reading in data excluding the dummy rows containing contig names:
ticker(1,1,txt="Reading Data file (be patient -again)")
dat.gen<-read.delim2(Data.name,comment.char=">")
attach(dat.gen);dim(dat.gen);names(dat.gen)

##################################################################################################################################################################
#Main analysis:
#########################################################################################################################################################################

xmat1<-data.frame(aggregate(Depth[Contig.length>=Length.minimum],by=list(Contig[Contig.length>=Length.minimum]),mean,na.rm=T))	# A raw estimate of mean read Depth per contig
plot(density(xmat1[,2]),col="black",lwd=2,xlab="Read Depth, Lambda[i]",ylab="Frequency",main="");abline(v=min.lambda,lwd=2,col="red2")	#Plotting the chosen cutoff on the lambda distribution
contigs<-xmat1[,1][xmat1[,2]>=min.lambda]			#A vector of unique contig names, excluding contigs below min.lambda from further analysis.

#A first estimate of mean read Depths for these contigs
meanDepth <-xmat1[,2][xmat1[,2]>=min.lambda]

#A first analysis
	#Assuming that the data vector "Contig" contains the name of each contig, and the same length as
	# the numeric data vector "Depth" which gives the read Depth at each location.
#-----------------
xmat<-matrix(NA,length(contigs),5);qmat<-matrix(NA,length(contigs),200)
ql<-round(QL*200);qu<-round(QU*200)
for(i in 1:length(contigs)){								# For each contig name, do the following:
m1<-glm(Depth[Contig==contigs[i]]~1,family="quasipoisson")			# Find the maximum likelihood estimate for the mean (lambda) and the standard error of that estimate
xmat[i,1:2]<-exp(summary(m1)$coef[1:2])						# Load lambda and SE in matrix
xmat[i,3]<-median(Depth[Contig==contigs[i]])					# Load median read Depth
xmat[i,4:5]<-quantile(Depth[Contig==contigs[i]],prob=c(seq(0,1,by=0.01)[2:101][c(ql,qu)]))	# Empirical quantiles for the read Depth of each contig
qmat[i,]<-quantile(Depth[Contig==contigs[i]],prob=c(seq(0,1,by=0.005)[2:201]),na.rm=T)	# Empirical quantiles for the read Depth of each contig
plot(0,0,type="n",xlab="Remaining iterations",ylab="",axes=F);text(0,0,length(contigs)-i)	#Show progress 
}

#Thus, we have an empirical density of lambdas, mean or median depending on parametric/non-parametric test:
if(Parametric==T)dens<-density(c(xmat[,1]),from=min.lambda,n=5000)
if(Parametric==F)dens<-density(c(xmat[,3]),from=min.lambda,n=5000)

	#Assuming that the largest peak represent read Depth of single-copy regions, mean read Depth is estimated as: 
baseN<-dens$x[which(dens$y==max(dens$y))]	
	#Finding a range of almost-as-likely mean read Depths:
baseN.CI<-baseN		#Starting a vector of likely values
baseN.P<-max(dens$y)	#Starting a vector of relative probabilities
i<-which(dens$y==max(dens$y))-1
repeat{i<-i+1
if(dens$y[i+1]>=(Prob*max(dens$y))){baseN.CI<-append(baseN.CI,dens$x[i+1]);baseN.P<-append(baseN.P,dens$y[i+1])}
if(dens$y[i+1]<=(Prob*max(dens$y)))break}
i<-which(dens$y==max(dens$y))
repeat{i<-i-1
if(dens$y[i-1]>=(Prob*max(dens$y))){baseN.CI<-append(baseN.CI,dens$x[i-1]);baseN.P<-append(baseN.P,dens$y[i-1])}
if(dens$y[i-1]<=(Prob*max(dens$y)))break}

	#Plot of read Depths:
plot(dens,lwd=2,cex.lab=1.25,main=paste("Read Depth approx.",round(baseN,1),"[",round(min(baseN.CI),2),",",round(max(baseN.CI),2),"]"),cex.axis=1.2,las=1);abline(v=baseN,col="blue",lty=2,lwd=1)
lines(baseN.CI[order(baseN.CI)],baseN.P[order(baseN.CI)],lwd=3,col="red2");abline(v=min.lambda)

	#IF we want to minimize the distance from integer values, i.e. integer values are more likely than non-integers, we modify baseN within its CI:
x<-numeric(250)	#Data holder
interval<-c(seq(min(baseN.CI)/baseN,1,length=125),seq(1,max(baseN.CI)/baseN,length=125))	#Vector of multiplication factors
for(i in 1:250){								#For each factor:
baseN2<-baseN*interval[i]						#modify mean rad Depth
if(Parametric==T) x[i]<-sum((xmat[,1]/baseN2-round(xmat[,1]/baseN2))^2)	#find sum of squares with observation
if(Parametric==F) x[i]<-sum((xmat[,3]/baseN2-round(xmat[,3]/baseN2))^2)	#find sum of squares with observation
}
m<-sample(rep(interval[x==min(x)],2),size=1)							#and choose the minimum value m:
plot(interval,x,type="o",pch=16,xlab="Multiplicatior of baseN",ylab="Sum of squares from integer values",main=round(m,4));abline(v=m,col="grey")

###############################################################################
###############################################################################
#Processing results:
	#Extracting contig numbers from text (i.e. names):
True_numbers<-F	#If this is "T", then the characters after "contig" is extracted for use as names. This may look strange when many numbers are missing, then use "F", as the true numbers in any case appear in the plot and table.
if(True_numbers==T)contig.numbers<-substring(contigs,9,11)
if(True_numbers==F)contig.numbers<-1:length(xmat[,1])

	#Plotting results and printing the contig number of those repeated more than show.number times:

baseN2<-baseN*m			#If distance from integers above has been used, multiply baseN with this
if(Parametric==T){							#If parametric test:
Ncopies<-xmat[,1]/baseN2						#Estimated number of copies of each contig is thus mean read Depth/singlecopy read Depth		
Ncopies.UCL<-(xmat[,1]+N.SE*xmat[,2])/min(baseN.CI)		#And LCL/UCLs are found from an appropriate number of standard errors of the mean lambda and the CI of the mean read Depth
Ncopies.LCL<-(xmat[,1]-N.SE*xmat[,2])/max(baseN.CI)
}
if(Parametric==F){							#If non-parametric test:
Ncopies<-xmat[,3]/baseN2						#Estimated number of copies of each contig is thus median read Depth/singlecopy read Depth	
Ncopies.UCL<-(qmat[,qu])/min(baseN.CI)				#And LCL/UCLs are found from the quantiles of each contig read Depth and the CI of the mean read Depth
Ncopies.LCL<-(qmat[,ql])/max(baseN.CI)
}

# Copy number estimate versus read Depth graph
plot(Ncopies,meanDepth,pch=16,las=1,xlab="Number of copies per contig",ylab="Per-contig read Depth",cex.lab=1.25,cex.axis=1.2,xlim=c(-0.1,max(Ncopies.UCL)+0.1))
axis(1,at=c(1:round(max(Ncopies.UCL))),cex.axis=1.2);abline(v=c(1:round(max(Ncopies.UCL+1))),lty=2,col="grey")
segments(Ncopies.UCL,meanDepth,Ncopies.LCL,meanDepth,lwd=2)

###############################################################################################################################
	#Plotting the results; number of copies for each contig, with a confidence interval shown and contig number for all contigs probably existing in two or more copies.

pdf(file="RepSeq_Figure.pdf")
plot(Ncopies,as.numeric(contig.numbers),pch=16,las=1,xlab="Number of copies per contig",ylab="Contig nr.",cex.lab=1.25,cex.axis=1.2,xlim=c(-0.1,max(Ncopies.UCL)+0.5))
axis(1,at=c(1:round(max(Ncopies.UCL))),cex.axis=1.2);abline(v=c(1:round(max(Ncopies.UCL+1))),lty=2,col="grey")
segments(Ncopies.UCL,as.numeric(contig.numbers),Ncopies.LCL,as.numeric(contig.numbers),lwd=2)
text(Ncopies.UCL[Ncopies>=show.numbers],as.numeric(contig.numbers)[Ncopies>=show.numbers]+0.4,substring(Contig.names[contigs],7)[Ncopies>=show.numbers],font=2,adj=0,cex=0.75,col="red2")
points(rep(0,max(as.numeric(contig.numbers)))[-as.numeric(contig.numbers)],c(1:max(as.numeric(contig.numbers)))[-as.numeric(contig.numbers)],col="grey",cex=3,pch=".")
dev.off()

X11()
plot(Ncopies,as.numeric(contig.numbers),pch=16,las=1,xlab="Number of copies per contig",ylab="Contig nr.",cex.lab=1.25,cex.axis=1.2,xlim=c(-0.1,max(Ncopies.UCL)+0.5))
axis(1,at=c(1:round(max(Ncopies.UCL))),cex.axis=1.2);abline(v=c(1:round(max(Ncopies.UCL+1))),lty=2,col="grey")
segments(Ncopies.UCL,as.numeric(contig.numbers),Ncopies.LCL,as.numeric(contig.numbers),lwd=2)
text(Ncopies.UCL[Ncopies>=show.numbers],as.numeric(contig.numbers)[Ncopies>=show.numbers]+0.4,substring(Contig.names[contigs],7)[Ncopies>=show.numbers],font=2,adj=0,cex=0.75,col="red2")
points(rep(0,max(as.numeric(contig.numbers)))[-as.numeric(contig.numbers)],c(1:max(as.numeric(contig.numbers)))[-as.numeric(contig.numbers)],col="grey",cex=3,pch=".")

#Make the results into a data frame to be written as a .csv or .txt file to a folder where it can be opened in excel or similar for visual editing.
results_table<-as.data.frame(cbind(as.character(Contig.names[contigs]),meanDepth,Ncopies,Ncopies.LCL,Ncopies.UCL));names(results_table)<-c("Contig","Mean read depth","Estimated number of copies","Estimated number of copies (LCL)","Estimated number of copies (UCL)")

if(tabletype=="txt")write.table(results_table,file="RepSeq.Table.txt", sep="\t", quote=F, row.names=F)
if(tabletype=="csv")write.csv2(results_table,file="RepSeq.Table.csv", quote=F, row.names=F)
#################################

