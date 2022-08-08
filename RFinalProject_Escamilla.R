#!/Usr/bin/env Rscript
#set working directory
setwd("/home/i2gds/dianae91/FinalProject_EscamillaD/00_rawdata")

#Step 2. Setup your experimental design.
rnames = c("Brain_1", "Brain_2", "Heart_1","Heart_2","Liver_1","Liver_2")
# create vectors of thre Meristem, Elongation and Mature stages
B2=rep("Brain",2)
H2=rep("Heart",2)
L2=rep("Liver",2)
# convert conditions to factors.
conditions = as.factor(c(B2,H2,L2))
# set order of the factor levels
conditions2 =factor(c(B2,H2,L2),
                    levels=c("Brain","Heart","Liver"))
# create "Design Matrix"
Samples = data.frame(row.names = rnames, 
                     condition = conditions)

# Step 3. Load data of raw data
B1 <- read.table('Brain_1.stats.txt',sep='\t',as.is=T, header=FALSE)
Brain_1<-B1[1,2]
B2 <- read.table('Brain_2.stats.txt',sep='\t',as.is=T, header=FALSE)
Brain_2<-B2[1,2]
H1 <- read.table('Heart_1.stats.txt',sep='\t',as.is=T, header=FALSE)
Heart_1<-H1[1,2]
H2 <- read.table('Heart_2.stats.txt',sep='\t',as.is=T, header=FALSE)
Heart_2<-H2[1,2]
L1 <- read.table('Liver_1.stats.txt',sep='\t',as.is=T, header=FALSE)
Liver_1<-L1[1,2]
L2 <- read.table('Liver_2.stats.txt',sep='\t',as.is=T, header=FALSE)
Liver_2<-L2[1,2]

InputDF<-data.frame(Brain_1=Brain_1,Brain_2=Brain_2,
                    Heart_1=Heart_1,Heart_2=Heart_2,
		    Liver_1=Liver_1,Liver_2=Liver_2)
#create the final dataframe(all other results with be output at this table)

Preprocessed_summary<-as.data.frame(t(InputDF))
colnames(Preprocessed_summary)<-c('Rawreads/replicate')

#change directories
setwd("/home/i2gds/dianae91/FinalProject_EscamillaD/01_processed")

#load the data of processed reads
B1 <- read.table('Brain_1.stats.txt',sep='\t',as.is=T, header=FALSE)
Brain_1<-B1[1,2]
B2 <- read.table('Brain_2.stats.txt',sep='\t',as.is=T, header=FALSE)
Brain_2<-B2[1,2]
H1 <- read.table('Heart_1.stats.txt',sep='\t',as.is=T, header=FALSE)
Heart_1<-H1[1,2]
H2 <- read.table('Heart_2.stats.txt',sep='\t',as.is=T, header=FALSE)
Heart_2<-H2[1,2]
L1 <- read.table('Liver_1.stats.txt',sep='\t',as.is=T, header=FALSE)
Liver_1<-L1[1,2]
L2 <- read.table('Liver_2.stats.txt',sep='\t',as.is=T, header=FALSE)
Liver_2<-L2[1,2]

InputDF1<-data.frame(Brain_1=Brain_1,Brain_2=Brain_2,
                    Heart_1=Heart_1,Heart_2=Heart_2,
                    Liver_1=Liver_1,Liver_2=Liver_2)

readsprocessed<-as.data.frame(t(InputDF1))
colnames(readsprocessed)<-c('Processedreads/replicate')

#add the number of processedreads/replicate to the Preprocessed_summary table
Preprocessed_summary$'Processedreads/replicate'<-readsprocessed$'Processedreads/replicate'

#Changing directories and loading data for mapped reads

setwd("/home/i2gds/dianae91/FinalProject_EscamillaD/02_mapping/Brain_1")
B1 <- read.table('Brain_1.stats.txt',sep='\t',as.is=T, header=FALSE)
Brain_1<-B1[,1]

setwd("/home/i2gds/dianae91/FinalProject_EscamillaD/02_mapping/Brain_2")
B2 <- read.table('Brain_2.stats.txt',sep='\t',as.is=T, header=FALSE)
Brain_2<-B2[,1]

setwd("/home/i2gds/dianae91/FinalProject_EscamillaD/02_mapping/Heart_1")
H1 <- read.table('Heart_1.stats.txt',sep='\t',as.is=T, header=FALSE)
Heart_1<-H1[,1]

setwd("/home/i2gds/dianae91/FinalProject_EscamillaD/02_mapping/Heart_2")
H2 <- read.table('Heart_2.stats.txt',sep='\t',as.is=T, header=FALSE)
Heart_2<-H2[,1]

setwd("/home/i2gds/dianae91/FinalProject_EscamillaD/02_mapping/Liver_1")
L1 <- read.table('Liver_1.stats.txt',sep='\t',as.is=T, header=FALSE)
Liver_1<-L1[,1]

setwd("/home/i2gds/dianae91/FinalProject_EscamillaD/02_mapping/Liver_2")
L2 <- read.table('Liver_2.stats.txt',sep='\t',as.is=T, header=FALSE)
Liver_2<-L2[,1]


InputDF2<-data.frame(Brain_1=Brain_1,Brain_2=Brain_2,
                    Heart_1=Heart_1,Heart_2=Heart_2,
                    Liver_1=Liver_1,Liver_2=Liver_2)

readsmapped<-as.data.frame(t(InputDF2))
colnames(readsmapped)<-c('Mappedreads/replicate')

#adding the number of mapped reads per replicate to the Preprocessed summary table
Preprocessed_summary$'Mappedreads/replicate'<-readsmapped$'Mappedreads/replicate'
write.table(Preprocessed_summary, "/home/i2gds/dianae91/FinalProject_EscamillaD/04_results/Preprocessed_summary.txt",row.names=TRUE,quote=FALSE)

##making graphs and summary tables of differential expressed genes 
##working directory
setwd("/home/i2gds/dianae91/FinalProject_EscamillaD/03_quantification/04_cuffdiff")
library("cummeRbund")
cuff<-readCufflinks()
cuff

##table of summary of differential expressed genes
BrainvsHeart<-getSig(cuff,"Brain","Heart",alpha=0.05,level="genes") 
TotalDEGenes<-c(length(BrainvsHeart))
f<-matrix(rnorm(9,0),1,3)
colnames(f)<-c('sample1','sample2','TotalDEGenes')
DE_summary<-data.frame(sample1=('Brain'),sample2=('Heart'),TotalDEGenes)

BrainvsLiver<-getSig(cuff,"Brain","Liver",alpha=0.05,level="genes")
TotalDEGenes<-c(length(BrainvsLiver))
f1<-matrix(rnorm(9,0),1,3)
colnames(f1)<-c('sample1','sample2','TotalDEGenes')
DF<-data.frame(sample1=('Brain'),sample2=('Liver'),TotalDEGenes)

##adding the new row to the DE_summary table
DE_summary<- rbind(DE_summary,DF)

HeartvsLiver<-getSig(cuff,"Heart","Liver",alpha=0.05,level="genes")
TotalDEGenes<-c(length(HeartvsLiver))
f2<-matrix(rnorm(9,0),1,3)
colnames(f2)<-c('sample1','sample2','TotalDEGenes')
DF1<-data.frame(sample1=('Heart'),sample2=('Liver'),TotalDEGenes)

##adding the new row to the DE_summary table
DE_summary<- rbind(DE_summary,DF1)

write.table(DE_summary, "/home/i2gds/dianae91/FinalProject_EscamillaD/04_results/DE_summary.txt",row.names=F,quote=FALSE)

##several graphs to see differential gene expression
##1 first graph-- density graph(1A)
dens<-csDensity(genes(cuff))
jpeg("/home/i2gds/dianae91/FinalProject_EscamillaD/04_results/figure1A.jpeg", width= 600, height= 400)
plot(dens)
dev.off()

##represent by replicates (1B)
dens_r<-csDensity(genes(cuff), replicates=T)
jpeg("/home/i2gds/dianae91/FinalProject_EscamillaD/04_results/figure1B.jpeg",width= 600, height= 400)
plot(dens_r)
dev.off()

##boxplot graph  (how homogenious are the values)(2A)
boxplot<-csBoxplot(genes(cuff))
jpeg("/home/i2gds/dianae91/FinalProject_EscamillaD/04_results/figure2A.jpeg",width= 600, height= 400)
plot(boxplot)
dev.off()

##boxplot graph (2B) with replicattes
boxplot1<-csBoxplot(genes(cuff),replicates=T)
jpeg("/home/i2gds/dianae91/FinalProject_EscamillaD/04_results/figure2B.jpeg",width= 600, height= 400)
plot(boxplot1)
dev.off()

##dendro plot graph (3A)
den<-csDendro(genes(cuff))
jpeg("/home/i2gds/dianae91/FinalProject_EscamillaD/04_results/figure3A.jpeg",width= 600, height= 400)
plot(den)
dev.off()

##dendro plot graph (3B) with reps
den1<-csDendro(genes(cuff),replicates=TRUE)
jpeg("/home/i2gds/dianae91/FinalProject_EscamillaD/04_results/figure3B.jpeg",width= 600, height= 400)
plot(den1)
dev.off()

##sig matrix 
mySigMat<-sigMatrix(cuff,level='genes',alpha=0.05)
jpeg("/home/i2gds/dianae91/FinalProject_EscamillaD/04_results/figure4.jpeg",width= 600, height= 400)
plot(mySigMat)
dev.off()

##DistHeat
myDistHeat<-csDistHeat(genes(cuff))
jpeg("/home/i2gds/dianae91/FinalProject_EscamillaD/04_results/figure5A.jpeg",width= 600, height= 400)
plot(myDistHeat)
dev.off()

##DistHeat with replicates
myDistHeat1<-csDistHeat(genes(cuff),replicates=TRUE)
jpeg("/home/i2gds/dianae91/FinalProject_EscamillaD/04_results/figure5B.jpeg",width= 600, height= 400)
plot(myDistHeat1)
dev.off()

##cluster plot
mySigGeneIds<-getSig(cuff,alpha=0.05,level='genes')
head(mySigGeneIds)
length(mySigGeneIds)
mySigGenes<-getGenes(cuff,mySigGeneIds)
mySigGenes
ic<-csCluster(mySigGenes,k=3)
head(ic$cluster)
icp<- csClusterPlot(ic)
jpeg("/home/i2gds/dianae91/FinalProject_EscamillaD/04_results/figure6.jpeg", width= 600, height=400)
plot(icp)
dev.off()

