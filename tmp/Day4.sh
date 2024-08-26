#Sample size
grep -F 'EUR' all_phase3.king.psam | wc -l

#Number of genetic variants
wc chr1-22.bim -l

./plink --bfile chr1-22 --snps-only --write-snplist


#Allele frequencies and missingness}
./plink --bfile chr1-22 --snps-only --freq --within pop_info.pheno


grep -F 'AFR' plink.frq.strat | wc -l


grep -F 'AFR' plink.frq.strat > freq_report.afr
grep -F 'AMR' plink.frq.strat > freq_report.amr
grep -F 'EUR' plink.frq.strat > freq_report.eur
grep -F 'EAS' plink.frq.strat > freq_report.eas
grep -F 'SAS' plink.frq.strat > freq_report.sas
grep -F 'AFR' plink.frq.strat | awk '$6 >0' freq_report.afr | wc -l
grep -F 'EUR' plink.frq.strat | awk '$6 >0' freq_report.eur | wc -l
grep -F 'EAS' plink.frq.strat | awk '$6 >0' freq_report.eas | wc -l
grep -F 'AMR' plink.frq.strat | awk '$6 >0' freq_report.amr | wc -l
grep -F 'SAS' plink.frq.strat | awk '$6 >0' freq_report.sas | wc -l


./plink --bfile chr1-22 --missing --within pop_info.pheno


awk '$4 > 0' plink.lmiss | wc -l


#Cross-population allele frequency comparisons
R
library(dplyr)
library(ggplot2)
freq <-read.table("plink.frq.strat", header =T)
plotDat <- freq %>%
  mutate(AlleleFrequency = cut(MAF, seq(0, 1, 0.25))) %>%
  group_by(AlleleFrequency, CLST) %>%
  summarise(FractionOfSNPs = n()/nrow(freq) * 100)

ggplot(na.omit(plotDat),
       aes(AlleleFrequency, FractionOfSNPs, group = CLST, col = CLST)) +
  geom_line() +
  scale_y_continuous(limits = c(0, 12)) +
  ggtitle("Distribution of allele frequency across genome")


#Calculation of Fst
#DO NOT PROCESS THIS NEXT COMMAND - DOES NOT WORK ON WORKSHOP COMPUTERS
#./plink2 --bfile chr1-22 --fst POP --pheno pop_info.pheno


#RESUME FROM HERE ONWARDS
#Derive information on pairwise R2 between all SNPs
./plink --bfile chr1-22 \
--keep-cluster-names AFR \
--within pop_info.pheno \
--r2 \
--ld-window-r2 0 \
--ld-window 999999 \
--ld-window-kb 2500 \
--threads 30 \
--out chr1-22.AFR

./plink --bfile chr1-22 \
--keep-cluster-names EUR \
--within pop_info.pheno \
--r2 \
--ld-window-r2 0 \
--ld-window 999999 \
--ld-window-kb 2500 \
--threads 30 \
--out chr1-22.EUR

cat chr1-22.AFR.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > chr1-22.AFR.ld.summary
cat chr1-22.EUR.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > chr1-22.EUR.ld.summary



R
install.packages("dplyr")
install.packages("stringr")
install.packages("ggplot2")
library(dplyr)
library(stringr)
library(ggplot2)


dfr<-read.delim("chr1-22.AFR.ld.summary",sep="",header=F,check.names=F, stringsAsFactors=F)
colnames(dfr)<-c("dist","rsq")
dfr$distc<-cut(dfr$dist,breaks=seq(from=min(dfr$dist)-1,to=max(dfr$dist)+1,by=100000))
dfr1<-dfr %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))


ggplot()+
  geom_point(data=dfr1,aes(x=start,y=mean),size=0.4,colour="grey20")+
  geom_line(data=dfr1,aes(x=start,y=mean),size=0.3,alpha=0.5,colour="grey40")+
  labs(x="Distance (Megabases)",y=expression(LD~(r^{2})))+
  scale_x_continuous(breaks=c(0,2*10^6,4*10^6,6*10^6,8*10^6),labels=c("0","2","4","6","8"))+
  theme_bw()


#Distribution of LD-block length

./plink --bfile chr1-22 --keep-cluster-names AFR --blocks no-pheno-req no-small-max-span --blocks-max-kb 250 --within pop_info.pheno  --threads 30 --out AFR
./plink --bfile chr1-22 --keep-cluster-names EUR --blocks no-pheno-req no-small-max-span --blocks-max-kb 250 --within pop_info.pheno  --threads 30 --out EUR

R
dfr.afr <- read.delim("AFR.blocks.det",sep="",header=T,check.names=F,stringsAsFactors=F)
colnames(dfr.afr) <- tolower(colnames(dfr.afr))

dfr.eur <- read.delim("EUR.blocks.det",sep="",header=T,check.names=F,stringsAsFactors=F)
colnames(dfr.eur) <- tolower(colnames(dfr.eur))

dfr.amr <- read.delim("AMR.blocks.det",sep="",header=T,check.names=F,stringsAsFactors=F)
colnames(dfr.amr) <- tolower(colnames(dfr.amr))

dfr.sas <- read.delim("SAS.blocks.det",sep="",header=T,check.names=F,stringsAsFactors=F)
colnames(dfr.sas) <- tolower(colnames(dfr.sas))

dfr.eas <- read.delim("EAS.blocks.det",sep="",header=T,check.names=F,stringsAsFactors=F)
colnames(dfr.eas) <- tolower(colnames(dfr.eas))

#Code for  one line (dfr.afr) is supplied. Others must be generated in order to complete the plot
plot (density(dfr.afr$kb), main="LD block length distribution", ylab="Density",xlab="LD block length (Kb)" )
lines (density(dfr.eur$kb), col="blue")
lines (density(dfr.eas$kb), col="red")
lines (density(dfr.amr$kb), col="purple")
lines (density(dfr.sas$kb), col="green")
legend("topright",c("AFR","EAS","EUR","SAS","AMR"), fill=c("black","red","blue","green","purple"))



#Principle Component Analysis

./plink --bfile chr1-22 --indep-pairwise 250 25 0.1 --maf 0.1 --threads 30 --out chr1-22.ldpruned_all_1kgv2
./plink --bfile chr1-22 --extract chr1-22.ldpruned_all_1kgv2.prune.in  --pca --threads 30


#R Plot
R
require('RColorBrewer')
options(scipen=100, digits=3)

eigenvec <- read.table('plink.eigenvec', header = F, skip=0, sep = ' ')
rownames(eigenvec) <- eigenvec[,2]
eigenvec <- eigenvec[,3:ncol(eigenvec)]
colnames(eigenvec) <- paste('Principal Component ', c(1:20), sep = '')


PED <- read.table("all_phase3.king.psam", header = TRUE, skip = 0, sep = '\t')
PED <- PED[which(PED$IID %in% rownames(eigenvec)), ]
PED <- PED[match(rownames(eigenvec), PED$IID),]



require('RColorBrewer')

PED$Population <- factor(PED$Population, levels=c(
  "ACB","ASW","ESN","GWD","LWK","MSL","YRI",
  "CLM","MXL","PEL","PUR",
  "CDX","CHB","CHS","JPT","KHV",
  "CEU","FIN","GBR","IBS","TSI",
  "BEB","GIH","ITU","PJL","STU"))

col <- colorRampPalette(c(
  "yellow","yellow","yellow","yellow","yellow","yellow","yellow",
  "forestgreen","forestgreen","forestgreen","forestgreen",
  "grey","grey","grey","grey","grey",
  "royalblue","royalblue","royalblue","royalblue","royalblue",
  "black","black","black","black","black"))(length(unique(PED$Population)))[factor(PED$Population)]

project.pca <- eigenvec
par(mar = c(5,5,5,5), cex = 2.0,
    cex.main = 7, cex.axis = 2.75, cex.lab = 2.75, mfrow = c(1,2))




plot(project.pca[,1], project.pca[,2],
     type = 'n',
     main = 'A',
     adj = 0.5,
     xlab = 'First component',
     ylab = 'Second component',
     font = 2,
     font.lab = 2)
points(project.pca[,1], project.pca[,2], col = col, pch = 20, cex = 2.25)
legend('bottomright',
       bty = 'n',
       cex = 3.0,
       title = '',
       c('AFR', 'AMR', 'EAS',
         'EUR', 'SAS'),
       fill = c('yellow', 'forestgreen', 'grey', 'royalblue', 'black'))

plot(project.pca[,1], project.pca[,3],
     type="n",
     main="B",
     adj=0.5,
     xlab="First component",
     ylab="Third component",
     font=2,
     font.lab=2)
points(project.pca[,1], project.pca[,3], col=col, pch=20, cex=2.25)



















































