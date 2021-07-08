# link to admixtools https://cran.r-project.org/web/packages/admixr/vignettes/tutorial.html
##packages required
library(pophelper) #for guide on installation http://www.royfrancis.com/pophelper/index.html
library(stringr)
library(tools)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(colorspace)
library(dplyr)

tbl=read.table("results/ancestry/Q_files/PCA_ready.5.Q")
breeds <-read.csv("results/ancestry/clusters.csv", fill=T, stringsAsFactors = F) #dataframe with individual id and Breeds
#fam <- read.table("results/ancestry/PCA_ready.fam")
Benson<-c("Benson",1,"unknown")
breeds_with_Benson <- rbind(breeds, Benson)
mergedAdmixtureTable = cbind(tbl, breeds_with_Benson)

tbl2=read.table("results/ancestry/Q_files/PCA_ready.6.Q")
breeds <-read.csv("results/ancestry/clusters.csv", fill=T, stringsAsFactors = F) #dataframe with individual id and Breeds
#fam <- read.table("results/ancestry/PCA_ready.fam")
Benson<-c("Benson",1,"unknown")
breeds_with_Benson <- rbind(breeds, Benson)

mergedAdmixtureTable2 = cbind(tbl2, breeds_with_Benson)
oneindperbreed <- mergedAdmixtureTable2 %>% 
  group_by(CLUSTER) %>% 
  filter(row_number()==1)
#oneindperbreed6<- oneindperbreed[, -c(7:9)]
#write.table(oneindperbreed6, 'results/ancestry/short6.Q')
#file with just one individual per breed
oneindperbreed5 <- mergedAdmixtureTable %>% 
  group_by(CLUSTER) %>% 
  filter(row_number()==1)
#oneindperbreed5<- oneindperbreed5[, -c(6:8)]
#write.table(oneindperbreed5, 'results/ancestry/short5.Q')
ordered_5 = oneindperbreed5[order(oneindperbreed5$CLUSTER),] 
barplot(t(as.matrix(oneindperbreed5[,1:5])), col=rainbow(5),
        xlab="Individuals", ylab="Ancestry", border=NA)
barNaming <- function(vec) {
  retVec <- vec
  for(k in 2:length(vec)) {
    if(vec[k-1] == vec[k])
      retVec[k] <- ""
  }
  return(retVec)
}

par(mar=c(10,4,4,4))
barplot(t(as.matrix(subset(ordered_5, select=V1:V5))), col=rainbow(6), border=NA,
        names.arg=barNaming(ordered_5$CLUSTER), las=2)

#does not work
qlist<- list.files(path='./results/ancestry', pattern=".Q", full.names = T, recursive= F)
p2 <- plotQ(qlist=qlist,imgoutput="join",returnplot=T, ordergrp=T,showlegend=T, exportplot=F,basesize=11)
grid.arrange(p2$plot[[1]])
class(qlist)

##add two column to the ancestry table (ID and Breeds)
tbl$ID<- as.vector(breeds_with_Benson$ID)
tbl$Breeds <- as.factor(as.vector(breeds_with_Benson$CLUSTER))

##plotting
library(ggplot2)
library(forcats)
library(ggthemes)
library(patchwork)

##very simple bar plot
barplot(t(as.matrix(tbl[,1:5])), col=rainbow(5),
        xlab="Individuals", ylab="Ancestry", border=NA)

mergedAdmixtureTable = cbind(tbl, indTable)
mergedAdmWithPopGroups = merge(mergedAdmixtureTable, popGroups, by="Pop")

##ordered by breeds ans axis named by breeds
ordered = tbl[order(tbl$Breeds),] ##ordered by Breeds
barplot(t(as.matrix(subset(ordered, select=V1:V5))), col=rainbow(6), border=NA)

barNaming <- function(vec) {
  retVec <- vec
  for(k in 2:length(vec)) {
    if(vec[k-1] == vec[k])
      retVec[k] <- ""
  }
  return(retVec)
}

par(mar=c(10,4,4,4))
barplot(t(as.matrix(subset(ordered, select=V1:V5))), col=rainbow(6), border=NA,
        names.arg=barNaming(ordered$Breeds), las=2)

##try pop helper
# read ADMIXTURE Q files 

###Q list to use for analysis and data exploration
qlist<- readQ(list.files(path='results/ancestry/Q_files', pattern=".Q", full.names = T, recursive= F))
tr1 <- tabulateQ(qlist=qlist)
sr1 <- summariseQ(tr1)
is.qlist(qlist)
p1 <- plotQ(qlist[c(3,4)],imgoutput="join",returnplot=T,exportplot=F,basesize=11)
grid.arrange(p1$plot[[1]])

slist1 <- alignK(qlist[c(3,4)])
p2 <- plotQ(slist1,imgoutput="join",returnplot=T, ordergrp=T,showlegend=T, exportplot=F,basesize=11)
grid.arrange(p2$plot[[1]])

###
sfiles <- list.files(path=system.file("files/structure",package="pophelper"), full.names=T)
slist <- readQ(files=sfiles,indlabfromfile=T)

##this is the meta deta 
threelabset <- read.delim(system.file("files/metadata.txt", package="pophelper"), header=T,stringsAsFactors=F)
twolabset <- threelabset[,2:3]

plotQ(slist[2:3],imgoutput="join",showindlab=T,grplab=twolabset,
      subsetgrp=c("Brazil","Greece"),selgrp="loc",ordergrp=T,showlegend=T,
      showtitle=T,showsubtitle=T,titlelab="The Great Structure",
      subtitlelab="The amazing population structure of your favourite organism.",
      height=1.6,indlabsize=2.3,indlabheight=0.08,indlabspacer=-1,
      barbordercolour="white",barbordersize=0,outputfilename="plotq",imgtype="png",
      exportpath=getwd())###
##this one looks good!!!
#try to keep one individual per breed
##this one is very good but the name of the inviduals are quite close together
meta <-read.csv("results/ancestry/clusters.csv", fill=T, stringsAsFactors = F)
B<-c("Benson",1,"unknown")
metadata <- rbind(breeds, Benson)
metadata<- metadata %>% select(CLUSTER)
colnames(metadata) <-"Breeds" #change name to cluster column from cluster to breeds

p1 <- plotQ(alignK(qlist[c(1,2,5)]), imgoutput="join", showindlab = F, grplab=metadata, 
            selgrp="Breeds", ordergrp=T, grplabangle = 60, grplabpos = 0.7 , grplabheight = 2,
            grplabspacer = -0.7,
            showlegend=T, 
            legendtextsize = 5,
            height=1, #indlabsize=8, panelspacer=1, indlabheight=0.08,indlabspacer=10
            barbordercolour="white",barbordersize=0.5, returnplot=T,exportplot=F,
            basesize=11)
grid.arrange(p1$plot[[1]])
#summariseQ(tr1, writetable=TRUE)
#look in grp function to sort by group