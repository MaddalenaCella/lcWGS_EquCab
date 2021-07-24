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
library(ggplot2)
library(hrbrthemes)

##packages to generate colour palette
##install them if they are missing
list.of.packages <- c("rstudioapi","fBasics","grDevices")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
##add packages to your library
library(rstudioapi)
library(fBasics)
library(grDevices)

##plot of K and CV
k_cv<- read.csv("analysis/Admix_CV.csv", fill=T, stringsAsFactors = F)
plot(k_cv)

## Load dataset from github
#data <- read.table("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/3_TwoNumOrdered.csv", header=T)
#data$date <- as.Date(data$date)

# Plot
pdf("results/ancestry/Cross_validation.pdf", width=7*0.5,height=7*0.5,pointsize=12*0.5)
ggplot(k_cv, aes(x=K., y=CV_error)) +
  geom_line(color="grey") +
  geom_point(shape=21, color="black", fill="#69b3a2", size=1) +
  xlab("K") + ylab("Cross-validation error") +
  theme_bw()
dev.off()
#ggsave(file="Cross_validation.pdf", plot = last_plot(), device = pdf, path = "results/ancestry")##save plot as pdf

breeds <-read.csv("results/ancestry/clusters.csv", fill=T, stringsAsFactors = F) #dataframe with individual id and Breeds
Benson<-c("Benson","BEN","unknown")
breeds_with_Benson <- rbind(breeds, Benson)

qlist1<- readQ(list.files(path='results/ancestry/Q_files/', pattern=".Q", full.names = T, recursive= F)) ##why can I not find all the Q files
meta <-read.csv("results/ancestry/clusters.csv", fill=T, stringsAsFactors = F)
B<-c("Benson","BEN","unknown")
metadata <- rbind(meta, B)
metadata<- metadata %>% select(CLUS)
colnames(metadata) <-"Breeds" #change name to cluster column from cluster to breeds

clist <- list(
  "shiny"=c("#1D4F9F","#BFF217","#60D5FD","#CC1577","#FF9326","#DF0101","#EC496F","#326397","#B26314","#027368","#A4A4A4","#610B5E"),
  "strong"=c("#11A4C8","#63C2C5","#1D4F9F","#0C516D","#2A2771","#396D35","#80C342","#725DA8","#B62025","#ED2224","#ED1943","#ED3995","#7E277C","#F7EC16","#F8941E","#8C2A1C","#808080"),
  "oceanfive"=c("#00A0B0", "#6A4A3C", "#CC333F", "#EB6841", "#EDC951"),
  "keeled"=c("#48B098", "#91CB62", "#FFEE3B", "#FB9013", "#FF3C28", "#000000"),
  "vintage"=c("#400F13", "#027368", "#A3BF3F", "#F2B950", "#D93A2B"),
  "muted"=c("#46BDDD","#82DDCE","#F5F06A","#F5CC6A","#F57E6A"),
  "teal"=c("#CFF09E","#A8DBA8","#79BD9A","#3B8686","#0B486B"),
  "merry"=c("#5BC0EB","#FDE74C","#9BC53D","#E55934","#FA7921"),
  "funky"=c("#A6CEE3", "#3F8EAA", "#79C360", "#E52829", "#FDB762","#ED8F47","#9471B4"),
  "retro"=c("#01948E","#A9C4E2","#E23560","#01A7B3","#FDA963","#323665","#EC687D"),
  "cb_paired"=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"),
  "cb_set3"=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F"),
  "morris"=c("#4D94CC","#34648A","#8B658A","#9ACD32","#CC95CC","#9ACD32","#8B3A39","#CD6601","#CC5C5B","#8A4500"),
  "wong"=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#006699","#D55E00","#CC79A7"),
  "krzywinski"=c("#006E82","#8214A0","#005AC8","#00A0FA","#FA78FA","#14D2DC","#AA0A3C","#FA7850","#0AB45A","#F0F032","#A0FA82","#FAE6BE"))

K4_6 <- plotQ(alignK(qlist1[c(6,7,8)], type="across"), imgoutput="join", showindlab = F, grplab=metadata, splab=c("K=4", "K=5", "K=6"),
            selgrp="Breeds", ordergrp = T, grplabangle = 270, grplabpos = 0.4 , grplabheight = 10,
            grplabspacer = -0.7, grplabsize = 2.5,
            showlegend=T, showdiv =F, clustercol=clist$shiny,
            legendtextsize = 6, legendkeysize = 6, legendpos="right",
            height=0.5, #indlabsize=8, panelspacer=1, indlabheight=0.08,indlabspacer=10
            barbordercolour="white",barbordersize=0.5, barsize = 0.7,
            returnplot=T,exportplot=F,
            outputfilename="plotq",imgtype="pdf",
            exportpath="results/ancestry",
            basesize=11)
grid.arrange(K4_6$plot[[1]])

K5_8 <- plotQ(alignK(qlist1[c(7,8,9,10)], type="across"), imgoutput="join", showindlab = F, grplab=metadata, splab=c("K=5", "K=6", "K=7", "K=8"),
              selgrp="Breeds", ordergrp = T, grplabangle = 90, grplabpos = 0.7 , grplabheight = 10,
              grplabspacer = -0.7, grplabsize = 2.5,
              showlegend=T, showdiv =F,
              legendtextsize = 5,
              height=0.5, #indlabsize=8, panelspacer=1, indlabheight=0.08,indlabspacer=10
              barbordercolour="white",barbordersize=0.5, barsize = 0.7,
              returnplot=T,exportplot=F,
              outputfilename="plotq",imgtype="pdf",
              exportpath="results/ancestry",
              basesize=11)
grid.arrange(K5_8$plot[[1]])
##all the breeds seem to change colour
##why dos the alignK not work :(
##look at best value of K

### PLOT CORRELATIONS TO SEE THE FIT
source("evalAdmix/visFuns.R") ## need to clone git repo locally

plot_evalAdmix <- function(Q_file, R_file){
pop <- as.vector(breeds_with_Benson$CLUSTER) # N length character vector with each individual population assignment
q <- as.matrix(read.table(Q_file)) # admixture porpotions q is optional for visualization but if used for ordering plot might look better
r <- as.matrix(read.table(R_file))
ord <- orderInds(pop=pop, q=q) # ord is optional but this make it easy that admixture and correlation of residuals plots will have individuals in same order
b<-plotAdmix(q=q, pop=pop, ord=ord)
a <- plotCorRes(cor_mat = r, pop = pop, ord=ord, title = "Admixture evaluation as correlation of residuals", max_z=0.25, min_z=-0.25,
                pop_labels = c(F,F), plot_legend = T)
return(a)
}
plot_evalAdmix("results/ancestry/Q_files/PCA_ready.5.Q","results/ancestry/corr/evaladmix5.corr")
plot_evalAdmix("results/ancestry/Q_files/PCA_ready.6.Q","results/ancestry/corr/evaladmix6.corr")
plot_evalAdmix("results/ancestry/Q_files/PCA_ready.9.Q","results/ancestry/corr/evaladmix9.corr")
plot_evalAdmix("results/ancestry/Q_files/PCA_ready.15.Q","results/ancestry/corr/evaladmix15.corr")

source("visFuns.R")

pop <- as.vector(breeds_with_Benson$CLUS) # N length character vector with each individual population assignment
q <- as.matrix(read.table("results/ancestry/Q_files/PCA_ready.5.Q")) # admixture porpotions q is optional for visualization but if used for ordering plot might look better
r <- as.matrix(read.table("results/ancestry/corr/evaladmix5.corr"))

ord <- orderInds(pop=pop, q=q) # ord is optional but this make it easy that admixture and correlation of residuals plots will have individuals in same order

plotAdmix(q=q, pop=pop, ord=ord)

graphics.off()
plotCorRes(cor_mat = r, pop = pop, ord=ord, title = "Admixture evaluation as correlation of residuals", max_z=0.25, min_z=-0.25, rotatelabpop=90, 
           cex.lab=0.4, cex.lab.2 = 0.4, plot_legend = T)

corr_list<- list.files(path='results/ancestry/corr', full.names = T, recursive= F)

#summariseQ(tr1, writetable=TRUE)
#look in grp function to sort by group