#!/usr/bin/Rscript
# Author: Maddalena Cella
# Email mc2820@ic.ac.uk
# Date: 01-07-2021
# Last Modified: 01-07-2021
# Desc: script that plots eigenvector values from plink PCA

pcs = read.table( "../analysis/PCA_analysis.eigenvec" )
View(pcs)
plot( pcs[,3], pcs[,4], xlab = "PC1", ylab = "PC2" )