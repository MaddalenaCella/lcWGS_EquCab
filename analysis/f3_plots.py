#!/usr/bin/env python3

# Author: Maddalena Cella mc2820@ic.ac.uk
# Script: f3_plots.py
# Description: python matplotlib plots of f3 statistics results
# Date: 27 July 2021

""" f3 plots with matplotlib """
__appname__="f3_plots.py"
__author__="Maddalena Cella mc2820@ic.ac.uk"
__version__="0.0.1"
__license__="None"

import pandas as p
import matplotlib.pyplot as plt

f3_table=pd.read_csv("results/ancestry/eigenstrat/modified_f3.csv", index_col= False)


d=f3_table.sort_values(by="F3")
y = range(len(d))
plt.figure(figsize=(14, 10))
c1= ["#50b041", "#50b041", "#50b041", "#50b041", "#50b041", "#aec7aa", "#aec7aa", "#aec7aa", "#aec7aa", "#aec7aa"]
plt.errorbar(d["F3"], y, xerr=d["SE"], fmt='o', mec='black', mfc= 'black', ecolor= c1)
plt.yticks(y, d["Populations"]);
plt.xlabel("F3(Pop A - Pop B - Target)");
#plt.show()
plt.savefig('results/ancestry/eigenstrat/horizontal_f3_plot.pdf')


f3_table_out=pd.read_csv("results/ancestry/eigenstrat/modified_f3_outgroup.csv", index_col= False)

p=f3_table_out.sort_values(by="F3")
y = range(len(p))
plt.figure(figsize=(14, 10))
c= ["#aec7aa", "#aec7aa", "#aec7aa", "#aec7aa", "#50b041", "#50b041", "#50b041", "#50b041", "#50b041", "#50b041" ]
plt.errorbar(p["F3"], y, xerr=p["SE"], fmt='o', mec='black', mfc= 'black', ecolor= c)
plt.yticks(y, p["Populations"]);
plt.xlabel("F3(Target - Pop B - Outgroup)");
#plt.show()
plt.savefig('results/ancestry/eigenstrat/horizontal_f3_outgroup.pdf')

