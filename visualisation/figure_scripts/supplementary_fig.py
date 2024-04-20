#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 01:27:13 2024

@author: ashleyzhou
"""
#supplementary figure 

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
data = pd.read_csv('md_unv.csv', header=None)
# Sample data
categories = ['Task-repeat','Within-domain', 'Between-domain', 'Restart', 'Rest']
means = data.T.values[0]
errordata =pd.read_csv('md_unv_error_bar.csv', header=None)
errors = errordata.values[0]
means = np.insert(means,0,0)
errors = np.insert(errors, 0, 0)

allmd=pd.read_csv('md_unv_allsubs.csv',header=None)
allmd=allmd.values

data = pd.read_csv('dmpfc_unv.csv', header=None)
means_dmpfc = data.T.values[0]
means_dmpfc= np.insert(means_dmpfc, 0, 0)
errordata =pd.read_csv('dmpfc_unv_error_bar.csv', header=None)
errors_dmpfc = errordata.values[0]
errors_dmpfc = np.insert(errors_dmpfc,0,0)
alldmpfc=pd.read_csv('dmpfc_unv_allsubs.csv',header=None)
alldmpfc=alldmpfc.values


width = 0.25
# Custom bar colors for each bar
bar_colors = ['b', 'g', 'r', 'orange']
bar_colors = [ "#f05b43", "#c62320","#c62320","#831818"]
bar_colors2 = [ "#ffd353", "#ffb242", "#ffb242","#ef8737"]
bar_colors = [ "#ffb242","#ffb242", "#ffb242", "#ffb242","#ef8737"]
bar_colors2 = ["#97c684","#97c684","#97c684","#97c684", "#6f9969"]
r= np.arange(len(categories))
# Create the bar plot with error bars and custom colors
plt.bar(r, means, width=width, label='MD',yerr=errors, capsize=2, color='orange', alpha=0.7, edgecolor='black')

for i in range(4):
    plt.scatter([i+1]*len(allmd[i]),allmd[i], color='grey',edgecolor='black',alpha=0.3,zorder=2,s=10)
    #plt.scatter([i]*len(allsubs)+np.random.uniform(-0.05,0.05,len(allsubs)),allsubs[:,i], color='grey',edgecolor='black',alpha=0.5,zorder=2)
    
plt.errorbar(r, means, yerr=errors, fmt='none',ecolor='black',capsize=4,capthick=1,zorder=3)


plt.bar(r+width, means_dmpfc, width=width, label='dmPFC',yerr=errors_dmpfc, capsize=2, color='blue', alpha=0.7, edgecolor='black',zorder=1)

for i in range(4):
    plt.scatter([i+width+1]*len(alldmpfc[i]),alldmpfc[i], color='grey',edgecolor='black',alpha=0.3,zorder=2,s=10)
   

plt.errorbar(r+width, means_dmpfc, yerr=errors_dmpfc, fmt='none',capsize=2,  ecolor='black',capthick=1,zorder=3)

#significant_bars = [0,1, 2,3]  # Assuming bars 2 and 3 are significant
#for i in significant_bars:
#    plt.text(i, means[i] + errors[i]+0.1, '*', ha='center', va='bottom', fontsize=12, color='black')
plt.xticks(r+width,categories,fontsize=8)
# Add labels and title
#plt.ylim([0,1.4])
#plt.xlabel('Condition')
plt.ylabel('BOLD contrast')
#plt.title('Core DMN response to regular trials, compared to task-stay', fontsize=12)

# Add a legend in the top left corner
#plt.legend(['Core DMN'], loc='upper left', fontsize='small')

# Make the style fit for a science paper publication
#plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.tight_layout()
plt.legend( loc='upper left', fontsize='large',frameon=False)


# Set the highest possible resolution (e.g., 300 DPI)
dpi = 600
fig1=plt.gcf()
#fig1.savefig('core_mtl_dmn_high_res.jpg', dpi=dpi)

plt.show()