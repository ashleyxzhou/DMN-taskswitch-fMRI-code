#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 14:16:55 2023
Decoding novelty for core and mtl, fig 4f
@author: ashleyzhou
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


#data = pd.read_csv('mtl_core_novelty_decoding.csv', header=None)
#means = data.values

errordata =pd.read_csv('novelty_mtl_core_error.csv', header=None)
errors = errordata.values

allsubs=pd.read_csv('allsubs_mtl_core_novelty_decoding.csv', header=None)
allsubs=allsubs.values

allmtl=allsubs[:,[0,2,4,6,8,10]]
#remove between subject variance
allmtl=allmtl.T - [np.nanmean(allmtl, axis=1)]*6 + np.ones((6,35))*np.nanmean(allmtl)
allmtl=allmtl.T
allmtl=allmtl[:,[0,1,2,4,3,5]]
allmtl[:,4]=np.mean(allmtl[:,4:6],1).tolist()

allcore=allsubs[:,[1,3,5,7,9,11]]
#remove between subject variance
allcore=allcore.T - [np.nanmean(allcore, axis=1)]*6 + np.ones((6,35))*np.nanmean(allcore)
allcore=allcore.T
allcore=allcore[:,[0,1,2,4,3,5]]
allcore[:,4]=np.mean(allcore[:,4:6],1).tolist()

means=[np.mean(allmtl,0).tolist(), np.mean(allcore,0).tolist()]
width = 0.3

#means = means[:,[0,1,2,4,3,5]]
#rest = np.mean(means[:][4:5],1)
#toplot=means[:][0:5]
#toplot[:,4]=rest.T


errors = errors[:,[0,1,2,4,3,5]]
reste = np.mean(errors[:,[4,5]],1)
errors_toplot=errors[:,0:5]
errors_toplot[:,4]=reste.T

bc2=["#007e2f","#007e2f","#007e2f","#007e2f"]
bc1=["#ffcd12","#ffcd12","#ffcd12","#ffcd12"]
bar_colors1 = [ "#ffb242","#ffb242", "#ffb242", "#ffb242","#ef8737"]
bar_colors2 = ["#97c684","#97c684","#97c684","#97c684", "#6f9969"]

bar_colors1 = [ "#ffd06f","#ffb242", "#ffb242", "#ffb242","#ef8737"]
bar_colors2 = ["#c2d6a4","#97c684","#97c684","#97c684", "#6f9969"]


categories = ['Task-repeat','Within-\ndomain', 'Between-\ndomain', 'Restart', 'Rest']
r = np.arange(len(categories))
# Create the bar plot with error bars and custom colors
plt.bar(r,means[1][0:5], width=width, label='Core DMN',yerr=np.abs(errors_toplot[1,:]), capsize=2, color=bar_colors1, alpha=1, edgecolor='black')

for i in range(5):
    plt.scatter([i]*len(allcore),allcore[:,i], color='grey',edgecolor='black',s=15,alpha=0.3,zorder=2)
    #plt.scatter([i]*len(allsubs)+np.random.uniform(-0.05,0.05,len(allsubs)),allsubs[:,i], color='grey',edgecolor='black',alpha=0.5,zorder=2)
    
plt.errorbar(categories, means[1][0:5], yerr=np.abs(errors_toplot[1,:]), fmt='none',ecolor='black',capsize=4,capthick=1,zorder=3)


plt.bar(r+width, means[0][0:5], width=width, label='MTL DMN',yerr=np.abs(errors_toplot[0,:]), capsize=2, color=bar_colors2, alpha=1, edgecolor='black')
for i in range(5):
    plt.scatter([i+width]*len(allmtl),allmtl[:,i], color='grey',edgecolor='black',s=15,alpha=0.3,zorder=2)
    #plt.scatter([i]*len(allsubs)+np.random.uniform(-0.05,0.05,len(allsubs)),allsubs[:,i], color='grey',edgecolor='black',alpha=0.5,zorder=2)
    
plt.errorbar(r+width, means[0][0:5], yerr=np.abs(errors_toplot[0,:]), fmt='none',ecolor='black',capsize=4,capthick=1,zorder=3)



plt.xticks(r+width,categories,fontsize=12)
plt.ylabel("Decoding accuracy for novelty (d')")
plt.legend( loc='upper right', fontsize='large',frameon=False)
plt.ylim([-1.8,2.5])
plt.xlim([-0.55,5])
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.tight_layout()

dpi = 600
fig1=plt.gcf()
fig1.set_size_inches(6,6) 
fig1.savefig('ndecode_novelty_mtl_core.jpg', dpi=dpi)

plt.show()