#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 13:12:48 2023
Context decoding for Core and MTL, fig 4c
@author: ashleyzhou
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data = pd.read_csv('mtl_core_context_decoding.csv', header=None)
means = data.values

errordata =pd.read_csv('context_mtl_core_error.csv', header=None)
errors = errordata.values

width = 0.35

means = means[:,[0,1,2,4,3,5]]
rest = np.mean(means[:,[4,5]],1)
toplot=means[:,0:5]
toplot[:,4]=rest.T


errors = errors[:,[0,1,2,4,3,5]]
reste = np.mean(errors[:,[4,5]],1)
errors_toplot=errors[:,0:5]
errors_toplot[:,4]=reste.T

bc2=["#007e2f","#007e2f","#007e2f","#007e2f"]
bc1=["#ffcd12","#ffcd12","#ffcd12","#ffcd12"]
bar_colors1 = [ "#ffd06f","#ffb242", "#ffb242", "#ffb242","#ef8737"]
bar_colors2 = ["#c2d6a4","#97c684","#97c684","#97c684", "#6f9969"]

categories = ['Task-repeat','Within-\ndomain', 'Between-\ndomain', 'Restart', 'Rest']
r = np.arange(len(categories))
# Create the bar plot with error bars and custom colors
plt.bar(r,toplot[1,:], width=width, label='Core DMN',yerr=np.abs(errors_toplot[1,:]), capsize=2, color=bar_colors1, alpha=0.7, edgecolor='black')

plt.bar(r+width, toplot[0,:], width=width, label='MTL DMN',yerr=np.abs(errors_toplot[0,:]), capsize=2, color=bar_colors2, alpha=0.7, edgecolor='black')

plt.xlim([-0.55,5.55])
plt.xticks(r+width,categories,fontsize=8)
plt.ylabel("Decoding accuracy for context (d')")
plt.legend( loc='upper right', fontsize='large',frameon=False)

plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.tight_layout()

dpi = 600
fig1=plt.gcf()
#fig1.savefig('decode_context_mtl_core.jpg', dpi=dpi)


