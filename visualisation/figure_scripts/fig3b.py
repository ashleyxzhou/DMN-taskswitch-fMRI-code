#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 18:37:51 2023
Scene memory recognition accuracy, figure 3
@author: ashleyzhou
"""


import matplotlib.pyplot as plt
import pandas as pd
data = pd.read_csv('recognition_accuracy.csv', header=None)
# Sample data
categories = ['Task-repeat','Within-domain', 'Between-domain', 'Restart', 'Rest']
means = data.values[0]
errordata =pd.read_csv('recognition_error.csv', header=None)
errors = errordata.values[0]
allsubs=pd.read_csv('allsubs_rec_accuracy.csv', header=None)
allsubs= allsubs.values

# Custom bar colors for each bar
bar_colors = ['pink','b', 'g', 'r', 'orange']
bar_colors = ["#c2d6a4", "#669d62","#669d62", "#669d62", "#1f5b25"]

# Create the bar plot with error bars and custom colors
plt.bar(categories, means, yerr=errors, capsize=2, color=bar_colors, alpha=0.7, edgecolor='black')

for i in range(5):
    plt.scatter([i]*len(allsubs),allsubs[:,i], color='grey',edgecolor='black',alpha=0.6,zorder=2)
    #plt.scatter([i]*len(allsubs)+np.random.uniform(-0.05,0.05,len(allsubs)),allsubs[:,i], color='grey',edgecolor='black',alpha=0.5,zorder=2)
    
plt.errorbar(categories, means, yerr=errors, fmt='none',ecolor='black',capsize=5,capthick=1,zorder=3)


# Add labels and title
plt.ylim([-1.25,2])
#plt.xlabel('Condition')
plt.ylabel("Scene Memory Recognition Accuracy (d')")
#plt.title(')
plt.xticks(fontsize=8)
# Add a legend in the top left corner
#plt.legend(['Core DMN'], loc='upper left', fontsize='small')

# Make the style fit for a science paper publication
#plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
# Set the highest possible resolution (e.g., 300 DPI)
dpi = 600

# bar_index1 = 0  # Index of the first bar to span
# bar_index2 = 3  # Index of the second bar to span
# vertical_line_y = max(means[bar_index1]+errors[bar_index1], means[bar_index2]+errors[bar_index2]) + 0.1
# plt.plot([bar_index1, bar_index2], [vertical_line_y, vertical_line_y], 'k-', alpha=0.5)

# bar_index1 = 2  # Index of the first bar to span
# bar_index2 = 4  # Index of the second bar to span
# vertical_line_y = max(means[bar_index1]+errors[bar_index1], means[bar_index2]+errors[bar_index2]) + 0.2
# plt.plot([bar_index1-0.5, bar_index2], [vertical_line_y, vertical_line_y], 'k-', alpha=0.5)
# plt.text((bar_index1+bar_index2)/2-0.25, vertical_line_y+0.05, '*', ha='center', va='bottom', fontsize=12, color='black')
fig1=plt.gcf()

# Export the figure to a .jpg format with the specified DPI
#fig1.savefig('rev_g_mem_high_res.jpg', dpi=dpi)

# Show the plot (optional)
plt.show()