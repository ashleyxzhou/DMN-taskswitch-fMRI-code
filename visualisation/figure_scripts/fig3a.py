
"""
Created on Mon Oct 23 18:43:05 2023
Context switch trials reaction time, figure3
@author: ashleyzhou
"""


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
data = pd.read_csv('rt_cs.csv', header=None)
# Sample data
categories = ['Task-repeat','Within-domain', 'Between-domain', 'Restart','Rest']
means = data.values[0]

errordata =pd.read_csv('rt_cs_error.csv', header=None)
errors = errordata.values[0]

allsubs=pd.read_csv('allsubs_rt_cs.csv', header=None)
allsubs= allsubs.values
#remove between subject variance
allsubs=allsubs.T - [np.nanmean(allsubs, axis=1)]*5 + np.ones((5,33))*np.nanmean(allsubs)
allsubs=allsubs.T

# Custom bar colors for each bar
bar_colors = ["#969bc7","#6b6ca3","#6b6ca3","#6b6ca3", "#434475"]

# Create the bar plot with error bars and custom colors
plt.bar(categories, means, yerr=errors, capsize=2, color=bar_colors,alpha=0.8, edgecolor='black',zorder=1)

#scatter individual points
for i in range(5):
    plt.scatter([i]*len(allsubs),allsubs[:,i], color='grey',edgecolor='black',alpha=0.3,zorder=2)
    #plt.scatter([i]*len(allsubs)+np.random.uniform(-0.05,0.05,len(allsubs)),allsubs[:,i], color='grey',edgecolor='black',alpha=0.5,zorder=2)
    
plt.errorbar(categories, means, yerr=errors, fmt='none',ecolor='black',capsize=5,capthick=1,zorder=3)
# Add labels and title
plt.ylim([0.45,1.2])
#plt.xlabel('Condition')
plt.ylabel('Mean RT of Background Detection Trials (s)')
#plt.title('Reaction time to context switch trials (jungle detection)', fontsize=12)
plt.xticks(fontsize=8)
plt.tight_layout()
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.axis()
# Set the highest possible resolution (e.g., 300 DPI)
dpi = 600

fig1 = plt.gcf()

# Show the plot (optional)
plt.show()