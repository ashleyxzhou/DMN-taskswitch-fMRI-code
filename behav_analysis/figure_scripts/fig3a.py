
"""
Created on Mon Oct 23 18:43:05 2023
Context switch trials reaction time, figure3
@author: ashleyzhou
"""


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
data = pd.read_csv('figure_values/rt_cs.csv', header=None)
# Sample data
categories = ['Task-repeat','Within-domain', 'Between-domain', 'Restart','Rest']
means = data.values[0]

errordata =pd.read_csv('figure_values/rt_cs_error.csv', header=None)
errors = errordata.values[0]

# Custom bar colors for each bar
bar_colors = ['b', 'g', 'r']
bar_colors = ["#7db0ea", "#447fdd", "#225bb2"]

bar_colors = ["#969bc7","#6b6ca3","#6b6ca3","#6b6ca3", "#434475"]

# Create the bar plot with error bars and custom colors
plt.bar(categories, means, yerr=errors, capsize=2, color=bar_colors, alpha=0.7, edgecolor='black')

# Add labels and title
plt.ylim([0.75,0.9])
#plt.xlabel('Condition')
plt.ylabel('Mean RT of Context Switch Trials(s)')
#plt.title('Reaction time to context switch trials (jungle detection)', fontsize=12)
plt.xticks(fontsize=8)
# Add a legend in the top left corner
#plt.legend(['Core DMN'], loc='upper left', fontsize='small')

# Make the style fit for a science paper publication
#plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.axis()
# Set the highest possible resolution (e.g., 300 DPI)
dpi = 300

fig1 = plt.gcf()

# Export the figure to a .jpg format with the specified DPI
#plt.savefig('b_rt_reg_high_res.jpg', dpi=dpi)
#fig1.savefig('revised_rt_cs_high_res.jpg', dpi=dpi)

# Show the plot (optional)
