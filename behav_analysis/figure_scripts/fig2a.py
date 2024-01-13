

"""
Created on Thu Aug  3 17:13:21 2023
Reaction times for regular trials, figure 2
@author: ashleyzhou
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
data = pd.read_csv('figure_values/alltaskrt_reg.csv', header=None)
# Sample data
categories = ['Task-repeat','Within-domain', 'Between-domain', 'Restart','Rest']
means = data.values[0]

errordata =pd.read_csv('figure_values/alltaskrt_reg_error.csv', header=None)
errors = errordata.values[0]

errordata =pd.read_csv('reg_rt_error.csv', header=None)
errors = errordata.values[0]
errors = np.insert(errors, 0, 0)
errors=np.insert(errors,4,0)

# Custom bar colors for each bar
bar_colors = ['b', 'g', 'r']
bar_colors = ["#7db0ea", "#447fdd", "#225bb2"]
bar_colors = ["#7db0ea","#447fdd", "#447fdd", "#447fdd"]

# Create the bar plot with error bars and custom colors
plt.bar(categories, means, yerr=errors, capsize=2, color=bar_colors, alpha=0.7, edgecolor='black')

#significant_bars = [0, 1, 2]  # Assuming bars 2 and 3 are significant
#for i in significant_bars:
#    plt.text(i, means[i] + errors[i]+0.02, '*', ha='center', va='bottom', fontsize=12, color='black')
    
# Add labels and title
plt.ylim([1,1.4])
#plt.xlabel('Condition')
plt.ylabel('Mean RT for Regular Trials(s)')
#plt.title('Reaction time for regular trials, compared to task-stay', fontsize=12)
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

bar_index1 = 0  # Index of the first bar to span
bar_index2 = 4  # Index of the second bar to span
#vertical_line_y = max(means[bar_index1]+errors[bar_index1], means[bar_index2]+errors[bar_index2]) + 0.25
vertical_line_y = means[0]
# plt.plot([bar_index1, bar_index2], [vertical_line_y, vertical_line_y], 'k-', alpha=0.5)
fig1 = plt.gcf()
# bar_index1 = 0  # Index of the first bar to span
# bar_index2 = 3  # Index of the second bar to span
# vertical_line_y = max(means[bar_index1]+errors[bar_index1], means[bar_index2]+errors[bar_index2]) + 0.25
plt.plot([bar_index1, bar_index2], [vertical_line_y, vertical_line_y], 'k-', linestyle='dashed',alpha=0.5)
# plt.text((bar_index1+bar_index2)/2, vertical_line_y+0.05, '*', ha='center', va='bottom', fontsize=12, color='black')
plt.show()

# Export the figure to a .jpg format with the specified DPI
#plt.savefig('b_rt_reg_high_res.jpg', dpi=dpi)
#fig1.savefig('s_withtsrevised_rt_reg_high_res.jpg', dpi=dpi)

# Show the plot (optional)
