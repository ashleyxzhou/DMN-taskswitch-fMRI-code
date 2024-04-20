
"""
Created on Fri Nov 24 16:12:19 2023
Decoding context, all significant parcellations, fig 4e

@author: ashleyzhou
"""
import matplotlib.pyplot as plt
import json
import numpy as np

with open('../Downloads/yeoDecodingResults.json', 'r') as json_file:
    data = json.load(json_file)
    
categories = ['Task-repeat','Within-\ndomain', 'Between-\ndomain', 'Restart', 'Rest']

colorbar = [ "#f9b4c9", "#d8527c","#d8527c","#d8527c","#9a133d",]

means=(data['Novelty']['combined'][0:6])
means=np.array(means)
means=means[[0,1,2,4,3,5]]
rest=np.mean(means[[4,5]])
toplot=means[0:5]
toplot[4]=rest.T

errors=data['Novelty']['combined_error'][0][0:6]
errors=np.array(errors)
errors=errors[[0,1,2,4,3,5]]
reste=np.mean(errors[[4,5]])
errors_toplot=errors[0:5]
errors_toplot[4]=reste.T

allsubs=data['Novelty']['all']
temp=np.array(allsubs)
allsubs=temp.reshape(-1,6)
#allsubs=temp_mat.tolist()
allsubs.tolist()

#remove between subject variance
allsubs=allsubs.T - [np.nanmean(allsubs, axis=1)]*6 + np.ones((6,35))*np.nanmean(allsubs)
allsubs=allsubs.T

plt.bar(np.arange(5),toplot, label='nov_co',yerr=errors_toplot, capsize=2, color =colorbar, alpha=0.8, edgecolor='black')

for i in range(5):
    plt.scatter([i]*len(allsubs),allsubs[:,i], color='grey',edgecolor='black',alpha=0.3,zorder=2)
    #plt.scatter([i]*len(allsubs)+np.random.uniform(-0.05,0.05,len(allsubs)),allsubs[:,i], color='grey',edgecolor='black',alpha=0.5,zorder=2)
    
plt.errorbar(categories, toplot, yerr=errors_toplot, fmt='none',ecolor='black',capsize=5,capthick=1,zorder=3)


bar_names =['Visual-C.',
'Visual-Peripheral',
'Somatosensory-Motor A',
'Somatosensory-Motor B',
'Dorsal-Attention A',
'Dorsal-Attention B',
'Ventral-Attention A',
'Ventral-Attention B',
'Limbic A',
'Limbic B',
'Control C',
'Control A',
'Control B',
'Temporal-Parietal',
'DMN-MTL',
'DMN-Core',
'DMN-dMPFC'
]
desired_order = data['Novelty']['index']
sorted_bar_names = [bar_names[i-1] for i in desired_order]

plt.xticks(np.arange(5),categories,fontsize=12)

plt.ylabel("Decoding accuracy for novelty (d')",fontsize=12)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.tight_layout()
plt.ylim([-1.5,1.5])


dpi = 600
fig1=plt.gcf()
fig1.set_size_inches(6,5) 
fig1.savefig('combo_novelty.jpg', dpi=dpi)
