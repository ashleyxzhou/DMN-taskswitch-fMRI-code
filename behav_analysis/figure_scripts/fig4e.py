
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

plt.bar(np.arange(5),data['Novelty']['combined'][0:5], label='novco',yerr=data['Novelty']['combined_error'][0][0:5], capsize=2, color =colorbar, alpha=0.7, edgecolor='black')

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



dpi = 600
fig1=plt.gcf()
fig1.set_size_inches(6,5) 
fig1.savefig('combo_novelty.jpg', dpi=dpi)
