#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 16:12:19 2023
Decoding context, all Yeo paracellations, fig 4a
@author: ashleyzhou
"""
import matplotlib.pyplot as plt
import json
import numpy as np

with open('../Downloads/yeoDecodingResults.json', 'r') as json_file:
    data = json.load(json_file)
    
#categories = ['Task-repeat','Within-\ndomain', 'Between-\ndomain', 'Restart', 'Rest']

colorbar = ["#9a133d", "#b93961", "#d8527c", "#f28aaa", "#f9b4c9",
            "#6996e3","#6996e3","#6996e3","#6996e3","#6996e3","#6996e3","#6996e3","#6996e3","#6996e3","#6996e3","#6996e3","#6996e3",]

plt.barh(np.arange(17),data['Context']['sorted'], label='context_17',xerr=data['Context']['error'][0], capsize=2, color =colorbar, alpha=0.7, edgecolor='black')

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
desired_order = data['Context']['index']
sorted_bar_names = [bar_names[i-1] for i in desired_order]

plt.xlabel("Decoding accuracy for context (d')",fontsize=12)

y_pos = np.arange(17)
plt.yticks(y_pos, labels=sorted_bar_names, fontsize=12)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.tight_layout()
plt.gca().invert_yaxis()

plt.xlim([0,2.25])


dpi = 600
fig1=plt.gcf()
fig1.set_size_inches(7.24, 8) 
#fig1.savefig('yeo_context.jpg', dpi=dpi)
