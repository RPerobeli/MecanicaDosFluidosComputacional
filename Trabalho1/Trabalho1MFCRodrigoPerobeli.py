# -*- coding: utf-8 -*-

import numpy as np
import matplotlib as plt
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits import mplot3d
    
    
#MAIN
'''
    Algoritmo que le um arquivo texto e plota a solução do trabalho
'''
delta_t = 0.1
deltaX = 0.1


arquivo =  open("MatrizSolucaoTrabalho1.txt",'r')
i= int(arquivo.readline())
j= int(arquivo.readline())

matriz = np.zeros([i,j])
matriz = np.array(matriz)
for c in range(matriz.shape[0]):
    for d in range(matriz.shape[1]):
        valor = float(arquivo.readline())
        matriz[c,d] = valor
    
arquivo.close()

t = np.zeros(i)
for c in range(i-1):
    t[c+1]= t[c]+delta_t
X = np.zeros(j)
for c in range(j-1):
    X[c+1] = X[c]+deltaX


fig=plt.figure.Figure()
ax.plot_surface(X, t, matriz, rstride=5, cstride=5,
cmap='viridis', edgecolor='none')
ax.set_title('surface');



    



    