"""
Práctica 3. TEOREMA DE LIOUVILLE

Alumna: Laura Cano Gómez
Subgrupo: U2 
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull, convex_hull_plot_2d
#https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html

#os.getcwd()




def deriv(q, dq0, d):  # Plantilla
    """
    ?????????????????
    
    Returns a ?????????????? object 

    Arguments:
        q   -> (variable de posición)
        dq0 -> (valor inicial de la derivada)
        d   -> (granularidad del parámetro temporal)
    """   
    #dq = np.empty([len(q)])
    dq = ( q[1 : len(q)] - q[0 : (len(q)-1)] ) / d
    dq = np.insert(dq, 0, dq0) #dq = np.concatenate(([dq0],dq))

    return dq

def F(q):  # Plantilla
    """
    Funcion oscilador no lineal pedida
    
    Returns a ????? object 

    Arguments:
        q   -> (variable de posición)
    """   
    k = 2
    ddq = - k * q * (q**2 -1)

    return ddq

#Resolución de la ecuación dinámica \ddot{q} = F(q), obteniendo la órbita q(t)
#Los valores iniciales son la posición q0 := q(0) y la derivada dq0 := \dot{q}(0)
def orb(n, q0, dq0, F, args=None, d=0.001):   # Plantilla
    """
    ?????????????????
    
    Returns a ?????????????? object 

    Arguments:
        q   -> (variable de posición)
        dq0 -> (valor inicial de la derivada)
        d   -> (granularidad del parámetro temporal)
    """  
    #q = [0.0]*(n+1)
    q = np.empty([n + 1])
    q[0] = q0
    q[1] = q0 + dq0 * d
    for i in np.arange(2, n + 1):
        args = q[i-2]
        q[i] = - q[i-2] + d**2 * F(args) + 2 * q[i-1]

    return q         




    
        
# Parte de areas, repasar
#Ejemplo de diagrama de fases (q, p) para un tiempo determinado



horiz = 0.251
# Parte del animate, y el horiz pasarlo como parametro
def make_gift(horiz):    # Plantilla sin nombre de funcion
    
    ax = fig.add_subplot(1,1, 1)
    seq_q0 = np.linspace(0.,1.,num=20)
    seq_dq0 = np.linspace(0.,2,num=20)
    q2 = np.array([])
    p2 = np.array([])
    for i in range(len(seq_q0)):
        for j in range(len(seq_dq0)):
            q0 = seq_q0[i]
            dq0 = seq_dq0[j]
            d = 10**(-3)
            n = int(horiz/d)
            #t = np.arange(n+1)*d
            q = orb(n,q0=q0,dq0=dq0,F=F,d=d)
            dq = deriv(q,dq0=dq0,d=d)
            p = dq/2
            q2 = np.append(q2,q[-1])
            p2 = np.append(p2,p[-1])
            plt.xlim(-2.2, 2.2)
            plt.ylim(-1.2, 1.2)
            plt.rcParams["legend.markerscale"] = 6
            ax.set_xlabel("q(t)", fontsize=12)
            ax.set_ylabel("p(t)", fontsize=12)
            ax.plot(q[-1], p[-1], marker="o", markersize= 10, markeredgecolor="red",markerfacecolor="red")
    
    return ax

        
        
        
        
        
        

#################################################################    
#  ESPACIO FÁSICO
################################################################# 

d = 10**(-4)

## Pintamos el espacio de fases
def simplectica(q0,dq0,F,col=0,d = 10**(-4),n = int(16/d),marker='-'): 
    q = orb(n,q0=q0,dq0=dq0,F=F,d=d)
    dq = deriv(q,dq0=dq0,d=d)
    p = dq/2
    plt.plot(q, p, marker,c=plt.get_cmap("winter")(col))


Horiz = 12
d = 10**(-3)

fig = plt.figure(figsize=(8,5))
fig.subplots_adjust(hspace=0.4, wspace=0.2)
ax = fig.add_subplot(1,1, 1)
#Condiciones iniciales:
seq_q0 = np.linspace(0.,1.,num=14)
seq_dq0 = np.linspace(0.,2,num=14)
for i in range(len(seq_q0)):
    for j in range(len(seq_dq0)):
        q0 = seq_q0[i]
        dq0 = seq_dq0[j]
        col = (1+i+j*(len(seq_q0)))/(len(seq_q0)*len(seq_dq0))
        #ax = fig.add_subplot(len(seq_q0), len(seq_dq0), 1+i+j*(len(seq_q0)))
        simplectica(q0=q0,dq0=dq0,F=F,col=col,marker='ro',d= 10**(-3),n = int(Horiz/d))
ax.set_xlabel("q(t)", fontsize=12)
ax.set_ylabel("p(t)", fontsize=12)
#fig.savefig('Simplectic.png', dpi=250)
plt.show()
