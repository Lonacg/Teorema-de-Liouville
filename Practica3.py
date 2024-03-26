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
from matplotlib import animation
import time


#os.getcwd()



#################################################################    
#  APARTADO II)
################################################################# 


def orb(n, q0, dq0, F, args=None, d=0.001):   # Plantilla
    """
    Obtiene la orbita de la ecuacion dinamica
    
    Returns a array object (q)

    Arguments:
        n   -> int (numero de variables de estado)
        q0  -> float (variable de s)
        dq0 -> float (valor inicial de la derivada)
        F   -> function (funcion oscilador no lineal)
        d   -> float (granularidad del parámetro temporal)
    """  
    q = np.empty([n + 1])
    q[0] = q0
    q[1] = q0 + dq0 * d
    for i in np.arange(2, n + 1):
        args = q[i-2]
        q[i] = - q[i-2] + d**2 * F(args) + 2 * q[i-1]

    return q         


def deriv(q, dq0, d):  # Plantilla
    """
    Returns an array object (dq)

    Arguments:
        q   -> array (variable de posición)
        dq0 -> array (valor inicial de la derivada)
        d   -> float (granularidad del parámetro temporal)
    """   
    dq = ( q[1 : len(q)] - q[0 : (len(q)-1)] ) / d
    dq = np.insert(dq, 0, dq0) 

    return dq


def F(q):  # Plantilla
    """
    Funcion oscilador no lineal pedida
    
    Returns a array object 

    Arguments:
        q   -> array (variable de posición)
    """   
    k = 2
    ddq = - k * q * (q**2 -1)

    return ddq


def area_aux(t, seq_q0, seq_dq0, d, title, show_it= True): # Plantilla
    """
    Calcula el area de la componente conexa

    Returns a list object:
        q2 ->
        p2 ->

    Arguments:
        seq_q0   -> list (numero de variables de estado)
        seq_dq0  -> float (variable de s)
        d -> float (granularidad del parámetro temporal)
        title   -> string (titulo de la grafica)

    """ 
    if show_it: 
        fig, ax = plt.subplots(figsize=(12,5))
        ax = fig.add_subplot(1, 1, 1)
    
    horiz = t
     
    q2 = np.array([])
    p2 = np.array([])
    for i in range(len(seq_q0)):
        for j in range(len(seq_dq0)):
            q0 = seq_q0[i]
            dq0 = seq_dq0[j]
            n = int(horiz/d)
            q = orb(n, q0= q0, dq0= dq0, F= F, d= d)
            dq = deriv(q, dq0= dq0, d= d)
            p = dq/2
            q2 = np.append(q2,q[-1])
            p2 = np.append(p2,p[-1])
            if show_it:
                plt.xlim(-2.2, 2.2)
                plt.ylim(-1.2, 1.2)
                plt.rcParams["legend.markerscale"] = 6
                ax.set_xlabel("q(t)", fontsize=12)
                ax.set_ylabel("p(t)", fontsize=12)
                plt.title(title)
                plt.plot(q[-1], p[-1], marker="o", markersize= 10, markeredgecolor="red",markerfacecolor="red")
    

    X = np.array([q2,p2]).T
    hull = ConvexHull(X)

    if show_it:
        #plt.show()   
        convex_hull_plot_2d(hull)

    return(hull.volume)


def calculate_area(d, t= 1/3, show_it= True):
    """
    Calcula el area total aplicando area_aux a cada componente convexa

    Returns a float object (the area)

    Arguments:
        d -> float (granularidad del parámetro temporal)

    """      
    polygonal_area = area_aux(t, np.linspace(0., 1., num=20), np.linspace(0., 2., num=20), d, "Polygonal Area", show_it)
    lower_area = area_aux(t, np.linspace(0., 1., num=40), [0], d, "Lower Area", show_it)
    right_area = area_aux(t, [1], np.linspace(0., 2., num=20), d, "Right Area", show_it)

    area = polygonal_area - (lower_area + right_area)
    
    return area
         


# Espacio de fases:
'''

d = 10**(-4)
## Pintamos el espacio de fases
def simplectica(q0, dq0, F, col= 0, d= 10**(-4), n = int(16/10**(-4)),marker='-'): 
    q = orb(n,q0=q0,dq0=dq0,F=F,d=d)
    dq = deriv(q,dq0=dq0,d=d)
    p = dq/2
    plt.plot(q, p, marker,c=plt.get_cmap("winter")(col))


Horiz = 12
d = 10**(-3)

fig = plt.figure(figsize=(8,5))
fig.subplots_adjust(hspace= 0.4, wspace= 0.2)
ax = fig.add_subplot(1, 1, 1)
#Condiciones iniciales:
seq_q0 = np.linspace(0., 1., num= 14)
seq_dq0 = np.linspace(0., 2, num= 14)
for i in range(len(seq_q0)):
    for j in range(len(seq_dq0)):
        q0 = seq_q0[i]
        dq0 = seq_dq0[j]
        col = (1+i+j*(len(seq_q0)))/(len(seq_q0)*len(seq_dq0))
        #ax = fig.add_subplot(len(seq_q0), len(seq_dq0), 1+i+j*(len(seq_q0)))
        simplectica(q0=q0, dq0=dq0, F=F, col=col, marker='ro',d= 10**(-3), n = int(Horiz/d))
ax.set_xlabel("q(t)", fontsize=12)
ax.set_ylabel("p(t)", fontsize=12)
#fig.savefig('Simplectic.png', dpi=250)
plt.show()





horiz = 0.251
# Parte del animate, y el horiz pasarlo como parametro
def make_gift(horiz):    # Plantilla sin nombre de funcion
    ax = fig.add_subplot(1, 1, 1)
    seq_q0 = np.linspace(0., 1., num=20)
    seq_dq0 = np.linspace(0., 2, num=20)
    q2 = np.array([])
    p2 = np.array([])
    for i in range(len(seq_q0)):
        for j in range(len(seq_dq0)):
            q0 = seq_q0[i]
            dq0 = seq_dq0[j]
            d = 10**(-3)
            n = int(horiz/d)
            #t = np.arange(n+1)*d
            q = orb(n, q0=q0, dq0=dq0, F=F, d=d)
            dq = deriv(q, dq0=dq0, d=d)
            p = dq/2
            q2 = np.append(q2, q[-1])
            p2 = np.append(p2, p[-1])
            plt.xlim(-2.2, 2.2)
            plt.ylim(-1.2, 1.2)
            plt.rcParams["legend.markerscale"] = 6
            ax.set_xlabel("q(t)", fontsize=12)
            ax.set_ylabel("p(t)", fontsize=12)
            ax.plot(q[-1], p[-1], marker="o", markersize= 10, markeredgecolor="red",markerfacecolor="red")
    
    return ax

        
        


def init():
    return make_gift(.1)


fig = plt.figure(figsize=(6,6))
ani = animation.FuncAnimation(fig, make_gift, np.arange(0.1,5.1,0.1), init_func=init, interval=60)
ani.save("Evolucion_D(t,inf).gif", fps = 10) 

        
'''


def main():



    '''
    APARTADO ii) 
    - Obtén el valor del área de D para t = 1/3 y una estimación del su intervalo 
    de error, presentando los valores de forma científicamente formal (puedes 
    estimar el error a partir de la sensibilidad al parámetro δ). 
    - ¿Se cumple el teorema de Liouville entre D0 y Dt? Razona la respuesta.
    '''
    
    d1 = 10**(-3)
    d2 = 10**(-4)
    t = 1/3

    print('APARTADO II)')
    a1 = calculate_area(d1, t)
    a2 = calculate_area(d2, t, show_it= False)


    print("El area para t= 1/3 es:")
    print(f'    Para d = 10**(-3): {a1}')
    print(f'    Para d = 10**(-4): {a2}')



    print(f'\nEstimacion del error: {abs(a1-a2)}')

    d_range = np.linspace(d1, d2, num=10)
    areas=[]
    for d in d_range:
        a = calculate_area(d, t, show_it= False)
        areas.append(a)
    
    print(areas)

    sum = 0
    for a in areas:
        sum += a
    print(1- sum/len(areas))
    


    '''

    
    print(f'Los codigos de Huffman son:\n - SEng:  {code_en}')
    print(f' - SEsp:  {code_es} \n')
    
    [mean_en , error_length_en] = average_length(distr_en, text_en)
    [mean_es , error_length_es] = average_length(distr_es, text_es)

    print('La longitud media es:\n')
    print(f' - L(SEng):  {mean_en} +- {error_length_en}\n')
    print(f' - L(SEsp):  {mean_es} +- {error_length_es}\n') 

    print('Para comprobar si satisface el Teorema de Shannon calculamos la entropia de cada idioma:')
    
    [entropy_en, error_entropy_en] = entropy(distr_en, text_en)
    [entropy_es, error_entropy_es] = entropy(distr_es, text_es)

    print(f'La entropia es:\n')
    print(f' - H(SEng):  {entropy_en} +- {error_entropy_en}\n')
    print(f' - H(SEsp):  {entropy_es} +- {error_entropy_es}\n')     
    
    print('Luego, para satisfacer el T. de Shannon, debe cumplir  H(S) <= L(S) < H(S) + 1.')
    print(f'Ingles:  ¿ {entropy_en} <= {mean_en} < {entropy_en + 1} ? -> { (entropy_en <= mean_en) and (mean_en < entropy_en + 1) }')
    print(f'Español: ¿ {entropy_es} <= {mean_es} < {entropy_es + 1} ? -> { (entropy_es <= mean_es) and (mean_es < entropy_es + 1) }')



    print('\n\nAPARTADO II) \nCodificacion de la palabra "Lorentz":')

    code_lorentz_en = encode('Lorentz', code_en)
    code_lorentz_es = encode('Lorentz', code_es)

    print(f'Ingles:  {code_lorentz_en}')
    print(f'Español: {code_lorentz_es}')

    print('\nComprobacion de la eficiencia de longitud frente al código binario usual para la palabra "Lorentz":')

    bits = 8 * len('Lorentz')                               # En el código ASCII cada caracter ocupa log2(256) = 8bits.

    print(f'Longitud del codigo binario usual: {bits}')
    print(f'Longitud del codigo de Huffman en el idioma ingles: {len(code_lorentz_en)}')
    print(f'Longitud del codigo de Huffman en el idioma español: {len(code_lorentz_es)}')



    print('\n\nAPARTADO III) \nDecodificacion del resultado anterior:')

    decode_lorentz_en = decode(code_lorentz_en, code_en)
    decode_lorentz_es = decode(code_lorentz_es, code_es)

    print(f'Ingles  ({code_lorentz_en}): {decode_lorentz_en}  ')
    print(f'Español ({code_lorentz_es}): {decode_lorentz_es}')    

    '''



if __name__ == '__main__':
    main()






