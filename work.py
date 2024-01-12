import intvalpy as ip
import numpy as np
import matplotlib.pyplot as plt
from math import isnan
from copy import deepcopy


def create_L(X):
    midL = np.zeros((2, 2))
    radL = np.zeros((2, 2))
    radL[0][0] = 0
    midL[0][0] = 1

    radL[0][1] = 0
    midL[0][1] = 1

    radL[1][0] = 0.5*(1/X.a[1] - 1/X.b[1])
    midL[1][0] =  1/X.b[1] + radL[1][0]

    
    interval = -X[0]/(X[1]**2)
    radL[1][1] = 0.5*(interval.b - interval.a)
    midL[1][1] = interval.a + radL[1][1]
    L = ip.Interval(midL, radL, midRadQ=True)
    return L

def inv_midL(L):
    midL = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            midL[i][j] = L.mid[i][j]
    return np.linalg.inv(midL)

def eig_midL(L):
    midL = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            midL[i][j] = L.mid[i][j]
    return np.linalg.eig(midL)

def F_x(x):
    B = ip.Interval([[1, 3], [1, 6]])
    F = ip.Interval([[0, 0], [0, 0]])

    F[0] = x[0] + x[1] - B[0]
    F[1] = x[0]/x[1] - B[1]
    return F

def kravchik(x_av, F, Lambda, I, X, L):
    return x_av - Lambda @ F + (I - Lambda @ L) @ (X - x_av)

def intersection(X, kr):
    #inter = ip.intersection(X, kr)
    # if isnan(inter.a[0]):
    #     inter[0] = X[0]
    # if isnan(inter.a[1]):
    #     inter[1] = X[1]
        
    return ip.intersection(X, kr)

if __name__ == "__main__":
    I = [[1,0], [0,1]]
    X = ip.Interval([[0.1, 2.7], [2, 8]])  
    iter = 0

    X_k = []
    X_k.append(X)
    i = 0

    kr_k = []

    num_iter = 50
    
    print("Количевство итераций:", i , "\n", X_k[i][0], X_k[i][1], )
    for i in range (1, num_iter):
        L = create_L(X)
        Lambda = inv_midL(L)
        q = eig_midL(I - Lambda @ L)
      #  print(I - Lambda @ L)
    

        x_av = X.mid
        F = F_x(x_av)
        kr = kravchik(x_av, F, Lambda, I, X, L)
        kr_k.append(kr)
        
        X_k.append(intersection(X, kr))
      

        print("Кравчик = ", kr_k[i-1][0], kr_k[i-1][1] )
        print("Количевство итераций:",i, "\n", "X_k = ", X_k[i][0], X_k[i][1], )
        X = X_k[i]
        iter += 1
        inter = ip.intersection(X, kr)
        if isnan(inter.a[0]):
            break
        if isnan(inter.a[1]):
            break
        

    kr_k.append(kravchik(x_av, F, Lambda, I, X, L))

    fig, ax = plt.subplots(figsize=(4, 4))


    for i in range(iter):
        kr_one = abs(kr_k[i][0].b - kr_k[i][0].a) 
        kr_two = abs(kr_k[i][1].b - kr_k[i][1].a) 
        one = abs(X_k[i][0].b - X_k[i][0].a)
        two = abs(X_k[i][1].b - X_k[i][1].a)
        Rect = plt.Rectangle((kr_k[i][0].a, kr_k[i][1].a), kr_one , kr_two, edgecolor='red', facecolor='none', linewidth=0.7)
        
        iveRect = plt.Rectangle((X_k[i][0].a, X_k[i][1].a), one , two, edgecolor='black', facecolor='none', linewidth=1.5)
      
        plt.gca().add_patch(iveRect)
        plt.gca().add_patch(Rect)

    x = np.arange(-10, 15)

    y = 1 - x
    plt.plot(x, y, color='g', linewidth=0.7, label = 'x + y = 1')

    y = 3 - x
    plt.plot(x, y, color='olive', linewidth=0.7, label = 'x + y = 3')

    y = x
    plt.plot(x, y, color='blue', linewidth=0.7, label = 'x/y = 1')
    
    y = x/6
    plt.plot(x, y, color='darkblue', linewidth=0.7, label = 'x/y = 6')


    plt.grid()
    # plt.xlim(-5, 3)
    # plt.ylim(-2, 5.5)
    plt.xlim(-4, 10)
    plt.ylim(-8, 10)
    plt.legend()
    plt.show()
    