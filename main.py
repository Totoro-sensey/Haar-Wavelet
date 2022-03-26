import matplotlib.pyplot as plt
import numpy as np


def result(x, y, x2, y2, N):
    Xr, Fk = get_Xr_Fk(y, N)
    Xr2, Fk2 = get_Xr_Fk(y2, N)
    fig, axs = plt.subplots(2, 2, figsize=(12, 7))
    axs[0][0].plot(x, y)
    axs[0][1].plot(np.repeat(x, 3), func_0(Xr))
    axs[1][0].plot(x2, y2)
    axs[1][1].plot(np.repeat(x, 3), func_0(Xr2))
    plt.show()

def map_color(mas_c21, mas_c22):
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(11, 6))
    plt.subplot(2, 1, 1)
    p1 = axs[0].pcolormesh(mas_c21, cmap='gray', edgecolors="face", shading='flat')
    plt.colorbar(p1, ax=axs[0])

    plt.subplot(2, 1, 2)
    p2 = axs[1].pcolormesh(mas_c22, cmap='gray', edgecolors="face", shading='flat')
    plt.colorbar(p2, ax=axs[1])

    plt.show()

def ito(xs1,ys1, ys2, yss1, yss2 ):
    fig, axs = plt.subplots(2, 2, figsize=(12, 7))
    axs[0][0].plot(xs1, ys1)
    axs[0][0].set_title('Исходный')
    axs[0][1].plot(xs1,ys2)
    axs[0][1].set_title('Восстановленный')
    axs[1][0].plot(xs1, yss1)
    axs[1][1].plot(xs1, yss2)
    plt.show()


def get_Xr_Fk(y, N):
    E = np.sum(y ** 2) / N
    k = np.arange(0, N)
    Re = np.arange(0, N)
    Im = np.arange(0, N)
    i = 0
    while i < N:
        Re[i] = np.sum(y*np.cos(2*np.pi*i*k/N))
        Im[i] = -np.sum(y*np.sin(2*np.pi*i*k/N))
        i += 1
    Re = Re/N
    Im = Im/N
    Xr = np.sqrt(Re * Re + Im * Im)
    Re[np.abs(Re) < E * 10**(-12)] = 0
    Im[np.abs(Im) < E * 10**(-12)] = 0
    Fk = np.arctan2(Im, Re)*180/np.pi
    return Xr, Fk

def func_0(x):
    return np.tile([0,1,0],len(x))*x.repeat(3)



def func_sigma2(x, b):
    mas = np.arange(0, x.size, 1)
    for i in range(x.size):
        if x[i] <= b:
            mas[i] = 0
        else:
            mas[i] = 1
    return mas


def res_func(y): #функция создает массив и заполняет каждый уровень передавая аргументы в функцию хаара
    mas_c = []
    mas_d = []
    iter = len(y)
    C = y
    mas_c.append(C)
    i=1
    while iter > 1:
        D = haara_D(C)
        C  = haara_func(C)
        mas_d.append(D)
        mas_c.append(C)
        i = i+1
        iter = iter / 2
    return mas_c, mas_d


def haara_func(y):
    C = np.arange(0, y.size/2, 1)
    j = 0
    i = 0
    while (i < len(y)):
        c1 = (y[i]+y[i+1])/(2**(0.5))
        C[j] = c1
        j = j+1
        i = i+2
    return C


def haara_D(y):
    D = np.arange(0, y.size/2)
    j = 0
    i = 0
    while (i < len(y)):
        d1 = (y[i]-y[i+1])/(2**(0.5))
        D[j] = d1
        j = j+1
        i = i+2
    return D

def createMas_c2(mas_c):
    mas_c2 = []
    C_2 = []
    for i in range(len(mas_c)):
        C = mas_c[i]
        for k in range(len(C)):
            el_c = C[k]
            j = 0
            while (j<2**i):
                C_2.append(el_c)
                j = j+1
        mas_c2.append(C_2)
        C_2 = []
    return mas_c2

def recovery_y(C, mas_d):
    mas_y = []
    mas_yr = []
    nowmas = [C]
    mas_y.append(float(C))
    iter = len(mas_d)
    while (iter > 0):
        iter = iter - 1
        row = mas_d[iter]
        i = 0
        while(i < (2**(len(mas_d) - iter -1))):
            c1 = (nowmas[i] + row[i])/(2**(0.5))
            c2 = (nowmas[i] - row[i]) / (2 ** (0.5))
            mas_yr.append(float(c1))
            mas_yr.append(float(c2))
            i = i+1
        mas_y.append(mas_yr)
        nowmas = mas_yr
        mas_yr = []
    return mas_y[7]

N = 128



#1
x = np.arange(0, N, 1)
y = np.sin((30 * np.pi * x)/N) + np.sin((20 * np.pi * x)/N)
y = np.asarray(y)
mas_c, mas_d = res_func(y)
mas_c21 = createMas_c2(mas_c)
mas_d21 = createMas_c2(mas_d)
y_now1 = recovery_y(mas_c[len(mas_c)-1], mas_d)

#2
b=100
x2 = np.arange(0, N, 1)
y2 = (func_sigma2(x2,0) * np.sin((30 * np.pi * x2)/N)) + (func_sigma2(x2,b) * np.sin((20 * np.pi * x2)/N))
result(x, y, x2, y2, N)
mas_c, mas_d = res_func(y2)
mas_c22 = createMas_c2(mas_c)
mas_d22 = createMas_c2(mas_d)
map_color(mas_c21, mas_c22)
map_color(mas_d21, mas_d22)
y_now2 = recovery_y(mas_c[len(mas_c)-1], mas_d)

ito(x, y, y_now1, y2, y_now2)

print(mas_d)


