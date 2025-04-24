import copy
import numpy as np
from matplotlib import pyplot as plt
import random as r
from numba import jit
import networkx as nx

city_list = []
m = 1000
epsilon = 0.01
oracle_temp = 7000


def energy_to_neighbour(one:np.array, two:np.array):
    return np.sqrt(np.sum(np.power(np.subtract(one,two),np.array[2,2])))

def deepcopy(from_list:list, to_list:list):
    for i in range(len(from_list)):
        to_list[i] = from_list[i]

def move1(cities:list,length:int):
    ind1 = r.randint(0,length-1)
    ind2 = r.randint(0,length-1)
    while ind1 == ind2 and abs(ind2-ind1) != 1:
        ind2 = r.randint(0,length-1)
    cities[ind2], cities[ind1] = cities[ind1], cities[ind2]

def move2(length:int,cities:list):
    ind1 = r.randint(0,length-1)
    ind2 = r.randint(0,length-1)
    while ind1 == ind2 : #no self or direct neighbour interaction
        ind2 = r.randint(0,length-1)
    ind1 = min(ind1,ind2)
    ind2 = max(ind1,ind2)
    temp = cities[ind1:ind2]
    temp.reverse()
    for i in range(ind1,ind2):
        cities[i] = temp[i-ind1]

def total_energy(length, cities):
    ret = 0.0
    for i in range(1,length):
        ret += energy_to_neighbour(one=city_list[cities[i-1]],two=city_list[cities[i]])
    ret += energy_to_neighbour(one=city_list[cities[length-1]],two=city_list[cities[0]])
    return ret

def metropolis_step(tour:list, temp:float):
    """TODO: Calculate Energy of Configuration. Select new configuration. By either: Move 1: Swapping 2 non-adjacent indices or Move 2: selecting a subsequence and inverting it. Calculate Energy of Configuration. Do balancing
    """
    length = len(tour)
    E_old = total_energy(cities=tour,length=length)
    tour_old = list(range(length))
    deepcopy(from_list=tour,to_list=tour_old)
    move1(tour) if np.random.randint(0,100)%2 else move2(tour)
    E_new = total_energy(cities=tour,length=length)
    delta = E_new - E_old
    factor = np.exp(-delta/temp)
    r = np.random.uniform(-3000.0,3000.0)
    if delta < 0 or (r > 0  and r < factor):
        return tour
    else:
        return tour_old

def evolution_step(temp:float, tour:list):
    config_0 = list(range(len(tour)))
    deepcopy(from_list=tour,to_list=config_0)
    for i in range(m):
        tour = metropolis_step(tour=tour,temp=temp)
        if total_energy(length=len(config_0),cities=config_0) > total_energy(length=len(tour),cities=tour):
            deepcopy(from_list=tour,to_list=config_0)
    return config_0

def evolution(temp:float, tour:list):
    while temp > oracle_temp:
        zwischenergebnis = evolution_step(temp=temp,tour=tour)
        if total_energy(length=len(tour),cities=tour) > total_energy(length=len(zwischenergebnis),cities=zwischenergebnis):
            tour = zwischenergebnis
        temp *= (1-epsilon)
    print(tour)

def init_temp(length:int, tour:list,temp:float):
    config_0 = list(range(len(tour)))
    deepcopy(from_list=tour,to_list=config_0)
    max_temp = total_energy(cities=tour,length=length)
    for i in range(m):
        if r.randint(0,100)%2 : 
            move1(cities=tour,length=length)  
        else:
            move2(cities=tour,length=length)
        max_temp = max(max_temp, total_energy(cities=tour,length=length))
        if total_energy(length=len(config_0),cities=config_0) > total_energy(length=len(tour),cities=tour):
            deepcopy(from_list=tour,to_list=config_0)    

if __name__ == "__main__":
    file = "./ch130.tsp"
    x, y = np.loadtxt(file, delimiter=' ', comments="EOF", skiprows=6, usecols=(1, 2), unpack=True)
    n = len(x)
    for i in range(n):
        city_list.append(np.array([x[i],y[i]]))
    temp = 0.0
    tour = list(range(n))
    init_temp(length=n,temp=temp,tour=tour)
    evolution(temp=temp,tour=tour)


 