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

class City:
    def __init__(self, pos:np.array):
        self.pos = pos
    def __repr__(self):
        return f"{self.pos}"
    def energy_to_neighbour(self,neighbour): # type: ignore
        return np.sqrt(np.sum(np.power(np.subtract(neighbour.pos,self.pos),np.array([2.0,2.0]))))

class Tour:
    def __init__(self, cities:list[int]):
        self.length = len(cities)
        self.cities = cities
    def __repr__(self):
        return f"{self.cities}"

    def move1(self):
        ind1 = r.randint(0,self.length-1)
        ind2 = r.randint(0,self.length-1)
        while ind1 == ind2 and abs(ind2-ind1) != 1:
            ind2 = r.randint(0,self.length-1)

        self.cities[ind2], self.cities[ind1] = self.cities[ind1], self.cities[ind2]

    def move2(self):
        ind1 = r.randint(0,self.length)
        ind2 = r.randint(0,self.length)
        while ind1 == ind2 : #no self or direct neighbour interaction
            ind2 = r.randint(0,self.length)
        ind1 = min(ind1,ind2)
        ind2 = max(ind1,ind2)
        temp = self.cities[ind1:ind2]
        temp.reverse()
        for i in range(ind1,ind2):
            self.cities[i] = temp[i-ind1]

    def total_energy(self):
        ret = 0.0
        for i in range(1,self.length):
            ret+=city_list[self.cities[i-1]].energy_to_neighbour(city_list[self.cities[i]])
        ret += city_list[self.cities[self.length-1]].energy_to_neighbour(city_list[self.cities[0]])
        return ret

def metropolis_step(tour:Tour, temp:float):
    """TODO: Calculate Energy of Configuration. Select new configuration. By either: Move 1: Swapping 2 non-adjacent indices or Move 2: selecting a subsequence and inverting it. Calculate Energy of Configuration. Do balancing
    """
    E_old = tour.total_energy()
    tour_old = Tour(tour.cities)
    tour.move1() if np.random.randint(0,100)%2 else tour.move2()
    E_new = tour.total_energy()
    delta = E_new - E_old
    factor = np.exp(-delta/temp)
    r = np.random.uniform(-3000.0,3000.0)
    if delta < 0 or (r > 0  and r < factor):
        return tour
    else:
        return tour_old

def evolution_step(temp:float, tour:Tour):
    config_0 = Tour(tour.cities)
    for i in range(m):
        tour = metropolis_step(tour=tour,temp=temp)
        config_0 = config_0 if config_0.total_energy() < tour.total_energy() else Tour(copy.deepcopy(tour.cities))
    return config_0

def evolution(temp:float, tour:Tour):
    while temp > oracle_temp:
        zwischenergebnis = evolution_step(temp=temp,tour=tour)
        tour = tour if tour.total_energy() < zwischenergebnis.total_energy() else zwischenergebnis
        temp *= (1-epsilon)
    print(tour)

def init_temp(length:int):    
    configuration = Tour(list(range(length)))
    config_0 = Tour(copy.deepcopy(configuration.cities))
    max_temp = configuration.total_energy()
    for i in range(m):
        configuration.move1() if r.randint(0,100)%2 else configuration.move2()
        max_temp = max(max_temp, configuration.total_energy())
        config_0 = config_0 if config_0.total_energy() < configuration.total_energy() else Tour(copy.deepcopy(configuration.cities))
    return max_temp, config_0

if __name__ == "__main__":
    file = "./ch130.tsp"
    x, y = np.loadtxt(file, delimiter=' ', comments="EOF", skiprows=6, usecols=(1, 2), unpack=True)
    n = len(x)
    for i in range(n):
        city_list.append(City(np.array([x[i],y[i]])))
    tour = Tour([])
    temp, tour = init_temp(n)
    evolution(temp=temp,tour=tour)


 