import numpy as np
from matplotlib import pyplot as plt
import random as r
from numba import jit

city_list = []
m = 1000
epsilon = 0.01

class City:
    def __init__(self, pos:np.array):
        self.pos = pos
    def __repr__(self):
        return f"{self.pos}"
    def energy_to_neighbour(self,neighbour): # type: ignore
        return np.sqrt(np.power(np.subtract(neighbour.pos,self.pos),np.array([2.0,2.0])))

class Tour:
    def __init__(self, cities:list[int]):
        self.length = len(cities)
        self.cities = cities
    def __repr__(self):
        return f"{self.cities}"
    @jit
    def move1(self):
        ind1 = r.randint(0,self.length)
        ind2 = r.randint(0,self.length)
        while ind1 == ind2 and abs(ind2-ind1) != 1:
            ind2 = r.randint(0,self.length)
        self.cities[ind2], self.cities[ind1] = self.cities[ind1], self.cities[ind2]
    @jit
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
    @jit
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
    tour.move1() if r.randint()%2 else tour.move2()
    E_new = tour.total_energy()
    delta = E_new - E_old
    factor = np.exp(-delta/temp)
    r = np.random.uniform(-3000.0,3000.0)
    if delta < 0 or (r > 0  and r < factor):
        return tour
    else:
        return tour_old

@jit
def evolution_step(temp:float, tour:Tour):
    config_0 = Tour(tour.cities)
    for i in range(m):
        tour = metropolis_step(tour=tour,temp=temp)
        config_0 = config_0 if config_0.total_energy() < tour.total_energy() else Tour(tour.cities)
    return config_0



@jit
def init_temp(length:int):    
    configuration = Tour(list(range(length)))
    config_0 = Tour(configuration.cities)
    max_temp = configuration.total_energy()
    for i in range(m):
        configuration.move1() if r.randint()%2 else configuration.move2()
        max_temp = max(max_temp, configuration.total_energy())
        config_0 = config_0 if config_0.total_energy() < configuration.total_energy() else Tour(configuration.cities)
    return max_temp, config_0

if __name__ == "__main__":
    n = 5 #number of cities
    for i in range(n):
        city_list.append(City(np.random.rand(2)))
    temp,tour = init_temp(n)

 