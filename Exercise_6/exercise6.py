import numpy as np
from matplotlib import pyplot as plt
import random as r

city_list = []
m = 1000
epsilon = 0.01

class City:
    def __init__(self, pos:np.ndarray[float]):
        self.pos = pos
    def energy_to_neighbour(self,neighbour:City): # type: ignore
        return np.sqrt(np.power(np.subtract(neighbour.pos,self.pos),np.array([2.0,2.0])))

class Tour:
    def __init__(self, cities:list[int]):
        self.length = len(cities)
        self.cities = cities
    def move1(self):
        ind1 = r.randint(0,self.length)
        ind2 = r.randint(0,self.length)
        while ind1 == ind2:
            ind2 = r.randint(0,self.length)
        self.cities[ind2], self.cities[ind1] = self.cities[ind1], self.cities[ind2]
    def move2(self):
        ind1 = r.randint(0,self.length)
        ind2 = r.randint(0,self.length)
        while ind1 == ind2:
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
        return ret

def metropolis_step(tour:Tour):
    """TODO: Calculate Energy of Configuration. Select new configuration. By either: Move 1: Swapping 2 non-adjacent indices or Move 2: selecting a subsequence and inverting it. Calculate Energy of Configuration. Do balancing
    """
    E_old = tour.total_energy()
    
