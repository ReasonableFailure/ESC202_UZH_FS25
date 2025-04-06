import numpy as np
from matplotlib import pyplot as plt

class City:
    def __init__(self, pos:np.ndarray[float]):
        self.x = pos[0]
        self.y = pos[1]

class Tour:
    def __init__(self, cities:list):
        self.length = len(cities)
        self.cities = cities




def metropolis_step():
    pass
