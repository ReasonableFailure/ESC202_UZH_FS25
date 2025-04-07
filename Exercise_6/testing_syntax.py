import random as r
import numpy as np

class City:
    def __init__(self, pos:np.array):
        self.pos = pos
    def __repr__(self):
        return f"{self.pos}"
ind1 = r.randint(0,10)
ind2 = r.randint(0,10)
print(abs(ind2-ind1) != 1)