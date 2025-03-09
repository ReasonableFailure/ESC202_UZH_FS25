from matplotlib import pyplot as plt
import numpy as np
from Exercise_2 import exercise_2 as ex2

class Particle:
    def __init__(self, r: np.ndarray, mass: float):
        self.r = r
        self.m = mass

# Density calculation:
# rho_particle[i] = sum over all N neighbours of their mass times the kernel (dependent on distance and on radius of neighbourhood.)
def tophat_kernel(r:float,h:float) -> float:
    """returns the 2D - tophat kernel based on 2 particles. r is position of point, r-particle_j is jth neighbour, h is radius of neighbourhood ball"""
    return 1/(np.pi * h**2)

def monaghan_kernel(r:float, h:float) -> float:
    """
    param r: length of vector to particle under consideration
    param h: current furthest neighbour radius
    """
    sigma = 40/(7 * np.pi)
    r_over_h = r/h

    if r >= 0 and r_over_h <= 0.5:
        return (sigma/h**2) * (1+6*(r_over_h**3 - r_over_h**2))
    elif r_over_h >= 0.5 and r_over_h <=1:
        return (sigma/h**2) * (2*(1 - r_over_h)**3)
    else:
        return 0.0

def density_calc(particles: np.array[Particle], root:ex2.Cell, Prio_Queue:ex2.prioq, neigh: int ,N:int, kernel:function)->list:
    densities = []
    x = []
    y = []
    tree = ex2.tree_builder(root=root,A=particles,dim=0)

    for particle in particles:
        
        for 

    return densities




if __name__ == "__main__":
    print("hello.")
    # prep data
    No_of_part = 1_000_000
    neighbours = 32
    A: np.ndarray = np.array([])
    for _ in range(No_of_part):
        p = Particle(r = np.random.rand(2), mass=np.random.random())
        A = np.append(A, np.array(p))
    root = ex2.Cell(rLow=np.array([0.0, 0.0]),rHigh=np.array([1.0, 1.0]),name="root",lo=0,hi=len(A) - 1)
    prioqueue = ex2.prioq(neighbours)
    #