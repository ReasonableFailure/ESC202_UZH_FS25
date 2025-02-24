import numpy as np
from matplotlib import pyplot as plt
import random

class Particle:
    def __init__(self, r: np.ndarray):
        self.r = r

class Cell:
    def  __init__(self,rHigh:np.ndarray,rLow:np.ndarray, lo: int, hi:int):
        self.leftChild = None 
        self.rightChild = None
        self.upperBound = rHigh
        self.lowerBound = rLow
        self.index_low = lo
        self.index_high = hi


def partition(A: np.ndarray[Particle],i: int, j:int, v:float,d:bool) -> int :
    if len(A) == 0:
        return None
    interval = A[i : j + 1]
    # Point to the last known particle for which r[d]>v.
    known = 0
    # Keeps track of the current index
    current = 0
    for particle in interval:
        # particle position smaller, swap needed
        if particle.r[d] < v:
            # only swap if not same index, otherwise not change made
            if known < current:
                # search through array until a larger value is found, in order to be swappable
                while interval[known].r[d] < v and known < current:
                    known += 1
                # make the swap
                interval[known], interval[current] = (interval[current], interval[known])
        # advance beyond swapped (now correct) index
        current += 1

    return known + i  #return first index of r[d] > v, accounting for interval starting at i.
    
def tree_builder(root: Cell, A: np.ndarray[Particle], dim: int):
    v = 0.5 * (root.lowerBound[dim] + root.upperBound[dim])
    s = partition(A, root.index_low, root.index_high, v, dim)

    # New cell bounds are set depending on the dimension.
    if dim == 0:
        rLow_Left = root.lowerBound
        rHigh_Left = np.array([v, root.upperBound[1]])
        rLow_Right = np.array([v, root.lowerBound[1]])
        rHigh_Right = root.upperBound
    else:
        rLow_Left = root.lowerBound
        rHigh_Left = np.array([root.upperBound[0], v])
        rLow_Right = np.array([root.lowerBound[0], v])
        rHigh_Right = root.upperBound

    # The left cell is generated if a left partition exists and the branching continued.
    if s > root.index_low:
        left_cell = Cell(rLow_Left, rHigh_Left, root.index_low, s - 1)
        root.leftChild = left_cell
        if len(A[root.index_low:s]) > 8:
            tree_builder(A, left_cell, 1 - dim) # alternate splitting dimensions

    # The right cell is generated if a right partition exists and then recursed into.
    if s <= root.index_high:
        right_cell = Cell(rLow_Right, rHigh_Right, s, root.index_low)
        root.upperCell = right_cell
        if len(A[s:root.index_high + 1]) > 8:
            tree_builder(A, right_cell, 1 - dim) # alternate splitting dimensions


def plot_particles(A:np.ndarray[Particle]):
    plt.scatter([p.r[0] for p in A], [p.r[1] for p in A], color="black")

def recursive_tree_plotter(root: Cell):
    if(root.leftChild):
        recursive_tree_plotter(root=root.leftChild)
    if(root.rightChild):
        recursive_tree_plotter(root=root.rightChild)
    xl = root.lowerBound[0]
    yl = root.lowerBound[1]
    xh = root.upperBound[0]
    yh = root.upperBound[1]
    plt.plot([xl, xh], [yl, yl], color="red")
    plt.plot([xl, xh], [yh, yh], color="red")
    plt.plot([xl, xl], [yl, yh], color="red")
    plt.plot([xh, xh], [yl, yh], color="red")

if __name__ == "__main__":
    A: np.ndarray = np.array([])
    for _ in range(1000):
        A = np.append(A, np.array(Particle(np.random.rand(2))))
    
    root = Cell(rLow=[0.0, 0.0],rHigh=[1.0, 1.0],lo=0,hi=len(A) - 1)
    plot_particles(A=A)
    tree_builder(root,A,0)
    recursive_tree_plotter(root=root)