import numpy as np
from matplotlib import pyplot as plt
import random

class Particle:
    def __init__(self, r: np.ndarray):
        self.r = r

class Cell:
    def  __init__(self,rHigh:np.array,rLow:np.array, name:str, lo, hi):
        self.leftChild = None 
        self.rightChild = None
        self.upperBound = rHigh
        self.lowerBound = rLow
        self.index_low = lo
        self.index_high = hi
        self.name = name
    def __repr__(self):
        return f"{id(self)}"


def partition(A: np.array,i: int, j:int, v:float,d:bool) -> int :
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
    
def tree_builder(root: Cell, A: np.array, dim: int):
    v = 0.5 * (root.lowerBound[dim] + root.upperBound[dim])
    pivot_index = partition(A, root.index_low, root.index_high, v, dim)
    # print(f"pivot_index = {pivot_index}")
    #initialise left child cell
    #split in x direction
    if dim == 0:
        ll_l = root.lowerBound
        ur_l = np.array([v,root.upperBound[1]])
    #split in y direction
    else:
        ll_l = root.lowerBound
        ur_l = np.array([root.upperBound[0],v])

    
    #initialise right child
    #split along x axis
    if dim == 0:
        ll_r = np.array([v,root.lowerBound[1]])
        ur_r = root.upperBound
    else:
        ll_r = np.array([root.lowerBound[0],v])
        ur_r = root.upperBound
    
    leftChild = Cell(ur_l,ll_l,root.name+"l",root.index_low,pivot_index)
    worth_it_left = pivot_index - root.index_low
    # print(f"worth_it_left = {worth_it_left}")
    
    
    rightChild = Cell(ur_r,ll_r,root.name+"r",pivot_index+1,root.index_high)
    worth_it_right = root.index_high - pivot_index -1
    # print(f"worth_it_right = {worth_it_right}")
    if root.index_high-root.index_low - 1 > 8:
        root.leftChild = leftChild
        root.rightChild = rightChild
        tree_builder(leftChild,A,(1-dim))
        tree_builder(rightChild,A,(1-dim))
    return root

    # # New cell bounds are set depending on the dimension.
    # if dim == 0:
    #     rLow_Left = root.lowerBound
    #     rHigh_Left = np.array([v, root.upperBound[1]])
    #     rLow_Right = np.array([v, root.lowerBound[1]])
    #     rHigh_Right = root.upperBound
    # else:
    #     rLow_Left = root.lowerBound
    #     rHigh_Left = np.array([root.upperBound[0], v])
    #     rLow_Right = np.array([root.lowerBound[0], v])
    #     rHigh_Right = root.upperBound

    # # print(f"rLowLeft = {rLow_Left}, rHighLeft = {rHigh_Left}, rLowRight = {rLow_Right}, rHighRight = {rHigh_Right}")
    # # The left cell is generated if a left partition exists and the branching continued.
    # if s > root.index_low:
    #     left_cell = Cell(rLow_Left, rHigh_Left, root.index_low, s - 1)
    #     root.leftChild = left_cell
    #     if len(A[root.index_low:s]) > 8:
    #         tree_builder(root=left_cell,A=A,dim=(1-dim)) # alternate splitting dimensions

    # # The right cell is generated if a right partition exists and then recursed into.
    # if s <= root.index_high:
    #     right_cell = Cell(rLow_Right, rHigh_Right, s, root.index_high)
    #     root.rightChild = right_cell
    #     if len(A[s:root.index_high + 1]) > 8:
    #         # print("here")
    #         tree_builder(root=right_cell,A=A,dim=(1-dim)) # alternate splitting dimensions ICH HAB DIH LIE DU MURPEL HIHIHIHIHIHIHIIH


def plot_particles(A:np.ndarray[Particle]):
    plt.scatter([p.r[0] for p in A], [p.r[1] for p in A], color="black", s=3)

def recursive_tree_plotter(root: Cell):
    if(not root.leftChild or not root.rightChild):
        xl = root.lowerBound[0]
        yl = root.lowerBound[1]
        xh = root.upperBound[0]
        yh = root.upperBound[1]
        plt.plot([xl, xh], [yl, yl], color="red")
        plt.plot([xl, xh], [yh, yh], color="red")
        plt.plot([xl, xl], [yl, yh], color="red")
        plt.plot([xh, xh], [yl, yh], color="red")
        return
    if(root.rightChild):
        recursive_tree_plotter(root=root.rightChild)
    if(root.leftChild):
        recursive_tree_plotter(root=root.leftChild)
    return


def treewalker(root:Cell):
    s = f"My name is {root.name}. My indices are [{root.index_low},{root.index_high}]."#\nMy children span:"
    if(not root.leftChild and not root.rightChild):
        s = s + " I am leaf."
        print(s)
        return
    
    if(root.leftChild):
        #s = s + f" left:[{root.leftChild.index_low},{root.leftChild.index_high}], "
        treewalker(root.leftChild)    
    if(root.rightChild):
        #s = s + f" right:[{root.rightChild.index_low},{root.rightChild.index_high}]\n"
        treewalker(root.rightChild)
    print(s)
    return


if __name__ == "__main__":
    A: np.ndarray = np.array([])
    for _ in range(1000):
        p =Particle(np.random.rand(2))
        A = np.append(A, np.array(p))
    
    root = Cell(rLow=np.array([0.0, 0.0]),rHigh=np.array([1.0, 1.0]),name="root",lo=0,hi=len(A) - 1)
    tree = tree_builder(root,A,0)
    treewalker(tree)
    plot_particles(A)
    recursive_tree_plotter(root=root)
    plt.show()
    