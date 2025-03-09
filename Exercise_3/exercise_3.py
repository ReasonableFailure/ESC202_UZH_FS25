from matplotlib import pyplot as plt
import numpy as np
from heapq import *
# from Exercise_2 import exercise_2 as ex2

class Particle:
    def __init__(self, r: np.ndarray, mass: float):
        self.r = r
        self.m = mass

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

# Priority queue
class prioq:
    def __init__(self, k):
        self.heap = []
        sentinel = (-np.inf, None, np.array([0.0,0.0]))
        for i in range(k):
            heappush(self.heap, sentinel)

    def replace(self, dist2:float, particle:Particle, dr):
        insert_node = (dist2,particle,dr)
        former_closest = heapreplace(self.heap,insert_node)
        return former_closest

    def key(self):
        # .... define key here
        minimum_node = self.heap[0]
        return minimum_node[0] #Distance in the tuple.

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
    if root.index_high-root.index_low > 8:
        root.leftChild = leftChild
        root.rightChild = rightChild
        tree_builder(leftChild,A,(1-dim))
        tree_builder(rightChild,A,(1-dim))
    return root

def celldist2(self:Cell, r:np.array):
    """Calculates the squared minimum distance between a particle
    position and this node."""
    if not self:
        return -np.inf 
    d1 = r - self.upperBound
    d2 = self.lowerBound - r
    d1 = np.maximum(d1, d2)
    d1 = np.maximum(d1, np.zeros_like(d1))
    return d1.dot(d1)

def isLeaf(cell:Cell):
    if cell:
        return not cell.leftChild and not cell.rightChild
    else:
        return False

def neighbour_search(pq:prioq, root:Cell, particles, r, rOffset): #this is ball_walk
    cnt = 0
    # print(f"\n\n\nStart new iteration. cnt = {cnt}")
    if not root:
        return 0
    if isLeaf(root):
        for i in range(root.index_low,root.index_high+1):
            part = particles[i]
            delta = ( part.r + rOffset) - r

            distance = -delta.dot(delta) # minus necessary because of min heap.
            # print(f"Particle[{i}] at {part.r} with offset {rOffset}.\nCentre is {r}\nDelta = {delta}, Distance = {distance}")
            if distance > pq.key():
                # print(f"Before: {[(l[0],l[2]) for l in pq.heap]}\n")
                formerly = pq.replace(distance,part,delta)
                # print(f"After: {[(l[0],l[2]) for l in pq.heap]}\n")
                cnt += 1
    else:  
        if root.rightChild:
            if -celldist2(root.rightChild,r-rOffset) > pq.key():
                c = neighbour_search(pq=pq,root=root.rightChild,particles=particles,r=r,rOffset=rOffset)
                cnt += c
        if root.leftChild:
            if -celldist2(root.leftChild,r-rOffset) > pq.key():
                c = neighbour_search(pq=pq,root=root.leftChild,particles=particles,r=r,rOffset=rOffset)
                cnt+=c
    return cnt

def neighbour_search_periodic(pq, root, particles, r, period):
    # walk the closest image first (at offset=[0, 0])
    for y in [0.0, -period[1], period[1]]:
        for x in [0.0, -period[0], period[0]]:
            rOffset = np.array([x, y])
            neighbour_search(pq, root, particles, r, rOffset)

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

def density_calc(particles, root:Cell, neigh: int ,N:int, kernel)->tuple:
    densities = []
    x = []
    y = []
    tree = tree_builder(root=root,A=particles,dim=0)
    period = np.array([1.0,1.0])

    for particle in particles:
        prio_Queue = prioq(neigh)
        neighbour_search_periodic(pq=prio_Queue,root=root,particles=A,r=particle.r,period=period)
        h = np.sqrt(-prio_Queue.key())
        rho_i = 0
        for j in range(neigh):
            particle_j = prio_Queue.heap[j]
            mass_j = particle_j[1].m
            r_i_r_j = np.sqrt(-particle_j[0])
            rho_j = mass_j * kernel(r_i_r_j,h)
            rho_i += rho_j
        densities.append(rho_i)
        x.append(particle.r[0])
        y.append(particle.r[1])
    return (x,y,densities)


if __name__ == "__main__":
    print("hello.")
    # prep data
    No_of_part = 10_000
    neighbours = 300
    A: np.ndarray = np.array([])
    for _ in range(No_of_part):
        p = Particle(r = np.random.rand(2), mass=1.0)
        A = np.append(A, np.array(p))
    root = Cell(rLow=np.array([0.0, 0.0]),rHigh=np.array([1.0, 1.0]),name="root",lo=0,hi=len(A) - 1)
    print("Calculating densities has started (get a refreshment while we wait)")
    x_th,y_th,densities_th = density_calc(particles=A,root=root,neigh=neighbours,N=No_of_part,kernel=tophat_kernel)
    x_m,y_m,densities_m = density_calc(particles=A,root=root,neigh=neighbours,N=No_of_part,kernel=monaghan_kernel)
    print("We're plotting!")
    fig,ax = plt.subplots(nrows=1,ncols=2)
    ax[0].scatter(x=x_th,y=y_th,s=40,c=densities_th,cmap="plasma")
    ax[1].scatter(x=x_m,y=y_m,s=40,c=densities_m,cmap="plasma")
    ax[0].set_title("Top Hat Kernel")
    ax[1].set_title("Monaghan Kernel") 
    plt.savefig(f"Density-for-{No_of_part}-2-Kernels.png")   
    plt.show()
    
    #