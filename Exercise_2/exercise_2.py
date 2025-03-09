from heapq import *
from matplotlib.patches import Circle
import numpy as np
from matplotlib import pyplot as plt
# from Exercise_1 import exercise_1 as ex

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

def queue_plotter(pq:prioq, r, period, axis, color = 'red'):
    allneighbours = np.array([p[1].r for p in pq.heap])
    # print(allneighbours)
    ax.scatter(x=allneighbours[:,0], y=allneighbours[:,1], c=color, s=4)
    for y in [0.0, -period[1], period[1]]:
        for x in [0.0, -period[0], period[0]]:
            rOffset = np.array([x, y])
            r2 = r.r - rOffset
            ax.add_patch(Circle(xy=(r2[0], r2[1]), radius= np.sqrt(-pq.key()), edgecolor = 'k', fill=False))

def plot_particles(fig,A:np.ndarray[Particle]):
    fig.scatter([p.r[0] for p in A], [p.r[1] for p in A], color="black", s=4)

def recursive_tree_plotter(axis:plt.axis,root: Cell):
    if(not root.leftChild or not root.rightChild):
        xl = root.lowerBound[0]
        yl = root.lowerBound[1]
        xh = root.upperBound[0]
        yh = root.upperBound[1]
        axis.plot([xl, xh], [yl, yl], color="red")
        axis.plot([xl, xh], [yh, yh], color="red")
        axis.plot([xl, xl], [yl, yh], color="red")
        axis.plot([xh, xh], [yl, yh], color="red")
        return
    if(root.rightChild):
        recursive_tree_plotter(axis = axis,root=root.rightChild)
    if(root.leftChild):
        recursive_tree_plotter(axis=axis,root=root.leftChild)
    return
    

if __name__ == "__main__":
    #initialise data
    A: np.ndarray = np.array([])
    for _ in range(100000):
        p = Particle(np.random.rand(2))
        A = np.append(A, np.array(p))
    root = Cell(rLow=np.array([0.0, 0.0]),rHigh=np.array([1.0, 1.0]),name="root",lo=0,hi=len(A) - 1)
    pq = prioq(300)
    pq2 = prioq(300)
    middle_point = Particle(np.array([0.5,0.5]))
    far_point =Particle(np.array([0.1,0.99]))  

    #operate on data
    tree_builder(root,A,0)
    neighbour_search_periodic(pq=pq,root=root,particles=A,r=middle_point.r,period=np.array([1.0,1.0]))
    neighbour_search_periodic(pq2,root,A,far_point.r,np.array([1.0,1.0]))
    
    #plotting
    fig, ax = plt.subplots()
    plot_particles(ax,A=A)
    recursive_tree_plotter(ax,root=root)
    queue_plotter(pq, middle_point, np.array([1.0,1.0]), ax)
    queue_plotter(pq2, far_point, np.array([1.0,1.0]), ax, color='cyan')
    ax.set_aspect('equal', 'box')
    ax.set_ylim(0,1)
    ax.set_xlim(0,1)
    plt.show()
