from heapq import *
import numpy as np
from Exercise_1 import exercise_1 as ex

# Priority queue
class prioq:
    def __init__(self, k):
        self.heap = []
        sentinel = (-np.inf, None, np.array([0.0,0.0]))
        for i in range(k):
            heappush(self.heap, sentinel)

    def replace(self, dist2:float, particle:ex.Particle, dr):
        insert_node = (dist2,particle,dr)
        former_closest = heapreplace(self.heap,insert_node)
        return former_closest

    def key(self):
        # .... define key here
        minimum_node = self.heap[0]
        return minimum_node[0] #Distance in the tuple.
def celldist2(self, r):
    """Calculates the squared minimum distance between a particle
    position and this node."""
    d1 = r - self.upperBound
    d2 = self.lowerBound - r
    d1 = np.maximum(d1, d2)
    d1 = np.maximum(d1, np.zeros_like(d1))
    return d1.dot(d1)
    
k = 32
pq = prioq(k)

def isLeaf(cell:ex.Cell):
    return (cell.leftChild==None) and (cell.rightChild==None)
        

def neighbour_search(pq:prioq, root:ex.Cell, particles:ex.Particle, r, rOffset): #this is ball_walk
    cnt = 0
    if isLeaf(root):
        for part in particles[root.index_low:root.index_high+1]:
            delta = part.r + rOffset - r
            distance = - delta.dot(delta)
            if distance > pq.key():
                formerly = pq.replace(distance,part,delta)
    else:
        if -celldist2(root.rightChild,r-rOffset)>pq.key():
            neighbour_search(pq=pq,root=root.rightChild,particles=particles,r=r,rOffset=rOffset)
        if -celldist2(root.leftChild,r)>pq.key():
            neighbour_search(pq=pq,root=root.leftChild,particles=particles,r=r,rOffset=rOffset)



def neighbour_search_periodic(pq, root, particles, r, period):
    # walk the closest image first (at offset=[0, 0])
    for y in [0.0, -period[1], period[1]]:
        for x in [0.0, -period[0], period[0]]:
            rOffset = np.array([x, y])
            neighbour_search(pq, root, particles, r, rOffset)