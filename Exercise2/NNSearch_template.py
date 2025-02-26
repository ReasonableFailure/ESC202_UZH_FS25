from heapq import *

# Priority queue
class prioq:
    def __init__(self, k):
        self.heap = []
        sentinel = (-np.inf, None, np.array([0.0,0.0]))
        for i in range(k):
            heappush(self.heap, sentinel)

    def replace(self, dist2, particle, dr):
        .... use heapreplace here

    def key(self):
        .... define key here

def celldist2(self, r):
    """Calculates the squared minimum distance between a particle
    position and this node."""
    d1 = r - self.rHigh
    d2 = self.rLow - r
    d1 = np.maximum(d1, d2)
    d1 = np.maximum(d1, np.zeros_like(d1))
    return d1.dot(d1)

pq = prioq(k)

def neighbor_search_periodic(pq, root, particles, r, period):
    # walk the closest image first (at offset=[0, 0])
    for y in [0.0, -period[1], period[1]]:
        for x in [0.0, -period[0], period[0]]:
            rOffset = np.array([x, y])
            neighbor_search(pq, root, particles, r, rOffset)