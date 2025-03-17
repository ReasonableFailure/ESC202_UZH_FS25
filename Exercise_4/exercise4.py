from matplotlib import pyplot as plt
import numpy as np
import matplotlib.animation as animation
from heapq import *

class Particle:
    def __init__(self, r: np.ndarray, mass: float, vel: np.ndarray, U:float):
        self.m = mass
        self.r = r
        self.velocity = vel
        self.internal_energy = U
        self.a = np.array([0.0,0.0])
        self.density = np.nan
        self.h = np.nan
        self.velocity_pred = np.array([0.0,0.0])
        # self.pressure = None # might not be necessary, but nice for visualisation, eventually.
        self.soundspeed = np.nan
        self.U_dot = np.nan
        self.U_pred = np.nan
        #lt operator
        def __lt__(self, other:Particle):
            np.sum(self.r**2) < np.sum(other.r**2)

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
    # print(type(particles))
    if not root:
        return 0
    if isLeaf(root):
        for i in range(root.index_low,root.index_high+1):
            part = particles[i]
            delta = ( part.r + rOffset) - r

            distance = -delta.dot(delta) - 1e-8 # minus necessary because of min heap.
            # print(f"Particle[{i}] at {part.r} with offset {rOffset}.\nCentre is {r}\nDelta = {delta}, Distance = {distance}")
            # print(type(distance))
            #if the distance here is exactly -0.0 it triggers the prioq library function to compare the 2nd element as the key. fix ?
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

# def tophat_kernel(r:float,h:float) -> float:
#     """returns the 2D - tophat kernel based on 2 particles. r is position of point, r-particle_j is jth neighbour, h is radius of neighbourhood ball"""
#     return 1/(np.pi * h**2)

def monaghan_kernel(r:float, h:float) -> float:
    """
    param r: length of vector to particle under consideration
    param h: current furthest neighbour radius
    """
    sigma = 40/(7 * np.pi)
    r_over_h = r/h

    if r >= 0 and r_over_h < 0.5:
        return (sigma/h**2) * (1+6*(r_over_h**3 - r_over_h**2))
    elif r_over_h >= 0.5 and r_over_h <=1:
        return (sigma/h**2) * (2*(1 - r_over_h)**3)
    else:
        return 0.0
    
def derivative_monaghan(r:float,h:float) -> float:
    sigma = 40/(7 * np.pi)
    r_over_h = r/h

    result = 6*sigma / (h**3)

    if r >= 0  and r_over_h < 0.5:
        result * (3*r-r_over_h**2 - 2*r_over_h)
    elif r_over_h >= 0.5 and r_over_h <=1:
        result * (-1)*((1-r_over_h)**2)
    else:
        result = 0.0
    return result

#deprecated. will remain for reference and visualisation
"""def density_calc(particles, root:Cell, neigh: int ,N:int, kernel)->tuple:
    densities = []
    x = []
    y = []
    tree = tree_builder(root=root,A=particles,dim=0)
    period = np.array([1.0,1.0])

    for particle in particles:
        prio_Queue = prioq(neigh)
        neighbour_search_periodic(pq=prio_Queue,root=root,particles=particles,r=particle.r,period=period)
        h = np.sqrt(-prio_Queue.key())
        particle.h = h
        rho_i = 0
        for j in range(neigh):
            particle_j = prio_Queue.heap[j]
            mass_j = particle_j[1].m
            r_i_r_j = np.sqrt(-particle_j[0])
            rho_j = mass_j * kernel(r_i_r_j,h)
            rho_i += rho_j
        #particle.density = rho_i
        densities.append(rho_i)
        x.append(particle.r[0])
        y.append(particle.r[1])
    return (x,y,densities)"""

def symmetrify_kernel(kernel, r:float, h_i:float, h_j:float):
    return 0.5*(kernel(r=r,h=h_i)+kernel(r=r,h=h_j))

def drift_one(particles:np.array, delta_t:float):
    """This function is potentially dangerous because it only works via side effects. Drift step one, calculates predicted speed and internal energy"""
    for particle_i in particles:
        particle_i.r += particle_i.velocity*delta_t
        particle_i.velocity_pred = particle_i.velocity + particle_i.a * delta_t
        particle_i.U_pred = particle_i.internal_energy + particle_i.U_dot * delta_t

def drift_two(particles:np.array, delta_t:float):
    """This function is potentially dangerous because it only works via side effects. Drift step two. Updates position"""
    for particle_i in particles:
        particle_i.r += particle_i.velocity*delta_t

def kick(particles:np.array, delta_t:float):
    """This function is potentially dangerous because it only works via side effects. Updates acceleration and U dot."""
    for particle_i in particles:
        particle_i.velocity += particle_i.a*delta_t
        particle_i.internal_energy += particle_i.U_dot*delta_t

def calc_pressure_term(gamma:float, particle_i:Particle): 
    #particle_i.pressure = (gamma-1)*particle_i.density*particle_i.internal_energy
    return (particle_i.soundspeed)**2 / (gamma*particle_i.density)

def calc_soundspeed(gamma:float, particle:Particle):
    """This function is potentially dangerous because it only works via side effects. Calculate sound speed given parameter gamma about gas"""
    particle.soundspeed = np.sqrt(gamma * (gamma - 1) * particle.U_pred)

def calc_viscosity_term(particle_i:Particle, particle_j:Particle) -> float:
    nu = 1e-6 #correction factor. prevents division by 0. Could be machine epsilon.
    c_ij = 0.5 * (particle_i.soundspeed + particle_j.soundspeed)
    rho_ij = 0.5 * (particle_i.density + particle_j.density)
    h_ij = 0.5 * (particle_i.h + particle_j.h)
    v_ij = particle_i.velocity - particle_j.velocity
    r_ij = particle_j.r - particle_i.r
    abs_r_ij = np.linalg.norm(r_ij)
    mu_ij = h_ij * (np.dot(v_ij,r_ij)) / (abs_r_ij + nu**2)

    alpha = 1 #from lecturer's notes. this constant influences how the fluid behaves, so will have to be filled with experimentally measured parameters. could also be 0.5
    beta = 2*alpha

    if np.dot(v_ij,r_ij) > 0: #convergence
        return (beta*mu_ij**2 - alpha*c_ij*mu_ij)/rho_ij
    else:
        return 0

def neighbour_forces_i(patricle_i:Particle, prio_Queue:prioq, neigh:int, gamma:float):
    for j in range(neigh):
        dist, particle_j, dr = prio_Queue.heap[j]
        Dvi_over_Dt(particle_i=patricle_i,particle_j=particle_j,gamma=gamma)
        Du_over_Dt(particle_i=patricle_i,particle_j=particle_j,gamma=gamma)

def Dvi_over_Dt(particle_i:Particle, particle_j:Particle, gamma:float):
    #auxiliary computation
    r_ij = np.linalg.norm(particle_j.r - particle_i.r)
    #formula from lecturer's notes i.e. loop body. separated out for clarity, may be reincorporated for less space intensive code.
    particle_i.a -= particle_j.m*(calc_pressure_term(gamma=gamma,particle_i=particle_i) + calc_pressure_term(gamma=gamma,particle_i=particle_j) + calc_viscosity_term(particle_i=particle_i,particle_j=particle_j))*symmetrify_kernel(kernel=derivative_monaghan,r=r_ij,h_i=particle_i.h,h_j=particle_j.h)

def Du_over_Dt(particle_i:Particle, particle_j:Particle, gamma:float):
    #auxiliary computation
    r_ij = np.linalg.norm(particle_j.r-particle_i.r)
    v_ij = np.linalg.norm(particle_i.velocity - particle_j.velocity)
    #formula from lecturer's notes i.e. loop body. separated out for clarity, may be reincorporated for less space intensive code.
    particle_i.U_dot += calc_pressure_term(gamma,particle_i)*particle_j.m*v_ij*symmetrify_kernel(derivative_monaghan,r_ij,particle_i.h,particle_j.h) + 0.5*particle_j.m*calc_viscosity_term(particle_i,particle_j)*v_ij*symmetrify_kernel(kernel=derivative_monaghan,r=r_ij,h_i=particle_i.h,h_j=particle_j.h)

def calc_density_i(particle_i :Particle, neigh:int, prio_Queue:prioq, particles:np.array, kernel) ->None:
    particle_i.h = np.sqrt(-prio_Queue.key())
    particle_i.density = 0
    for j in range(neigh):
        particle_j = prio_Queue.heap[j]
        mass_j = particle_j[1].m
        r_i_r_j = np.sqrt(-particle_j[0])
        rho_j = mass_j * kernel(r=r_i_r_j,h=particle_i.h)
        particle_i.density += rho_j

def calc_forces(particles:np.array,neigh:int) -> None:
    # for periodic boundary conditions
    gamma = 5/3 #this holds for any gas with degree of freedom f = 3. gamma = (f+2)/f. Number taken from lecturer's notes
    period = np.array([1.0,1.0])
    root = Cell(rLow=np.array([0.0, 0.0]),rHigh=np.array([1.0, 1.0]),name="root",lo=0,hi=len(particles) - 1)
    root = tree_builder(root=root,A=particles,dim=0)
    prio_Queue = prioq(neigh)
    for particle in particles:
        neighbour_search_periodic(pq=prio_Queue,root=root,particles=particles,r=particle.r,period=period)
        calc_density_i(particle_i=particle,neigh=neigh,particles=particles,prio_Queue=prio_Queue,kernel=monaghan_kernel)
        calc_soundspeed(particle=particle,gamma=gamma)
        neighbour_forces_i(patricle_i=particle,prio_Queue=prio_Queue,neigh=neigh,gamma=gamma)
    

def sph_leapfrog(iterations : int, delta_t : float,part_num:int,neigh:int) -> None:
    #init section
    #populate particles    
    A: np.ndarray = np.array([])
    for _ in range(No_of_part):
        p = Particle(r=np.random.rand(2),mass=1.0,vel=np.array([0,0]),U=10.0)
        A = np.append(A, np.array(p)) 
    #initialise upred, vpred, density, acceleration,...   
    drift_one(particles=A,delta_t=0)
    calc_forces(particles=A,neigh=neigh)
    for step in range(iterations):
        drift_one(particles=A,delta_t=0.5*delta_t)
        calc_forces(particles=A,neigh=neigh)
        kick(particles=A,delta_t=delta_t)
        drift_two(particles=A,delta_t=0.5*delta_t)




if __name__ == "__main__":
    delta_t = 0.0003 # from exercise description
    neighbours = int(input("neigbourhood points\n"))
    No_of_part = int(input("Number of points\n"))
    repetition = int(input("How many iterations?\n"))
    sph_leapfrog(iterations=repetition,delta_t=delta_t,part_num=No_of_part,neigh=neighbours)

