from matplotlib import pyplot as plt
from matplotlib import animation as an
import numpy as np
from heapq import *
# import matplotlib.animation as animation

class Particle:
    def __init__(self, r: np.ndarray, mass: float, vel: np.ndarray, U:float, id:int):
        self.m = mass
        self.r = r
        self.velocity = vel
        self.internal_energy = U
        self.a = np.array([0.0,0.0])
        self.density = 0.0
        self.h = 0.0
        self.velocity_pred = np.array([0.0,0.0])
        # self.pressure = None # might not be necessary, but nice for visualisation, eventually.
        self.name = id
        self.soundspeed = 0.0
        self.U_dot = 0.0
        self.U_pred = 0.0
        #lt operator
    def __lt__(self, other):
        np.sum(self.r**2) < np.sum(other.r**2)
    def __repr__(self):
        # return f"{self.name} at r = {self.r}, with v = {self.velocity}, a = {self.a}, soundspeed = {self.soundspeed}\nWeight = {self.m}, Density = {self.density}, Internal Energy = {self.internal_energy}\nNeighbourshood size = {self.h}"
        return f"{self.name} Density = {self.density}"

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
        return f"[{self.index_low},{self.index_high}]"

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
    def cleanup(self,k: int):
        self.__init__(k)




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

def neighbour_search(pq:prioq, root:Cell, particles, r, rOffset, id): #this is ball_walk
    cnt = 0
    # print(f"\n\n\nStart new iteration. cnt = {cnt}")
    # print(type(particles))
    if not root:
        return 0
    if isLeaf(root):
        for i in range(root.index_low,root.index_high+1):
            part = particles[i]
            if(part.name == id):
                # with open("0_density_particles.txt","a") as file:
                #     file.write(f"Self-interaction {id}\n")
                #     file.write(f"Location: {r}\n")
                #     file.write(f"Location of particle: {part.r}\n")
                #     file.write(f"delta = ( {part.r} + {rOffset}) - {r} = {( part.r + rOffset) - r}\n")
                continue
            delta = ( part.r + rOffset) - r
            if(delta[0] == 0.0 and delta[1] == 0.0):
                with open("0_density_particles.txt","a") as file:
                    file.write(f"0 Density particle {id}\n")
                    file.write(f"Location: {r}\n")
                    file.write(f"Location of particle: {part.r}\n")
                    file.write(f"delta = ( {part.r} + {rOffset}) - {r} = {( part.r + rOffset) - r}\n")
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
                c = neighbour_search(pq=pq,root=root.rightChild,particles=particles,r=r,rOffset=rOffset, id=id)
                cnt += c
        if root.leftChild:
            if -celldist2(root.leftChild,r-rOffset) > pq.key():
                c = neighbour_search(pq=pq,root=root.leftChild,particles=particles,r=r,rOffset=rOffset, id=id)
                cnt+=c
    return cnt

def neighbour_search_periodic(pq, root, particles, r, period, particle_id):
    # walk the closest image first (at offset=[0, 0])
    for y in [0.0, -period[1], period[1]]:
        for x in [0.0, -period[0], period[0]]:
            rOffset = np.array([x, y])
            neighbour_search(pq, root, particles, r, rOffset, particle_id)

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
    result = sigma/h**2
    if r >= 0 and r_over_h < 0.5:
        result *= (1+6*(r_over_h**3 - r_over_h**2))
    elif r_over_h >= 0.5 and r_over_h <=1:
        result *= (2*(1 - r_over_h)**3)
    else:
        result = 0.0
    # print(result)
    return result
    
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
    # print(particle_i.density)
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
    # print("Calculation of a based on interaction:")
    # print(f"Particle i = {particle_i}")
    # print(f"Particle j = {particle_j}")
    #formula from lecturer's notes i.e. loop body. separated out for clarity, may be reincorporated for less space intensive code.
    particle_i.a -= particle_j.m*(calc_pressure_term(gamma=gamma,particle_i=particle_i) + calc_pressure_term(gamma=gamma,particle_i=particle_j) + calc_viscosity_term(particle_i=particle_i,particle_j=particle_j))*symmetrify_kernel(kernel=derivative_monaghan,r=r_ij,h_i=particle_i.h,h_j=particle_j.h)
    # particle_i.a -= 0.5 * particle_j.m * derivative_monaghan(r=r_ij,h=particle_j.h)* (
    #     calc_pressure_term(gamma=gamma,particle_i=particle_i)+
    #     calc_pressure_term(gamma=gamma,particle_i=particle_j)+
    #     calc_viscosity_term(particle_i=particle_i,particle_j=particle_j))        
    # particle_j.a -= 0.5 * particle_i.m * derivative_monaghan(r=r_ij,h=particle_i.h) * (
    #     calc_pressure_term(gamma=gamma,particle_i=particle_i)+
    #     calc_pressure_term(gamma=gamma,particle_i=particle_j)+
    #     calc_viscosity_term(particle_i=particle_i,particle_j=particle_j)
    

def Du_over_Dt(particle_i:Particle, particle_j:Particle, gamma:float):
    #auxiliary computation
    r_ij = np.linalg.norm(particle_j.r-particle_i.r)
    v_ij = np.linalg.norm(particle_i.velocity - particle_j.velocity)
    # print("Calculation of Udot based on interaction of")
    # print(f"Particle i = {particle_i}")
    # print(f"Particle j = {particle_j}")
    #formula from lecturer's notes i.e. loop body. separated out for clarity, may be reincorporated for less space intensive code.
    particle_i.U_dot += calc_pressure_term(gamma,particle_i)*particle_j.m*v_ij*symmetrify_kernel(derivative_monaghan,r_ij,particle_i.h,particle_j.h) + 0.5*particle_j.m*calc_viscosity_term(particle_i,particle_j)*v_ij*symmetrify_kernel(kernel=derivative_monaghan,r=r_ij,h_i=particle_i.h,h_j=particle_j.h)

def calc_density_i(particle_i :Particle, neigh:int, prio_Queue:prioq, particles:np.array, kernel) ->None:
    particle_i.h = np.sqrt(-prio_Queue.key())
    # print(f"Particle i = {particle_i}")
    for j in range(neigh):
        distance, particle_j, dr = prio_Queue.heap[j]
        # print(f"Particle j = {particle_j}")
        mass_j = particle_j.m
        r_i_r_j = np.sqrt(-distance)
        rho_j = mass_j * kernel(r=r_i_r_j,h=particle_i.h)
        particle_i.density += rho_j
    # print(f"Particle i = {particle_i}")
    

def calc_forces(particles:np.array,neigh:int) -> None:
    # for periodic boundary conditions
    gamma = 5/3 #this holds for any gas with degree of freedom f = 3. gamma = (f+2)/f. Number taken from lecturer's notes
    period = np.array([1.0,1.0])
    root = Cell(rLow=np.array([0.0, 0.0]),rHigh=np.array([1.0, 1.0]),name="root",lo=0,hi=len(particles) - 1)
    root = tree_builder(root=root,A=particles,dim=0)
    prio_Queue = prioq(neigh)
    for particle in particles:
        prio_Queue.cleanup(neigh)
        neighbour_search_periodic(pq=prio_Queue,root=root,particles=particles,r=particle.r,period=period, particle_id = particle.name)
        calc_density_i(particle_i=particle,neigh=neigh,particles=particles,prio_Queue=prio_Queue,kernel=monaghan_kernel)
        # if(particle.density == 0.0):
        #     with open("0_density_particles.txt","w") as file:
        #         file.write(f"{particle}")
        #         file.write(f"{prio_Queue.heap}")
        #         file.write("\n\n\n\n\n\n\n\n")                
        calc_soundspeed(particle=particle,gamma=gamma)
    for particle in particles:
        neighbour_forces_i(patricle_i=particle,prio_Queue=prio_Queue,neigh=neigh,gamma=gamma)
    
def sph_init(part_num:int, neigh:int)->np.ndarray:
    #init section
    #populate particles    
    A = []
    for _ in range(No_of_part):
        p = Particle(r=np.random.rand(2),mass=1.0,vel=np.array([0.0,0.0]),U=10.0, id=_)
        if _ == 225:
            p.internal_energy = 1000.0
        A.append(p)
    #initialise upred, vpred, density, acceleration,...   
    drift_one(particles=A,delta_t=0)
    calc_forces(particles=A,neigh=neigh)
    return A

def sph_leapfrog( delta_t : float,A:np.array,neigh:int) -> np.ndarray:

    drift_one(particles=A,delta_t=0.5*delta_t)
    calc_forces(particles=A,neigh=neigh)
    kick(particles=A,delta_t=delta_t)
    drift_two(particles=A,delta_t=0.5*delta_t)
    return A


delta_t = 0.0003
No_of_part = int(input("Number of points\n"))
neighbours = int(input("neigbourhood points\n"))
repetition = int(input("How many iterations?\n"))
particle_trace = np.ndarray(shape=(No_of_part,repetition+1),dtype=Particle)
particle_trace[:,0]=np.array(sph_init(part_num=No_of_part,neigh=neighbours))
fig, axs = plt.subplots(nrows=1,ncols=2)
artists=[]
x = []
y = []
rho = []
U = []
for particle in particle_trace[:,0]:        
        x.append(particle.r[0])
        y.append(particle.r[1])
        rho.append(particle.density)
        U.append(particle.internal_energy)
axs[0].set(xlim = [0,1],ylim=[0,1],xlabel="Density Plot")
axs[1].set(xlim = [0,1],ylim=[0,1],xlabel="Internal Energy Plot")
density_plot = axs[0].scatter(x=x,y=y,c=rho,cmap="plasma",s=2)
energy_plot = axs[1].scatter(x=x,y=y,c=U,cmap="seismic",s=2)
artists.append([density_plot,energy_plot])
for rep in range(repetition):
    print(rep)
    particle_trace[:,rep+1] = sph_leapfrog(delta_t=delta_t,A=particle_trace[:,rep],neigh=neighbours) 
    x = []
    y= []
    rho = []
    U = []   
    for particle in particle_trace[:,rep+1]:        
        x.append(particle.r[0])
        y.append(particle.r[1])
        rho.append(particle.density)
        U.append(particle.internal_energy)
    density_plot = axs[0].scatter(x=x,y=y,c=rho,cmap="plasma",s=2)
    energy_plot = axs[1].scatter(x=x,y=y,c=U,cmap="seismic",s=2)
    artists.append([density_plot,energy_plot])

ani=an.ArtistAnimation(fig=fig,artists=artists,interval=300)
ani.save(writer="ffmpeg", filename="/home/faye/UZH/6Sem/ESC202/animation-ex4.mp4")




    



