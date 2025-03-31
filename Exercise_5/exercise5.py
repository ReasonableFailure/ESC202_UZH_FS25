import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation as ani
from numba import jit
from random import randint, random
import functools

@jit
def energy_calculation(location:tuple, lattice_length:int,lattice:np.ndarray) -> int:
    neighbours = [
        (location[0],(location[1]-1)%lattice_length),
        (location[0],(location[1]+1)%lattice_length),
        ((location[0]-1)%lattice_length,location[1]),
        ((location[0]+1)%lattice_length,location[1])
    ]
    energy = 0.0
    for n in neighbours:
        energy+= -1*lattice[location]*lattice[n]
    return energy

@jit
def metropolis_step(lattice:np.ndarray, lattice_length:int, temperature:float) -> np.ndarray:
    location = (randint(0,lattice_length-1),randint(0,lattice_length-1))
    e_before = energy_calculation(location=location,lattice_length=lattice_length,lattice=lattice)
    lattice[location] *= -1
    e_after = energy_calculation(location=location,lattice_length=lattice_length,lattice=lattice)
    delta = e_after-e_before
    r = random()
    factor = np.exp(-(delta/temperature))
    if  delta < 0 or r < factor:
        return lattice
    else: 
        lattice[location] *= -1 #reject flip
        return lattice

@jit
def temperature_step(lattice:np.ndarray, lattice_length:int, temperature:float):
    for i in range(40*lattice_length**2):
        lattice = metropolis_step(lattice=lattice,lattice_length=lattice_length,temperature=temperature)
    return lattice
@jit
def mean_magnet(lattice:np.ndarray, lattice_length:int):
    return np.sum(lattice)/(lattice_length**2)
@jit
def evolution(lattice: np.ndarray,lattice_length):
    lattice_list = []
    temp_list = []
    mag_list = []
    temp = 4.0 # start at 4.0
    temp_min = 0.1
    temp_step = 0.98
    while temp > temp_min:
        temp_list.append(temp)
        lattice = temperature_step(lattice=lattice,lattice_length=lattice_length,temperature=temp)
        lattice_list.append(np.copy(lattice))
        mag_list.append(mean_magnet(lattice=lattice,lattice_length=lattice_length))       
        temp*=temp_step   
    return (lattice_list,temp_list,mag_list)    

def image_creator(frame,results:list,axis):
    axis.clear()
    im = axis.imshow(X=results[frame],origin="lower",cmap="RdBu",vmin=-1,vmax=1)
    return im


if __name__ == "__main__":
    lattice_length = 150
    lattice = np.zeros(shape=(lattice_length,lattice_length),dtype=int)
    for i in range(lattice_length):
        for j in range(lattice_length):
            r = random()
            if r >= 0.5: lattice[i,j] = 1
            else: lattice[i,j] = -1
    lattice_list, temp_list,mag_list = evolution(lattice=lattice,lattice_length=lattice_length)
    fig,axs = plt.subplots()
    axs.scatter(x=temp_list,y=mag_list,s=2)
    axs.set_xlabel("Temperature")
    axs.set_ylabel("Mean Magentisation")
    fig.savefig(fname="magnetisation.png")
    fig1,axs1 = plt.subplots()
    animator = functools.partial(image_creator,results=lattice_list,axis=axs1)
    animation = ani(fig=fig1,func=animator,frames=len(lattice_list),interval=100)
    animation.save(filename="ising_animation.mp4", fps=10,writer="ffmpeg")

    