from exercise_1 import Cell, Particle, partition
import numpy as np
def no_particle_test() -> bool: 
    A = np.array([])
    s = partition(A, 0, 0, 0.5, 0)
    return s == None 

def all_particles_right() -> bool:
    A = np.array([Particle(np.array([1.0, 1.0]))
                  ,Particle(np.array([0.9,0.3]))
                  ,Particle(np.array([0.7,0.4]))
                  ,Particle(np.array([0.6,0.3]))])
    s = partition(A, 0, len(A)-1, 0.5, 0)
    return s == 0

def inverse_order() -> bool:
    A=np.array([Particle(np.array([0.99,0.99])),
                Particle(np.array([0.88,0.88])),
                Particle(np.array([0.66,0.66])),
                Particle(np.array([0.44,0.44])),
                Particle(np.array([0.33,0.33])),
                Particle(np.array([0.11,0.11]))])
    s = partition(A, 0, len(A)-1,0.5,0)
    return s==2

def already_ordered() -> bool:
    A = np.array([Particle(np.array([0.1,0.9])),
    Particle(np.array([0.2,.9])),
    Particle(np.array([0.7,0.9])),
    Particle(np.array([0.8,0.9]))])
    s = partition(A, 0, len(A)-1,0.5,0)
    return s == 1

def testrunner(): 
    print("running 4 tests:")
    print(f"no particle test: {no_particle_test()}")
    print(f"all particles x > 0.5 test: {all_particles_right()}")
    print(f"inverse order test: {inverse_order()}")
    print(f"already ordered test: {already_ordered()}")
    
if __name__ == "__main__":
    testrunner()