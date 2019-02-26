from time import sleep

# List of parameters to be used
numberOfCells = 5
gr = 396.41 #solar radius^3/solar mass*seconds
initVelocity = 0
initDensity = 0.000000000000001
initEnergy = 2.1395 * (10 ** (-64)) # 3/2 kT = 2.071 E-16 joules
x0 = 0.5
dx = 0.005
dt = 0.0000001

# List of functions to be used

def calc_pressure(ncells,rho,e): #Calculating pressure
    pcold = 1
    ecold = 1
    p = []
    for i in range(0,ncells):
        p.append(pcold + (2*rho[i]/3)*(e[i] - ecold))

    return(p)

def calc_force(ncells,area,p,x,x0): # Calculating force
    area = [1] + area
    xnew = [x0] + x
    force = []
    
    for i in range(0,ncells):
        if xnew[i] < 1:
            force.append(p[i-1]*area[i] - p[i]*area[i] - gr*(xnew[i]+xnew[i+1]/2))
        else:
            force.append(p[i-1]*area[i] - p[i]*area[i] - gr/((xnew[i]+xnew[i+1])/2)**2)
    print(force)
    return(force)

def calc_velocity(ncells,vel,force,m,dt): # Calculating velocity
    velocity = []
    for i in range(0,ncells):
        velocity.append(vel[i]+force[i]/m[i]*dt)

    return(velocity)

def calc_position(ncells,x,vel,dt): # Calculating pressure
    position = []
    for i in range(0,ncells):
        position.append(x[i] + vel[i]*dt)

    return(position)

def calc_area():
    a=1

    return(a)

def calc_volume(ncells,x,area):
    volume = [x[0]*area[0]]
    for i in range(1,ncells):
        volume.append(x[i]*area[i]-x[i-1]*area[i-1])

    return(volume)

def calc_density(ncells,v):
    density = []
    for i in range(0,ncells):
        density.append(1/v[i])

    return(density)

def calc_energy(ncells,vol_old,vol,p):
    energy = []
    for i in range(0,ncells):
        energy.append((vol[i]-vol_old[i])*p[i])

    return(energy)

#Initialization
def initCells(ncells, velocity0, density0, energy0, xstep,x0):
    area = [1.0001] * ncells
    velocity = [velocity0] * ncells
    density = [density0] * ncells
    energy = [energy0] * ncells
    pressure = calc_pressure(ncells,density,energy)
    position = []
    for i in range(1,ncells+1):
        position.append(x0 + dx*i)
    volume = calc_volume(ncells,position,area)
    mass = []
    for i in range(0,ncells):
        mass.append(volume[i]*density[i])
    return([area,mass,velocity,volume,density,pressure,energy,position])

def iterate(ncells,area,m,vel,vol,rho,p,e,x,dt,x0):
    pressure = calc_pressure(ncells,rho,e)
    force = calc_force(ncells,area,pressure,x,x0)
    velocity = calc_velocity(ncells,vel,force,m,dt)
    position = calc_position(ncells,x,vel,dt,)
    area = [1] * ncells
    volume = calc_volume(ncells,position,area)
    density = calc_density(ncells,volume)
    energy = calc_energy(ncells,vol,volume,pressure)

    return([area,m,velocity,volume,density,pressure,energy,position])

def run(ncells,vel,rho,e,dx,dt,x0):
    init = initCells(ncells,vel,rho,e,dx,x0)
    cont = 1
    while cont:
        area = init[0]
        mass = init[1]
        velocity = init[2]
        volume = init[3]
        density = init[4]
        pressure = init[5]
        energy = init[6]
        position = init[7]
        print('Area: ',init[0], '\nMass: ',init[1], '\nVelocity: ',init[2], '\nVolume: ',init[3], '\nDensity: ',init[4], '\nPressure: ',init[5], '\nInternal Energy: ',init[6], '\nPosition: ',init[7])

        init = iterate(ncells,area,mass,velocity,volume,density,pressure,energy,position,dt,x0)
        cont = input('Type anytihng to contine. \n')
    
    return

run(numberOfCells,initVelocity,initDensity,initEnergy,dx,dt,x0)
