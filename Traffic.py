# traffic - Program to solve the generalized Burger  
# equation for the traffic at a stop light problem

### GITHUB link:###
# https://github.com/abain42/BainALexander_Lab9

# Set up configuration options and special features
import numpy as np
import matplotlib.pyplot as plt


#* Select numerical parameters (time step, grid spacing, etc.).

N = 600
L = 1200      # System size (meters)
h = 2      # Grid spacing for periodic boundary conditions
v_max = 25.    # Maximum car speed (m/s)

tau = h/v_max
print('Last car starts moving after ', (L/4)/(v_max*tau), 'steps')
nstep = 1500
coeff = tau/(2*h)          # Coefficient used by all schemes
xgrid = np.linspace(-600,600,N)

#* Set initial and boundary conditions
rho_max = 1.0                   # Maximum density
Flow_max = 0.25*rho_max*v_max   # Maximum Flow
Flow = np.empty(N)
cp = np.empty(N);  cm = np.empty(N)
# Initial condition is a square pulse from x = -300 to x = 0
rho = np.zeros(N)
x = np.linspace(-600, 600, N)
rho[(x >= -300) & (x <= 0)] = rho_max    # Max density in the square pulse

rho[int(N/2)] = rho_max/2   # Try running without this line
#based on some thought (and googling), I think this is supposed to create a smoother transition from the abrupt change of density at
#x = 0 to x->0+ from rho max to 0,by adding rho = rhomax/2 at x == 0. However, it doesn't seem to visually change the graphs at all

# Use periodic boundary conditions
ip = np.arange(N) + 1  
ip[N-1] = 0          # ip = i+1 with periodic b.c.
im = np.arange(N) - 1  
im[0] = N-1          # im = i-1 with periodic b.c.

#* Initialize plotting variables.
iplot = 1
xplot = (np.arange(N)-1/2.)*h - L/2.    # Record x scale for plot
rplot = np.empty((N,nstep+1))
tplot = np.empty(nstep+1)
rplot[:,0] = np.copy(rho)   # Record the initial state
tplot[0] = 0                # Record the initial time (t=0)


#* Loop over desired number of steps.
for istep in range(nstep) :

    #* Compute the flow = (Density)*(Velocity)
    Flow[:] = rho[:] * (v_max*(1 - rho[:]/rho_max))
  
    #* Compute new values of density using Lax method: 
     
    rho[:] = .5*( rho[ip] + rho[im] ) - coeff*( Flow[ip] - Flow[im] )
    

    #* Record density for plotting.
    rplot[:,iplot] = np.copy(rho)
    tplot[iplot] = tau*(istep+1)
    iplot += 1

#* Graph contours of density versus position and time.
levels = np.linspace(0., 1., num=11) 
ct = plt.contour(xplot, tplot, np.flipud(np.rot90(rplot)), levels) 
plt.clabel(ct, fmt='%1.2f') 
plt.xlabel('x')
plt.ylabel('time')
plt.title('Density contours')
plt.show()

times = [0,250,500,750,1000,1250,1500]
for z in times:
    plt.plot(xplot, rplot[:,z], label =f't={z}')
plt.xlabel('Position x')
plt.ylabel('Density rho(x, t)')
plt.legend()
plt.title('Density Snapshot at Different Times')
plt.show()

### ANSWER TO QUESTION ###
#yes a shock form does appear it occurs at x = 0+ and is because 
#of the abrupt change in density from x = 0 to x > 0 
#the equation put in above was supposed to smooth it out but it 
#seems to have no effect on the appearance of the graphs when commented out
