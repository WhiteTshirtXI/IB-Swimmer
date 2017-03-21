# Simulation of infinite sheet problem for comparison with asymptotic results
#Last updated Jan 17,2017

# Import Packages
import sys, petsc4py, os
arg = sys.argv
arg.extend(['-fluid_solver', 'Minev_Penalty', '-momentum_solver','Thomas','-penalty_solver', 'Thomas', '-transpose_plan','FFTW_PATIENT'])
petsc4py.init(arg)

if os.environ["PETSCIBM_DIR"]+'/solvers/' not in sys.path:
   sys.path.append( os.environ["PETSCIBM_DIR"]+'/solvers/' )
   
from petsc4py import PETSc
from mpi4py import MPI
import numpy as np
import PeriodicIBSolver as ibSolver
import PlottingUtils as utils
import pylab as plt

# Get number of procs in each direction
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# Parameters
rho = 1.0           # Density
mu = 0.01           # Viscosity

Omega = 8 * np.pi 
Tfinal = 1.5                      # Final Time

DomainLength = np.array([ 0.2, 0.2]) # Domain Dimensions

N = 64
dx =  1.0/N
dt = 1e-6                   # Time-step
#Tfinal = 0.5
#Tfinal =  25*dt
NN = np.array([N,N])      # Fluid Grid-points in each direction
Ns = 128   # Number of IB points
   
sigma_spring = 1e4           # Spring Constant
sigma_bending = 0.01        # Bending Constant


output_folder = 'TuckFormula0.005'
output_dt = round(.05/dt)*dt
#output_dt = dt
#print output_dt
output_times = np.arange(0,Tfinal, output_dt)


#
# Construct Membrane
#

Kappa =  10* np.pi 
L =    0.2

#IC_ChiX = lambda (S): S + 6.0
IC_ChiX = lambda (S): S 
IC_ChiY = lambda (S): DomainLength[1]/2 + 0.0053 * np.sin(Kappa * S)


DSS_TARGET_X = lambda (S): 0.0 * np.cos(S)
DSS_TARGET_Y = lambda (S): -0 * 5 * 0.2 * (np.pi/DomainLength[0]/2)**2 * np.sin(S*np.pi/DomainLength[0]/2)
ds = L/Ns
    
Data = {}
Data["Dimensions"] = 2
Data["Points"] = []
Data["ForceConnections"] = []
    
S = np.arange(0.0, L, ds)
IC_ChiX = IC_ChiX(S)
IC_ChiY = IC_ChiY(S)
print IC_ChiX
print IC_ChiY
DSS_TARGET_X = DSS_TARGET_X(S)
DSS_TARGET_Y = DSS_TARGET_Y(S)
print S
#A =  0.1 * (np.pi/DomainLength[0]/2)**2
K =   10 * np.pi
A = 0.0053 * K**2
print IC_ChiX, IC_ChiY, S
for i in range(Ns):
   Data["Points"].append({"PointID" : i, "X1" : IC_ChiX[i], "X2" : IC_ChiY[i]})
   
   R1 = i+1
   R2 = i+2
   L1 = i-1
   L2 = i-2
   if i == Ns-2:
      R1 = Ns-1 
      R2 = 0
   
   if i == Ns-1: 
      R1 = 0 
      R2 = 1

   if i == 0: 
      L1 = Ns-1
      L2 = Ns -2

   if i == 1:
      L1 = 0 
      L2 = Ns -1

   #Data["ForceConnections"].append({"ForceID" : "BENDING", "FCID" : 2*i, "PointID" : i, "L2PointID" : i-2, "LPointID" : i-1, "R2PointID" : R2, "RPointID" : R1, "sigma" : sigma_bending, "D2X_TARGET_1" : DSS_TARGET_X[i], "D2X_TARGET_2" : DSS_TARGET_Y[i], "D2X_LTARGET_1" : DSS_TARGET_X[i-1], "D2X_LTARGET_2" : DSS_TARGET_Y[i-1], "D2X_RTARGET_1" : DSS_TARGET_X[(i+1)%Ns], "D2X_RTARGET_2" : DSS_TARGET_Y[(i+1)%Ns], "ds" : ds, "dv" : ds})
   #if i<Ns/2: Data["ForceConnections"].append({"ForceID" : "BENDING", "FCID" : 2*i, "PointID" : i, "L2PointID" : i-2, "LPointID" : i-1, "R2PointID" : R2, "RPointID" : R1, "sigma" : sigma_bending, "A" : A, "K":K, "Omega":Omega,  "ds" : ds, "dv" : ds})
   #if i>=Ns/2: Data["ForceConnections"].append({"ForceID" : "BENDING", "FCID" : 2*i, "PointID" : i, "L2PointID" : i-2, "LPointID" : i-1, "R2PointID" : R2, "RPointID" : R1, "sigma" : sigma_bending, "A" : A/2, "K":K, "Omega":Omega,  "ds" : ds, "dv" : ds})
   Data["ForceConnections"].append({"ForceID" : "BENDING", "FCID" : 2*i, "PointID" : i, "L2PointID" : L2, "LPointID" : L1, "R2PointID" : R2, "RPointID" : R1, "sigma" : sigma_bending, "A" : A, "K":K, "Omega":Omega,  "ds" : ds, "dv" : ds})
   Data["ForceConnections"].append({"ForceID" : "ELASTIC", "FCID" : 2*i+1, "PointID" : i, "LPointID" : L1, "RPointID" : R1, "sigma" : sigma_spring, "L" : 1.0, "ds" : ds, "dv" : ds})

#diff = (jellyfish_width - jellyfish_width_contracted)

#Data["ForceConnections"].append({"ForceID" : "ELASTIC_OSCILLATINGSTIFFNESS_TWOPOINT", "FCID" : 2*Ns, "PointID1" : int(Ns/8), "PointID2" : int(7*Ns/8), "sigma" : 0.0, "L" : jellyfish_width_contracted, "OscillatingTypeID" : 2, "OscillatingFrequency" : 1.0/(contraction_time+expand_time), "OscillatingAmplitude" : sigma_contraction, "onoff_ratio" : contraction_time/(contraction_time+expand_time)})

#Data["ForceConnections"].append({"ForceID" : "ELASTIC_OSCILLATINGSTIFFNESS_TWOPOINT", "FCID" : 2*Ns+1, "PointID1" : 0, "PointID2" : Ns-1, "sigma" : 0.0, "L" : jellyfish_width_contracted, "OscillatingTypeID" : 2, "OscillatingFrequency" : 1.0/(contraction_time+expand_time), "OscillatingAmplitude" : sigma_contraction, "onoff_ratio" : contraction_time/(contraction_time+expand_time)})
#Data["ForceConnections"].append({"ForceID" : "ELASTIC_OSCILLATINGLENGTH_TWOPOINT", "FCID" : 2*Ns, "PointID1" : 0, "PointID2" : Ns-1, "sigma" : sigma_contraction, "L" : jellyfish_width_contracted, "OscillatingTypeID" : 3, "OscillatingFrequency" : 1.0/(contraction_time+expand_time), "OscillatingAmplitude" : diff, "onoff_ratio" : contraction_time/(contraction_time+expand_time)})
  
Re = rho * Omega/(mu*K**2)

# Run IB Problem
ibSolver.Run2d(rho, mu, NN, DomainLength, Data, dt, output_times, output_folder)

print Re
# Plot Movie
if rank == 0:
   movie = utils.IBMovie(output_folder+"/movie")
   max_i = utils.MaxIndexInIBData( output_folder )
   plt.clf()

   fig_width_pt = 400.0  # Get this from LaTeX using \showthe\columnwidth
   inches_per_pt = 1.0/72.27               # Convert pt to inch
   fig_width = fig_width_pt*inches_per_pt  # width in inches
   fig_size =  [3*fig_width,fig_width]
   params = {'backend': 'ps',
             'font.size': 8,
             'axes.labelsize': 22,
             'text.fontsize': 22,
             'font.family':'Times',
             'legend.fontsize': 10,
             'xtick.labelsize': 22,
             'ytick.labelsize': 22,
             'lines.markersize': 3,
             'lines.linewidth':3,
             'text.usetex': True,
             'savefig.dpi':200,
             'figure.figsize': fig_size
   }
   plt.rcParams.update(params)
   plt.clf()
   # plt.hold(True)
   # X, Y =[], []
   # for i in range(1,max_i+1):
   #     [ t, x, y, vx, vy, p, fx, fy, ibx, iby] = utils.Retrieve2dIBData( i, output_folder )
   #     X.append(ibx[Ns-1])
   #     Y.append(iby[Ns-1])
   
   # plt.plot(X,Y)
   # #plt.axis([0 ,DomainLength[1],0,DomainLength[0]])
   # plt.title('Swimmer head position Tfinal=10, Period of beating = 1/3,\n Initial position=(%f,%f)'%(X[1],Y[1]))
   # movie.plotframe()

   # for i in range(1,max_i+1):
   #    [ t, x, y, vx, vy, p, fx, fy, ibx, iby] = utils.Retrieve2dIBData( i, output_folder )
   #    plt.plot(ibx,iby,'.',color='#a64527')
   #    plt.axis([0 ,DomainLength[1],0,DomainLength[0]])
   #    movie.plotframe()
   sum = 0   
   for i in range(1,max_i+1):
      [ t, x, y, vx, vy, p, fx, fy, ibx, iby] = utils.Retrieve2dIBData( i, output_folder )
      [ vx, vy] = utils.RetrieveMAC2dIBData( i, output_folder )
      vorticity = utils.VelocityMACToCenterCurl( vx, vy, np.array([dx,dx]) )
      #vorticity = (np.abs(vorticity)>10)*vorticity
      #print "%e "%(np.max(np.abs(vorticity)))
#     avg_vel = np.average(vx)
#     Re = avg_vel*rho*L/mu
      #print Re
      #sum = sum + ibx[0]
#sum = sum + vx[0]
      #print ibx[Ns-1]
      center_iby = np.average(iby)

      #[ x, y, vorticity] = [ utils.Refine2dVariable( x, 256 ), utils.Refine2dVariable( y, 256 ), utils.Refine2dVariable( vorticity, 256 )] 
      Q = plt.quiver( x, y, vx, vy, headlength=2, headwidth=2, scale=250.0, color='0.5')
      # #Q = plt.quiver( x, y, fx, fy, headlength=1, headwidth=1, scale=5000.0, color='0.5')

      # #print "%e %e"%(t,(np.max(iby)-np.min(iby))/(np.max(ibx)-np.min(ibx)))
      # #plt.pcolor(y, x, vorticity,cmap=plt.cm.Blues,vmin=-5.0e+02,vmax=5.0e+02)
      plt.imshow(vorticity, interpolation='bilinear', cmap=plt.cm.Blues, origin='lower',vmin=-5.0e+02, vmax=5.0e+02,extent=[y[0,0],y[-1,0],x[0,0],x[0,-1]])
      plt.clim()
      plt.plot(ibx,iby,'.',color='#a64527')
      #plt.axis([center_iby-9.0,center_iby+3.0,DomainLength[0]/3,2*DomainLength[0]/3])
      plt.axis([0 ,DomainLength[1],0,DomainLength[0]])
      plt.title('t=%f'%(t))
      plt.rcParams.update(params)

   #movie.plotframe()
      movie.plotframe()
   #print sum   
   movie.render(output_folder+"/movie",10)


