import numpy as np
import matplotlib.pyplot as plt
from scipy import special
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix' # math fontの設定
plt.rcParams['font.size'] = 25

rhoInf = 3.058e-5 # kg/m^3
veloInf = 6880.0 # m/s
radiStar = 1.8e-3 # m
tempInf = 219. # K
tempWall = 1.0e3 # K
omega = 2./3. # mu \propto Temp^omega
qInf = 0.5*rhoInf*(veloInf**3)
gamma = 1.403
mAir = 28.966e-3 # kg/mol
gasConst = 8.314462 # gas constant
boltz = 1.380e-23 # boltzman cont. J/K
#mu = 1.458e-6*(tempInf**1.5)/(tempInf + 110.4) # Sutherland
mu = 1.716e-5*((tempInf/273.0)**(2./3.)) # power law by Maxwell, Rayleigh (for rarefied gas)
veloSound = np.sqrt(gamma*gasConst*tempInf/mAir) # a = sqrt(gamma*R*T/M)
mach = veloInf/veloSound
reynolds = rhoInf*2.0*radiStar*veloInf/mu
w_r = (mach**(2.*omega))/reynolds # rarefied parameter

alpha = -0.476
c2 = 1.10
c3 = 1.8

qNS = (7.455E-5)*(rhoInf**0.4705)*(veloInf**3.089)*(radiStar**(-0.52)) # W/m^2
hct = qNS/(qInf*(1.-alpha*w_r))
s_0 = np.sqrt(0.5*gamma)*mach
qfm = rhoInf*(np.sqrt(2.*boltz*tempInf/mAir)**3.)*0.25/np.sqrt(np.pi)*((s_0**2.+gamma/(gamma-1.)-0.5*(gamma+1)/(gamma-1)*tempWall/tempInf)*(np.exp(-s_0**2.)+np.sqrt(np.pi)*s_0*(1+special.erf(np.sqrt(np.pi)*s_0)))-0.5*np.exp(-s_0**2.))
hfm = qfm/qInf

h_0 = hct+(hfm-hct)*max(((w_r-c2)/(w_r+c3)),0.0)

def htheta(theta):
    htheta = (0.55+0.45*np.cos(2.*np.deg2rad(theta))+(w_r/3.)*np.cos(np.deg2rad(theta)))/(1+(w_r/3.))
    #htheta = htheta*h_0*qInf # coeff -> heat flux
    #htheta = htheta*h_0 # heat flux coefficient
    return htheta

surfData = np.loadtxt('surf.50000',delimiter=' ',skiprows=9,usecols=[0,1,2,3])
#cfdData = np.loadtxt('qw',delimiter=' ',skiprows=1)
cfdData = np.loadtxt('qw',delimiter=' ',skiprows=1)
imax = len(surfData)
cfdImax = len(cfdData)
surfData = surfData.T
cfdData = cfdData.T
heatIndex = 1 # use index
heatMax = max(surfData[heatIndex])
cfdHeatMax = max(cfdData[heatIndex])
cfdXMax = max(cfdData[0])
print(cfdImax)
cfdX = np.linspace(0,180,cfdImax)

for i in range(imax):
    surfData[0,i] = 180.*float(imax-surfData[0,i])/float(imax-1)
    #surfData[heatIndex,i] = surfData[heatIndex,i]/heatMax
    surfData[heatIndex,i] = surfData[heatIndex,i]/1e6 # MW/m^2
    #surfData[heatIndex,i] = surfData[heatIndex,i]/qInf

for i in range(cfdImax):
    cfdData[0,i] = 180.*cfdData[0,i]/cfdXMax
    #cfdData[heatIndex,i] = cfdData[heatIndex,i]/cfdHeatMax
    cfdData[heatIndex,i] = cfdData[heatIndex,i]/1e6 # MW/m^2
    #cfdData[heatIndex,i] = cfdData[heatIndex,i]/qInf

print("h0 = "+str(h_0)+", qNS = "+str(qNS)+", qInf = "+str(qInf)+", Wr = "+str(w_r))

#print(surfData)

x = np.linspace(0,180.,1e2)

#plt.style.use('dark_background')
plt.figure(figsize=(12,8))
plt.xlim(0,180)
plt.ylim(0,8.5)
plt.xticks(np.arange(0, 180 + 1, 30))
plt.yticks(np.arange(0, 9 + 1, 1))
plt.xlabel('Surface Angle [deg.]')
plt.ylabel('Heat Flux [$\\mathrm{MW/m^2}$]')
plt.hlines(8.057, 0, 180, 'k', linestyles='dashed', label='Stag. Heat Flux (DKR)')
#plt.plot(x,htheta(x),label='Sighn et. al.',linewidth=4)
plt.plot(surfData[0],surfData[heatIndex],marker='x',linestyle='None',label='DSMC result',color='red')
plt.plot(cfdX,cfdData[heatIndex],marker='.',linestyle='None',label='Navier–Stokes result',color='#3180b6')
#plt.tick_params(colors='white')
plt.legend()
plt.savefig('heat-flux-magni.pdf',format='pdf',dpi=1200,transparent=True,bbox_inches='tight',pad_inched=0)
plt.show()
