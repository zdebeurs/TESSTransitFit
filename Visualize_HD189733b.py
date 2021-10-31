import numpy as np
import matplotlib.pyplot as plt
import batman

import matplotlib as mpl
mpl.rc('font',**{'sans-serif':['Helvetica'], 'size':15,'weight':'bold'})
mpl.rc('axes',**{'labelweight':'bold', 'linewidth':1.5})
mpl.rc('ytick',**{'major.pad':22, 'color':'k'})
mpl.rc('xtick',**{'major.pad':10,})
mpl.rc('mathtext',**{'default':'regular','fontset':'cm','bf':'monospace:bold'})
mpl.rc('text', **{'usetex':True})
mpl.rc('contour', **{'negative_linestyle':'solid'})

#Read LC for now
Time, Flux = np.loadtxt("HD189733bLC.txt", delimiter=",", unpack=True)

Period, T0, RpRs, aRs, b, u1, u2 = np.loadtxt("Figures/HD189733theta.txt")

#Evaluate the Model
paramsBatman = batman.TransitParams()

#Convert b in to inclination
Inc = np.arccos(b/aRs)
Inc_Deg = np.rad2deg(Inc)


# Evaluate a batman model
paramsBatman.t0 = T0                        #time of conjunction
paramsBatman.per = Period                   #orbital period
paramsBatman.rp = RpRs                      #planet radius (in units of stellar radii)
paramsBatman.a = aRs                        #semi-major axis (in units of stellar radii)


paramsBatman.inc = Inc_Deg                  #orbital inclination (in degrees)
paramsBatman.ecc = 0                        #eccentricity
paramsBatman.w = 90.0                       #longitude of periastron (in degrees)
paramsBatman.limb_dark = "quadratic"        #limb darkening model
paramsBatman.u = [u1,u2]                    #limb darkening parameters

#Generate the light curve for batman
mTransit = batman.TransitModel(paramsBatman, Time, supersample_factor = 5, exp_time = 0.0014495)#initializes model
ModelFlux =  mTransit.light_curve(paramsBatman)

mTransit = batman.TransitModel(paramsBatman, Time, supersample_factor = 5, exp_time = 0.0014495)#initializes model
ModelFlux =  mTransit.light_curve(paramsBatman)

#Subtract the offset
Offset = np.mean(Flux-ModelFlux)
ModelFlux+= Offset

All_T0_Values = T0+np.arange(-10000,10000)*Period

FirstIndex = np.argmin(np.abs(Time[0]-All_T0_Values))
FirstT0 = All_T0_Values[FirstIndex]

fig, ax = plt.subplots(figsize=(16,8),nrows=2, ncols=6)


for i in range(12):
    Row = i//6
    Column = i%6

    #ax[Row, Column].axvline(x=FirstT0, color="green")
    ax[Row, Column].set_xlim(-0.063*24, 0.063*24)
    ax[Row, Column].set_ylim(0.961, 1.002)
    ax[Row, Column].text(0.40, 0.970, "Epoch "+str(i+1), color="blue")

    if Column ==0:
        ax[Row, Column].set_ylabel("Normalized Flux", fontsize=18)
    else:
        ax[Row, Column].set_yticks([])

    if Row ==1:
        ax[Row, Column].set_xlabel("Time since \n Mid-Transit [Hours]", fontsize=18)
        #ax[Row, Column].set_xticks(rotation=90)
    else:
        ax[Row, Column].set_xticks([])

    if Row==0 and Column==5:
        continue
    ax[Row, Column].text(-1.30, 0.967, "Residual")
    ax[Row, Column].plot((Time-(FirstT0+i*Period))*24, Flux, "ko", markersize=2)
    ax[Row, Column].plot((Time-(FirstT0+i*Period))*24, ModelFlux, "r-")

    ax[Row, Column].plot((Time-(FirstT0+i*Period))*24, Flux-ModelFlux+0.965, "ko", markersize=2)

plt.tight_layout()
plt.subplots_adjust(hspace=0, wspace=0)
plt.savefig("HD189733b_Transits.pdf")
plt.savefig("HD189733b_Transits.png")
