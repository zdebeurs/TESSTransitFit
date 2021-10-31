#This will contain batman fitting wrapped with emcee
import batman
import numpy as np
import matplotlib.pyplot as plt
import emcee


def fold_data(t,y, period):
  # simple module to fold data based on period

  folded = t%period
  inds = np.array(folded).argsort()
  t_folded = folded[inds]
  y_folded = y[inds]
  return t_folded, y_folded



def logLikelihood(theta, Time, Flux):

    #unpacking the values
    Period, T0, RpRs, aRs, b, u1, u2 = theta


    #Convert b in to inclination
    Inc = np.arccos(b/aRs)
    Inc_Deg = np.rad2deg(Inc)

    global ErrorSumLeast

    if RpRs<0.001 or RpRs>0.3:
        print("RpRs prior")
        return -np.inf

    if aRs<2 or aRs>50.0:
        print("aRs prior")
        return -np.inf

    if b<0 or b>1.2:
        print("b prior")
        return -np.inf



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

    #Subtract the offset
    Offset = np.mean(Flux-ModelFlux)
    ModelFlux+= Offset

    ResidualTransit = ModelFlux - Flux
    ErrorSum = np.sum(ResidualTransit**2)

    if ErrorSum<ErrorSumLeast:
         print("Saved the best figure")
         np.savetxt("Figures/HD189733theta.txt", theta)
         ErrorSumLeast = ErrorSum
         if 1==1:
            #phase fold the transit
            FoldedTime, FoldedFlux = fold_data(Time-T0, Flux, Period)
            _, FoldedModel = fold_data(Time-T0, ModelFlux, Period)
            _, FoldedResidual = fold_data(Time-T0, ResidualTransit, Period)

            SelectIndex = FoldedTime>1.0
            FoldedTime[SelectIndex]-=Period

            ArrangeIndex = np.argsort(FoldedTime)

            FoldedTime = FoldedTime[ArrangeIndex]
            FoldedFlux = FoldedFlux[ArrangeIndex]
            FoldedModel = FoldedModel[ArrangeIndex]



            plt.figure()
            plt.subplot(211)
            plt.plot(FoldedTime, FoldedFlux, "ko")
            plt.plot(FoldedTime, FoldedModel, "r-")
            plt.xlabel("Time")
            plt.ylabel("Flux")
            plt.subplot(212)
            plt.plot(FoldedTime, FoldedFlux-FoldedModel, "ko")
            plt.ylabel("Residual")
            plt.savefig("Figures/HD_189733b_BestFit.png")
            plt.close()


    STDVal = 0.003601548
    ChiSq = ErrorSum/(0.5*STDVal**2.0)
    return -(ChiSq)


def StartSampler():


    #Read LC for now
    Time, Flux = np.loadtxt("HD189733bLC.txt", delimiter=",", unpack=True)

    nWalkers = 30

    #Parameters for Phase curve fitting
    Period_Init = np.random.normal(2.218575200 ,1e-6,nWalkers)
    T0_Init = np.random.normal(2453955.525551, 1e-3, nWalkers)               #planet radius (in units of stellar radii)
    RpRs_Init = np.random.normal(0.15, 0.01, nWalkers)      #planet radius (in units of stellar radii)
    aRs_Init = np.random.normal(8.84, 0.1, nWalkers)       #semi-major axis (in units of stellar radii)
    b_Init = np.random.normal(0.656,0.1,nWalkers)                  #impact parameter
    u1_Init = np.random.normal(0.5,1e-2, nWalkers)               #u1
    u2_Init = np.random.normal(0.2, 1e-2, nWalkers)              #u2


    StartingGuessSelected = np.column_stack((Period_Init, T0_Init, RpRs_Init, aRs_Init, \
                                             b_Init, u1_Init, u2_Init))
    global paramsBatman, ErrorSumLeast
    paramsBatman = batman.TransitParams()     #object to store transit parameters
    ErrorSumLeast = np.inf

    _, nDim = np.shape(StartingGuessSelected)
    sampler = emcee.EnsembleSampler(nWalkers, nDim, logLikelihood, args=[Time,Flux])
    pos, prob, state = sampler.run_mcmc(StartingGuessSelected, 30000, progress=True)

    #Now save the files
    np.save("Figures/HD189733Sample.npy", sampler.chain)



StartSampler()
