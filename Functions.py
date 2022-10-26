from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter, medfilt
from scipy.interpolate import interp1d
#Need flattening function

#Read the file
def Load_LC(AllLoc):
    '''
    Parameters:
                AllLoc: list of string
                        provide locations as an list


    Returns:
            Time and flux as arrays
    '''

    AllTime = []
    AllFlux = []

    #Loop over the file locations
    for Loc in AllLoc:
        File = fits.open(Loc)
        CurrentTime = File[1].data['TIME']
        CurrentFlux = File[1].data['SAP_FLUX']


        AllTime.extend(CurrentTime)
        AllFlux.extend(CurrentFlux/np.nanmedian(CurrentFlux))

    AllTime = np.array(AllTime)
    AllFlux = np.array(AllFlux)

    NanIndex = np.logical_or(np.isnan(AllTime), np.isnan(AllFlux))

    AllTime = np.ascontiguousarray(AllTime[~NanIndex])
    AllFlux = np.ascontiguousarray(AllFlux[~NanIndex])



    return AllTime+2457000, AllFlux


def NormalizeFlux(AllTime, AllFlux, TransitIndex):

    #ModelFlux = savgol_filter(AllFlux[~TransitIndex], 101, 3)
    ModelFlux = medfilt(AllFlux[~TransitIndex], 31)
    Interpolator  = interp1d(AllTime[~TransitIndex], ModelFlux)
    InterpolatedFlux = Interpolator(AllTime)

    NormalizedFlux = AllFlux/InterpolatedFlux


    plt.figure()
    plt.subplot(211)
    plt.plot(AllTime, AllFlux, "ko")
    plt.plot(AllTime, InterpolatedFlux, "r-")
    plt.ylabel("Flux")
    plt.xlabel("Time")
    plt.subplot(212)
    plt.plot(AllTime, AllFlux/InterpolatedFlux, "ko")
    plt.ylabel("Normalized Flux")
    plt.xlabel("Time")
    plt.show()

    return NormalizedFlux
