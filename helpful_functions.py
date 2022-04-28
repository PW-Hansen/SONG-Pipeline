# Various helpful functions.
from scipy.signal import correlate
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt

#%% Plotting functions.
# Plots a specific order of a spectrum.
def PlotOrder(data,headers,spec_num,order):
    spec,blaze,wave=GetSpecData(data,headers,spec_num,order)
    plt.figure()
    plt.plot(wave,spec/blaze)    

# Plots a full spectrum.
def PlotSpectrum(data,headers,spec_num):
    plt.figure()
    for order in range(len(data[spec_num][0])):
        spec,blaze,wave=GetSpecData(data,headers,spec_num,order)
        plt.plot(wave,spec/blaze)
        
#%% Various Gaussians.
# Gaussian with offset.
def GaussianOffset(x,a,mu,sig,c):
    return a*np.exp(-np.power(x-mu,2.)/(2*np.power(sig,2.)))+c

# Pure Gaussian.
def Gaussian(x,a,mu,sig):
    return a*np.exp(-np.power(x-mu,2.)/(2*np.power(sig,2.)))

#%% Standard crosscorrelation function.
def CrossCor(flux,template_flux):
    return correlate(1-flux,1-template_flux)


#%% Wavelength-related functions.
# Outputs a vector representing a series of wavelengths between two intervals,
# where the distance between the vector elements is equivalent in terms of the
# radial velocities.
# Radial velocity is assumed to be km/s.
def GetWavelengthIntervals(lam_1,lam_2,rad_vel):
    rad_vel_adj=rad_vel
    c=299792.458
    wavelengthvector=[lam_1]
    while wavelengthvector[-1]<lam_2:
        new_lam=(1+rad_vel_adj/c)*wavelengthvector[-1]
        wavelengthvector.append(new_lam)
    return wavelengthvector[:-1]

# Creates a new set of wavelenghts with distances corresoinding to a specific
# jump in radial velocity through interpolation.
def InterpolatRadVel(rad_vel,lam,flux):
    new_lam=GetWavelengthIntervals(lam[0],lam[-1],rad_vel)
    new_flux=np.interp(new_lam,lam,flux)
    return new_lam,new_flux


# Function to compute the wavelength. Compares the wavelength caliberations
# and does a linear fit based on the times of the before and after 
# calibrations, rather than a simple mean value of them.
def CalcWave(headers,spec_num,order,w0,w1):
    times=[]
    string_keys=['wavebefo','waveafte']
    for key in string_keys:
        tempt=headers[spec_num][key].split('s1_')[1].split('.fits')[0].split('T')
        times.append(tempt[0]+'T'+tempt[1].replace('-',':'))
    t=Time(times,format='isot',scale='utc').mjd
    t=[t[0],headers[spec_num]['mjd-mid'],t[1]]
    t-=t[0]
    t/=t[2]
    
    wdiff=w1-w0
    
    return w0+t[1]*wdiff
    
#%% Functions relating to observational concerns.
# Function to get the "position" from 0 to 1 of where the midpoint of an
# observation happens, relative to the before and after wavelength calibration.
def GetMidScaled(header):
    times=[]
    string_keys=['wavebefo','waveafte']
    for key in string_keys:
        tempt=header[key].split('s1_')[1].split('.fits')[0].split('T')
        times.append(tempt[0]+'T'+tempt[1].replace('-',':'))
    t=Time(times,format='isot',scale='utc').mjd
    t=[t[0],header['mjd-mid'],t[1]]
    t-=t[0]
    t/=t[2]
    
    return t[1]

#%% Extracting the data for a specific order in a specific spectra.
# Gets the flux, blaze, and wavelenght for a specific spectrum and order.
def GetSpecData(data,headers,spec_num,order):
    mid_scaled=[]
    for header in headers:
        try: 
            mid_scaled.append(GetMidScaled(header))
        except:
            mid_scaled.append(0.5)

    spectrum=data[spec_num]
    
    spec0=spectrum[0]
    spec1=spectrum[1]
    blaze=spectrum[2]
    wave0=spectrum[3]
    wave1=spectrum[4]
    
    spec=np.zeros_like(spec0)
    wave=np.zeros_like(wave0)
    
    for spec0val,spec1val,i in zip(spec0[order],spec1[order],range(len(spec0[0]))):
        if np.mean(spec0val)<1:
            spec[order][i]=spec1val
        else:
            spec[order][i]=spec0val
            
    for wave0val,wave1val,i in zip(wave0[order],wave1[order],range(len(wave0[0]))):
        if wave0val!=0 and wave1val!=0:
            wdiff=wave1val-wave0val

            wave[order][i]=wave0val+mid_scaled[spec_num]*wdiff
        elif wave0val!=0:
            wave[order][i]=wave0val
        elif wave1val!=0:
            wave[order][i]=wave1val
        else:
            wave[order][i]=0.0
            
    return spec[order],blaze[order],wave[order]
