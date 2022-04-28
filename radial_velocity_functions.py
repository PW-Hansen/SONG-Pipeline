import numpy as np
import helpful_functions as hf
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#%% Radial velocity determination functions.
# Function to fit a Gaussian to a zoomed-in section of the cross-correlation
# function. This function is preparatory, and gives a rough estimate for where
# the the peak of the CCF is.
def RadVelDeterminationBroad(deltav,lag_max,crosscor,noise='silent'):
    # Sets up the values that are examined from the CCF.
    lag_vector=np.arange(-lag_max,lag_max+deltav,deltav)
    crosscor_middle_index=int((len(crosscor)-1)/2)
    crosscor_central_values=crosscor[crosscor_middle_index-int(lag_max/deltav):crosscor_middle_index+int(lag_max/deltav)+1]

    # For debugging purposes, the CCF can be plotted.
    if noise!='silent':
        plt.figure()
        plt.grid()
        plt.xlabel('Displacement (km/s)')
        plt.ylabel('Cross-correlation signal')
        plt.plot(lag_vector,crosscor_central_values,'k-')
    
    # First pass. Fits a Gaussian with an offset to the CCF. 
    init_guesses=[max(crosscor_central_values),lag_vector[np.where(crosscor_central_values==max(crosscor_central_values))][0],lag_max/20,np.mean(crosscor_central_values)]
    limits=([0,init_guesses[1]-1,1,min(crosscor_central_values)],\
             [init_guesses[0],init_guesses[1]+1,init_guesses[2]*2,max(crosscor_central_values)])

    # Uncertainty is set to bias the Gaussian to fit to the largest peak above
    # smaller vlaues.
    unc=crosscor_central_values-min(crosscor_central_values)+1
    unc=1/unc
    
    popt,pcov = curve_fit(hf.GaussianOffset,lag_vector,crosscor_central_values,p0=init_guesses,bounds=limits,sigma=unc)
    
    # Plotting still for debugging purposes.
    if noise!='silent':
        plt.figure()
        plt.grid()
        plt.xlabel('Displacement (km/s)')
        plt.ylabel('Cross-correlation signal')
        plt.plot(lag_vector,crosscor_central_values,'k-')
        plt.plot(lag_vector,hf.GaussianOffset(lag_vector,*popt),'b-')

    # Second pass. Examines values within one sigma from the expected peak.
    second_indices_mask=np.logical_and(lag_vector<popt[1]+popt[2]*2, lag_vector>popt[1]-popt[2]*2)
    
    second_lag=lag_vector[second_indices_mask]
    second_crosscor=crosscor_central_values[second_indices_mask]
    
    init_guesses=[max(second_crosscor),lag_vector[np.where(crosscor_central_values==max(second_crosscor))][0],len(second_lag)/4]
    
    popt,pcov = curve_fit(hf.Gaussian,second_lag,second_crosscor,p0=init_guesses)
    
    # Plotting still possible for debugging purposes.
    if noise!='silent':
        plt.figure()
        plt.grid()
        plt.xlabel('Displacement (km/s)')
        plt.ylabel('Cross-correlation signal')
        plt.plot(second_lag,second_crosscor,'k-')
        plt.plot(second_lag,hf.Gaussian(second_lag,*popt),'b-')
    
    return popt[1]

# Function to fit a Gaussian to a zoomed-in section of the cross-correlation
# function. 
def RadVelDeterminationNarrow(deltav,lag_max,crosscor,guess_peak,guess_width,noise='silent'):
    # Sets up values.
    lag_vector=np.arange(-lag_max,lag_max+deltav,deltav)
    crosscor_middle_index=int((len(crosscor)-1)/2)
    crosscor_central_values=crosscor[crosscor_middle_index-int(lag_max/deltav):crosscor_middle_index+int(lag_max/deltav)+1]

    # Still debugging.
    if noise!='silent':
        plt.figure()
        plt.plot(lag_vector,crosscor_central_values,'k-')

    # Finds all values within two sigma of the predicted Gauss peak, takes
    # the highest value, and creates a new mask within one sigma of that 
    # value, then fits a Gaussian to those values.
    mask_initial=[guess_peak-2*guess_width <= lag_value <= guess_peak+2*guess_width for lag_value in lag_vector]
    new_peak=lag_vector[np.where(max(crosscor_central_values[mask_initial])==crosscor_central_values)[0][0]]
    mask_local=[new_peak-guess_width <= lag_value <= new_peak+guess_width for lag_value in lag_vector]

    if noise!='silent':
        plt.figure()
        plt.plot(lag_vector[mask_local],crosscor_central_values[mask_local],'k-')

    init_guesses=[max(crosscor_central_values[mask_local]),new_peak,5]   

    popt,pcov = curve_fit(hf.Gaussian,lag_vector[mask_local],crosscor_central_values[mask_local],p0=init_guesses)
    
    return crosscor_central_values,popt
