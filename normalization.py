import numpy as np
import helpful_functions as hf
import write_functions as wf
import robust
import os
import matplotlib.pyplot as plt

#%% Dataset normalization.
def NormalizeAllData(data,headers,program,paths,template,wavelength,noise,method):
    # Grabbing data paths.
    normpath=paths[0]
    ccfspath=paths[2]

    # Grabbing constants.
    deltav=float(program.settings['deltav'])
    lag_max=float(program.settings['lag_max'])

    # Normalizing. Always checks whether the data has already been processed,
    # and if that is not the case, performs normalization and asves the data.
    print('Normalizing the spectra and crosscorrelating')
    os.chdir(ccfspath)
    for spec_num in range(len(data)):
        # Gets the file name from the header.
        file_name=headers[spec_num]['file']
        file_name=file_name.replace('.fits','')+'_ccfs.fits'
        print('Processing '+file_name+', spec_num '+str(spec_num))
        if file_name in os.listdir():
            print(file_name+' already processed.')
    
        # If data has not been processed, do so and save it.
        else:       
            # Sets up a matrix for the spectrum containing order-wise lists of 
            # flux, lambda, template flux, and crosscorrelation values.
            # Simply set to have shapes (1,1) at the start, as the size of the
            # length will be modified later. This later length will vary from
            # spectrum to spectrum based on the wavelength scaled.
            fluxs       =np.zeros((1,1))
            lams        =np.zeros((1,1))
            templates   =np.zeros((1,1))
            crosscors   =np.zeros((1,1))
            matrices=[fluxs,lams,templates,crosscors]
                
            # Runs through the orders for the spectrum starting with the
            # highest order first. This is because the higher order will have
            # the longest lenghts of variables, as having the wavelength
            # differences being equal in terms of radial velocity rather than
            # wavelength makes for shorter steps at higher wavelenghts, 
            # leading to more values across an otherwise equally large span
            # in wavelength.
            for order in range(len(data[0][0]))[::-1]:
                matrices=NormalizeSpectrum(data,headers,spec_num,order,deltav,lag_max,template,wavelength,matrices,noise,method)
            
            wf.NormalizationWriteFits(headers[spec_num],matrices,normpath,ccfspath)


#%% Normalization procedure for a spectrum.
# Normalizes the data and outputs interpolated flux, wavelength, and template
# values according to a specific radial velocity stepsize, as well as the
# crosscorrelation function.
# The method used is left to the user, but as of writing this comment, there
# is only a single one available.
def NormalizeSpectrum(data,headers,spec_num,order,deltav,lag_max,template,wavelength,matrices,noise,method):
    if method=='median':
        lam,flux=NormalizationMethodMedianThreshold(data,headers,spec_num,order,deltav,noise)
    
    # Extracts the values from the matrix.
    fluxs=matrices[0]
    lams=matrices[1]
    templates=matrices[2]
    crosscors=matrices[3]

    # Removing the edges of the order, as these are generally more volatile
    # then more central values in the order.
    flux=flux[lam>min(lam)+2]
    lam=lam[lam>min(lam)+2]
    flux=flux[lam<max(lam)-2]
    lam=lam[lam<max(lam)-2]
    
    # Checks whether the extracted matrix values has already been transformed. 
    # If they has not, which should only occur for the highest order in the 
    # spectrum, then tthey are ransformed so that they can to contain the 
    # flux, lambda, and template wavelength.
    # +10 is just to be on the safe side, as there have been occassions during
    # testing where a lower order may have more values than a higher order, 
    # due to different cutoffs, but this effect should never result in a 
    # difference of more than 10.
    if np.shape(fluxs)==(1,1):
        fluxs       =np.empty((len(data[0][0]),len(flux)+10)); fluxs[:]    =np.NaN
        lams        =np.empty((len(data[0][0]),len(flux)+10)); lams[:]     =np.NaN
        templates   =np.empty((len(data[0][0]),len(flux)+10)); templates[:]=np.NaN
    
    # Interpolates flux values from the template to match onto the wavelength
    # values for the analyzed data.
    template_flux=np.interp(lam,wavelength,template)

    # Takes the order-wise values for flux, lambda, and template flux values
    # and throws them into the matrix for the spectrum.
    fluxs[order,:len(lam)]=flux
    lams[order,:len(lam)]=lam
    templates[order,:len(lam)]=template_flux
        
    # Gets the cross-correlation values as well.
    crosscor=hf.CrossCor(flux,template_flux)
    if np.shape(crosscors)==(1,1):
        crosscors   =np.empty((len(data[0][0]),len(crosscor)+10)); crosscors[:]=np.NaN
    crosscors[order,:len(crosscor)]=crosscor

    # Gathers the new matrix values and returns them to the parent function.
    matrices[0]=fluxs
    matrices[1]=lams
    matrices[2]=templates
    matrices[3]=crosscors        
        
    return matrices


#%% Normalization methods.

# Normalization median threshold method.
# This method attempts to find a suitable second-order polynomial which will
# result in most of the absorption lines not being considered for further
# normalization. It does this by first naievely normalizing the entire
# blaze-corrected data, then taking the median of these normalized values, 
# and then taking the median values of the values above the initial median 
# values.
# Then, this threshold is applied as a mask and a second-order polynomial is
# fitted to the remaining values.
def NormalizationMethodMedianThreshold(data,headers,spec_num,order,deltav,noise='silent'):
    # Grabs the data.
    spec,blaze,wave=hf.GetSpecData(data,headers,spec_num,order)    
    blazecorrected=spec/blaze

    wave=wave[~np.isnan(blazecorrected)]
    blazecorrected=blazecorrected[~np.isnan(blazecorrected)]
    
    # Getting rid of values which are entirely too high.
    mean=robust.mean(blazecorrected)
    std=robust.std(blazecorrected)
    
    wave=wave[blazecorrected<mean+5*std]
    blazecorrected=blazecorrected[blazecorrected<mean+5*std]
    
    # As a first pass, performs a normalization with the assumption that the
    # data fits a second degree polynomial.
    x=np.polyfit(wave,blazecorrected,2)
    p=np.poly1d(x)
    if min(p(wave))<0:
        p[0]-=2*min(p(wave))
    normalized=blazecorrected/p(wave)
            
    # Determines the median value of the quick-and-dirty normalization, then
    # takes a new median of values above the first median and uses this as a
    # threshold value.
    median=np.median(normalized)
    threshold=np.median(normalized[normalized>median])
    
    # Another second-degree polynomial is fit to the data, but this time, only
    # values which exceeded the threshold after being normalized above are
    # included a fitting points.
    x2=np.polyfit(wave[normalized>threshold],blazecorrected[normalized>threshold],2)
    p2=np.poly1d(x2)
    new_normalized=blazecorrected/p2(wave)
    
    # The newly normalized data are interpolated to specific steps 
    # corresponding to equal shifts in radial velocity.
    new_lam,new_flux=hf.InterpolatRadVel(deltav,wave,new_normalized)

    # Potentially interesting plots for debugging purporses.
    if noise!='silent':
        plt.figure()
        plt.plot(wave,blazecorrected,'b--')
        plt.plot(wave,p(wave),'k-')
        plt.grid()
        
        plt.figure()
        plt.plot(wave,normalized,'b--')
        plt.plot(wave,[threshold]*len(wave),'k-')
        plt.grid()

        plt.figure()
        plt.plot(wave,blazecorrected,'b-')
        plt.plot(wave[normalized>threshold],blazecorrected[normalized>threshold],'k.')
        plt.grid()
        
        plt.figure()
        plt.plot(new_lam,new_flux,'b-')
        plt.grid()
    
    
    return np.array(new_lam),new_flux
