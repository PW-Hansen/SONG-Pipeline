import numpy as np
import os
import robust
import astropy.io.fits as pf
import helpful_functions as hf
from scipy.signal import correlate
from scipy.optimize import curve_fit

#%% Importing Arcturus data.
def ArcturusTemplate(program):
    os.chdir(program.settings_path)
    # Importing template for crosscorrelation as well as telluric lines.
    atlas_file = 'ARCTURUS_atlas.fits'   # Add futemll path if needed.
    atlas      = pf.getdata( atlas_file, 1 )
    wavelength = atlas[ 'wavelength']
    sun        = atlas[ 'solarflux']
    arcturus   = atlas[ 'arcturus']
    telluric   = atlas[ 'telluric']
    
    # Correcting a mistake in the Atlas Sun spectrum.
    i_5000=np.where(min(abs(wavelength-5000))==abs(wavelength-5000))[0][0]
    sun[i_5000+3:]=sun[i_5000+2:-1]
    sun[i_5000+2]=0.5*(sun[i_5000+1]+sun[i_5000+3])

    return wavelength,sun,arcturus,telluric

#%% Empirical template.
# Mostly finished, but I never finished implementing it.
def EmpiricalTemplate(data,paths,counts):
    normpath=paths[0]
    
    os.chdir(normpath)
    
    # Finds the spectrum with the highest count.
    spec_max=np.where(max(counts)==counts)[0][0]
    
    norm_fluxes=[]
    norm_lams=[]
    for file in os.listdir():
        if '.fits' in file:
            norm_fluxes.append(pf.getdata(file,1))
            norm_lams.append(pf.getdata(file,2))
        
    xcorr_peak=np.zeros((len(data),20))
    xcorr_peak[spec_max]=0
    
    for i in range(20):
        order=i+10
        for spec_num in range(len(data)):
            if spec_num!=spec_max:
                
                lam_start=norm_lams[spec_num][order][0]
                lam_end=norm_lams[spec_num][order][~np.isnan(norm_fluxes[spec_num][order])][-1]
                lam_mid=(lam_start+lam_end)/2
                lam_num_values=len(norm_fluxes[spec_num][order][~np.isnan(norm_fluxes[spec_num][order])])
                
                lag_vector=np.arange(lam_start-lam_mid,lam_end-lam_mid,(lam_end-lam_start)/lam_num_values/2)
    
                crosscor=correlate(1-norm_fluxes[spec_num][order][~np.isnan(norm_fluxes[spec_num][order])],\
                                    1-norm_fluxes[spec_max][order][~np.isnan(norm_fluxes[spec_max][order])])
                crosscor_middle_index=int((len(crosscor)-1)/2)
                crosscor_central_values=crosscor[crosscor_middle_index-200:\
                                                  crosscor_middle_index+200+1]
                len_diff=len(lag_vector)-len(crosscor)
                    
                lag_vector_middle=lag_vector[crosscor_middle_index-200-len_diff:crosscor_middle_index+200+1-len_diff]
                
                cccv=crosscor_central_values
    
                guess=[max(cccv),lag_vector_middle[np.where(max(cccv)==cccv)[0][0]],0.25,min(cccv)]
                unc=crosscor_central_values-min(crosscor_central_values)+1
                unc=1/unc
                    
                popt,pcov = curve_fit(hf.GaussianOffset,lag_vector_middle,cccv,p0=guess,sigma=unc)
                
                mask=[popt[1]-popt[2]<= lag_value <= popt[1]+popt[2] for lag_value in lag_vector_middle]
                
                popt,pcov = curve_fit(hf.Gaussian,lag_vector_middle[mask],cccv[mask],p0=guess[:3])
                
                xcorr_peak[spec_num,i]=popt[1]
    
    spectrum_shift=[robust.mean(xcorr_peak[spec_num,:]) for spec_num in range(len(data))]    
    min_val=min(spectrum_shift)
    max_val=max(spectrum_shift)
    
    i=0
    order=i+20
    
    spectrum_shift=xcorr_peak[:,i]    
    min_val=min(spectrum_shift)
    max_val=max(spectrum_shift)
    
    min_lam_val=min(norm_lams[spec_max][order])-min_val+0.5
    max_lam_val=max(norm_lams[spec_max][order])-max_val-0.5
    
    interp_lam=norm_lams[spec_max][order]
    interp_lam=interp_lam[interp_lam>min_lam_val]
    interp_lam=interp_lam[interp_lam<max_lam_val]
    
    norm_interp=[]
    for spec_num in range(len(data)):
        lam=norm_lams[spec_num][order]-2*spectrum_shift[spec_num]
        flux=norm_fluxes[spec_num][order]
    
        flux=flux[lam>min_lam_val]
        lam=lam[lam>min_lam_val]
        flux=flux[lam<max_lam_val]
        lam=lam[lam<max_lam_val]    
            
        interp_val=np.interp(interp_lam,lam,flux)
        norm_interp.append(interp_val[~np.isnan(interp_val)])
        
    mask=[min_lam_val<= val <= max_lam_val for val in norm_lams[spec_max][order]]
    
    template_flux_matrix=np.zeros((len(data),len(norm_interp[0])))
    for spec_num in range(len(data)):
        template_flux_matrix[spec_num,:]=norm_interp[spec_num]
    
    
    template_lam=interp_lam
    template_flux=[robust.mean(template_flux_matrix[:,i]) for i in range(np.shape(template_flux_matrix)[1])]
