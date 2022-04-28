import radial_velocity_functions as rvf
import helpful_functions as hf
import write_functions as wf
import numpy as np
import os
import astropy.io.fits as pf
from scipy.optimize import curve_fit
from astropy.time import Time
from barycorrpy import get_BC_vel


# Grabs the center values of the crosscorrelations, then divides them
# by the maximum value, and combines everything to provide an initial
# guess for more detailed fitting.
def ScaleCCF(crosscors,lag_max,deltav):
    cc_scaled=[]
    for cc in crosscors:
        crosscor_temp=cc[~np.isnan(cc)]
        midindex=int((len(crosscor_temp)-1)/2)
        cc_mid=crosscor_temp[midindex-int(lag_max/deltav):midindex+int(lag_max/deltav)+1]
        cc_scaled.append(cc_mid/max(cc_mid))
    
    cc_scaled_combined=[0]*len(cc_scaled[0])
    for cc_scale in cc_scaled:
        cc_scaled_combined+=cc_scale
    return cc_scaled_combined


#%%
# Computes the radial velocity for every order in every spectrum.
def GetRadialVels(data,headers,program,paths,noise):
    # Grabbing data paths.
    gauspath=paths[1]
    ccfspath=paths[2]

    # Grabbing constants.
    deltav=float(program.settings['deltav'])
    lag_max=float(program.settings['lag_max'])

    # Checks whether fitting to CCF has already occurred. If not, done here using
    # methods outlined above.
    print('Fitting to crosscorrelation')
    vel_diff_matrix=np.zeros((len(data),len(data[0][0])))
    
    BC_corrs=[]
    
    for spec_num in range(len(data)):
        # Gets the file name from the header.
        os.chdir(gauspath)
        file_name=headers[spec_num]['file']
        print('Processing '+file_name+', spec_num '+str(spec_num))
        
        file_name_split=file_name.split('.')
    
        if file_name_split[0]+'_gaussparams.'+file_name_split[1] in os.listdir():
            print(file_name+' already processed, loading prior results')
            gaus_data=pf.getdata(file_name.replace('.fits','')+'_gaussparams.fits',2)
            vel_diff_matrix[spec_num]=gaus_data[:,1]
        # If data has not been processed, do so and save it.
        else:
            print('Crosscorrelating '+file_name)
            os.chdir(ccfspath)
            crosscors=pf.getdata(file_name.replace('.fits','')+'_ccfs.fits')
            
            # First pass. Gets an approximate value for the CCF peak.
            for order in range(len(data[0][0])):
                crosscor=crosscors[order][~np.isnan(crosscors[order])]            
                vel_diff=rvf.RadVelDeterminationBroad(deltav,lag_max,crosscor,noise='silent')
                vel_diff_matrix[spec_num,order]=vel_diff
            
            # Scales crosscorrelation functions.
            cc_scaled_combined=ScaleCCF(crosscors,lag_max,deltav)
    
            lag_vector=np.arange(-lag_max,lag_max+deltav,deltav)
            cc_peak_index=np.where(max(cc_scaled_combined)==cc_scaled_combined)[0][0]
            gauss_guess=[max(cc_scaled_combined)-min(cc_scaled_combined),lag_vector[cc_peak_index],10,min(cc_scaled_combined)]
            popt,pcov = curve_fit(hf.GaussianOffset,lag_vector,cc_scaled_combined,p0=gauss_guess)
                    
            ccfs=np.zeros((len(data[0][0]),int(1+lag_max/deltav*2)))
            gaus=np.zeros((len(data[0][0]),3))
            
            # Second pass. Zooms in on the peak found above and fits only to
            # values close to it.
            for order in range(len(data[0][0])):
                crosscor=crosscors[order][~np.isnan(crosscors[order])]      
                
                cross_cor_data=rvf.RadVelDeterminationNarrow(deltav,lag_max,crosscor,popt[1],popt[2],noise)
                vel_diff_matrix[spec_num,order]=cross_cor_data[1][1]
                ccfs[order,:]=cross_cor_data[0]
                gaus[order,:]=cross_cor_data[1]
            
            wf.CrossCorWriteFits(headers[spec_num],ccfs,gaus,gauspath)
        
        # Barycentric correction of the radial velocity.
        header=headers[spec_num]
        t_obs=Time(header['JD-mid'],format='jd',scale='utc')
    
        BC_corr=get_BC_vel(t_obs,starname=header['S-OBJECT'],lat=header['sitelat'],longi=header['sitelong'],alt=header['siteelev'])[0][0]
        BC_corrs.append(BC_corr)
        
        vel_diff_matrix[spec_num]+=BC_corr/1000
    
    
    vel_diff_matrix*=1000
    
    return vel_diff_matrix
