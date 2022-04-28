import numpy as np
import robust

#%% Errors.
def CalcErrors(data,vel_diff_matrix,counts):
    vel_differences_matrix=np.zeros((len(data),len(data[0][0])))
    
    # Comparing other spectra to the one with the highest count.
    spec_num=np.where(max(counts)==counts)[0][0]
    
    spec_subtract=[robust.mean(vel_diff_matrix[:,order]) for order in range(len(data[0][0]))]
    
    for i in range(len(data)):
        vel_differences_matrix[i]=spec_subtract-vel_diff_matrix[i]-robust.mean(spec_subtract-vel_diff_matrix[i])
    vel_differences_matrix[spec_num]=vel_diff_matrix[spec_num]-robust.mean(vel_diff_matrix[spec_num])
    
    vel_differences_mean=[robust.mean(vel_differences_matrix[:,order]) for order in range(len(data[0][0]))]
    vel_differences_std=[robust.std(vel_differences_matrix[:,order]) for order in range(len(data[0][0]))]
    
    vel_difference_weights=[val**-2 for val in vel_differences_std]
    
    errors=[1/np.sqrt(sum(vel_difference_weights))]*len(data)
    
    return errors,vel_differences_std,vel_differences_mean