import os
import astropy.io.fits as pf
import numpy as np

#%% Grabbing the data.
def ImportData(program):

    spectrum_cutoff=float(program.settings['spectrum_cutoff'])
    datapath=program.data_path
    
    os.chdir(datapath)
    
    # Creating folders in which data will later be saved.
    sub_folders=['ccfs','gaus','norm']
    for sub_folder in sub_folders:
        if sub_folder not in os.listdir():
            os.mkdir(sub_folder)
        
    if '\\' in datapath:
        normpath=datapath+'\\norm'
        gauspath=datapath+'\\gaus'    
        ccfspath=datapath+'\\ccfs'
    else:
        normpath=datapath+'/norm'
        gauspath=datapath+'/gaus'
        ccfspath=datapath+'/ccfs'
    
    paths=[normpath,gauspath,ccfspath]
    
    # Importing data.
    headers=[]
    data=[]
    files=os.listdir()
    files.sort()
        
    for file in files:
        print(file)
        if '.fits' in file:
            with pf.open(file) as hdul:
                # Quick and dirty filter: demands that there be a certain number of average
                # counts at order 30.
                new_data=hdul[0].data
                header=hdul[0].header
                JDtime=header['JD-DATE']
                # Chopping off edges of the data, as well as adjusting it depending on
                # whether or not it was taken before or after the change on Teide.
                if JDtime>2458210.00000:
                    if np.mean(new_data[0][30])>spectrum_cutoff:
                        if np.mean(new_data[0][30]/new_data[2][30])<1:
                            data.append(new_data[:,2:51,:])
                            headers.append(header)
                elif JDtime<2458198.00000:
                    if np.mean(new_data[0][29])>spectrum_cutoff:
                        if np.mean(new_data[0][29]/new_data[2][29])<1:
                            data.append(new_data[:,1:50,:])
                            headers.append(header)
    
    counts=[]
    for spec_num in range(len(data)):
        counts.append(np.mean(data[spec_num][0][28]))
        
    return data,headers,paths,counts
