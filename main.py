import prompt as pt
import normalization as norm
import radial_velocity as rv
import data_import as di
import template as tmp
from errors import CalcErrors

#%% 
# To use this script, simply run this file and enter the required input.
# When using this script on a regular PC with multiple drives, this script and
# the data should never be on different drives, unless the user has plenty of
# time to burn.
# Many thanks to the developers of the barycorrpy, which was used to compute
# the barycentric corrections, whose work is detailed in
# http://iopscience.iop.org/article/10.3847/2515-5172/aaa4b7.


#%% Run the prompt/program.
program=pt.Program()
program.start()

#%% New main.
# Importing a template. By default, the Sun is simply used as a template, 
# though Arcturus is also available. If the user wishes to use another
# template, simply change the line below so that the template variable is
# equal to the normalized flux values for the template in question as well
# as adjusting the wavelength to match the new template.
wavelength,sun,arcturus,telluric=tmp.ArcturusTemplate(program)
template=sun

# Importing the data and headers from the data path, setting up paths for
# saving data, and determining the counts of each spectrum.
# Counts is simply the mean flux value at order 28 in the spectrum.
data,headers,paths,counts=di.ImportData(program)

# Takes the raw data and normalizes it. The normalized spectrum is then 
# cross-correlated with the template spectrum, and the results of this is 
# saved to data sub-folders.
# Noise: dictates whether or not whether plots will be produced for each 
# normalized spectrum. This is meant as a debugging tool and with individual or
# a small number of spectra. To do so, simply change the value of the noise 
# variable to literally anything but 'silent'.
# Method: determines which method is used. Currently the only supported
# method is 'median'. For more information about this method, see 
# normalization.py.
norm.NormalizeAllData(data,headers,program,paths,template,wavelength,noise='silent',method='median')

# Uses the cross-correlation functions produced above to determine the radial
# velocity of each order in each spectrum.
# Noise: as above.
vel_diff_matrix=rv.GetRadialVels(data,headers,program,paths,noise='silent')

# Computation of errors.
# At least I think so. I'm honestly not 100% certain that this is the correct
# method of calculating the errors.
errors,vel_differences_std,vel_differences_mean=CalcErrors(data,vel_diff_matrix,counts)
