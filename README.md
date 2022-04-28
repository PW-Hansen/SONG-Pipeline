# SONG-Pipeline
A pipeline created to process raw data from the SONG telescope.

Files and their contents:

main.py:
As the name indicates, this is the primary file for the Pipeline. When this
file is run, it will bring up a prompt asking where the data that one wishes
to reduce is located.

ARCTURUS_atlas:
Contains normalized flux values for the Sun and Arcturus, with matching
wavelenghts. Used as the default template.

data_import.py:
The function which takes the datapath provided in the main file and
returns the data.

errors.py:
The error calculation function.

helpful_functions.py:
As the name suggests, this file contains a variety of functions that are used
either by other files but kept separate for the sake of neatness, as well as
a few functions for data visualization (PlotOrder and PlotSpectrum).

normalization.py:
The functions for normalizing the dataset. The first function runs through the
data and checks whether the .fits function currently being processed has 
already been normalized, and if not, it calls goes through the entire
spectrum order-by-order using the second function, then saves the results.
This second function calls a normalization method - which are found at the
of the normalization.py file - and uses that to provide the first function
with the normalized values for an order.
The remaining function is the actual method used for normalizing the data,
but users who wish to utilize this pipeline can use their own function instead.

prompts.py:
Prompt which is brought up when running the file. Asks the user to supply
the following:
- Location of the settings file if not found in the dictionary.
- Whether any of the settings should be changed.
- Where the data to be processed is located.

radial_velocity_functions.py:
The two functions which are used to determine the radial velocity by locating
the peak in the cross-correlation function. The first takes a broader look and
gives a rough estimate for the location of the peak, while the second takes 
a much narrower look and provides a more precise result.

radial_velocity.py:
Runs through the data, applying the functions in radial_velocity_functions
to them, applies barycentric correction, and so on, in order to determine
the radial velocity for each order in each spectrum.

robust.py:
See actual file for more information, as I did not write this but merely 
imported it.

template.py:
Contains one function to load the Arcturus Atlas and another largely finished
but not wholly implemented one to produce empirical templates.

write_functions.py:
Functions that writes normalization, cross-correlation functions, and 
Gaussians fit to said function to subfolders created in the datapath.
