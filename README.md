Python-based computation of momentum maps and velocity gradients.pdf

Context. Students, professors, and scientists who work in astrophysics actively use .fits files, which show maps of the star clouds,
for their own purposes. However, the calculation of 0-2 moment maps from these files is a difficult procedure which requires a lot
of time to calculate and an understanding of the particular features of the .fits file they use. Some .fits files use non-SI units, which
makes it difficult to work with them, especially when calculating 0-2 momentum maps. The calculation of velocity gradients is a
complex task with a large number of calculations that require the work of all the factors described. At the moment there are good
solutions for calculating velocity gradients in Julia programming language, but the task was set to calculate momentum maps 0-2
and velocity gradients in Python.
Aims. Develop code in Python that will adapt or will ask to help with adapting the .fits file to the SI system, correctly and quickly
calculate moment maps 0-2 and velocity gradients in a given particular dataset and investigate their consequences for the observed
system’s dynamics and magnetic field configuration.
Methods. This work uses Python 3, astropy.io.fits, numpy, and matplotlib.pyplot libraries, which are used to read .fits files,
process information in 3D ranges, and calculate totals to produce 0-2 moment maps and velocity gradients map. Also the study made
use of a 3D PPV cube generated from FITS data to illustrate the Position-Position-Velocity distribution. The VGT QU mv error
function, which takes the 3D PPV cube as input and generates the pseudo Stokes parameters for a particular velocity range, was
used to derive the pseudo Stokes parameters (Qi, Ui).Using the Proj function, the Qi and Ui cubes were projected onto 2D maps (Qa,
Ua). To improve visualization, the imfilter gaussian function was used to apply Gaussian smoothing to the maps. The polarization
angle (psi) was calculated by applying the atan function to the ratio of Ui and Qi, followed by a 90-degree rotation to identify
the direction of the magnetic field. Using matplotlib’s quiver function, the intensity map (Ii) from another FITS file was read and
overlaid with gradient vectors denoting the polarization angle.
Results. The result of this work is the finished code, which can be used on a web site or desktop application to calculate moments
0-2 and velocity gradients, entering data about the structure of the file and its features, and getting pictures of moment maps 0-
2. The ”Velocity Gradients” that plot shows the observed system’s intensity map with superimposed gradient vectors indicating
the direction of velocity gradients. The figure provides useful visual insights into the dataset’s velocity structure and polarization
features.
Conclusion. A unified calculation of 0-2 moment maps is not an easy task based on the complexity of adapting the code to all
variants of .fits files, because these files are used with different structures for different needs of researchers. The different structure
refers to additional information in the file header, file data type, file size and dimensions, which can also be confusing to the person
reading the file and calculating 0-2 momentum maps and velocity gradients. The analysis of velocity gradients and polarization
qualities contributes to a better understanding of the observed system’s dynamics and magnetic field structure. The findings of this
study show the importance of velocity gradients in investigating the underlying physical processes and their impact on the observed
phenomena.

	The MomentMaps class was created to visualise moment maps 0, 1 and 2. This class has methods:
1) Constructor, which takes the file path as an argument and stores it in an object of the class
2) get_fits_file method, which opens the file using the astropy library
3) get_fits_data method, which takes a 3D array from the file fits using astropy.fits.get_data.
Next, the matplotlib library was used to create the maps themselves:
In order to calculate Momentum map 0 next steps should be used:
- 3D data array should be used
- Sizes of all arrays should be contained in the fields vrad_dim, ra_dim and dec_dim
- Arrays with data from 1 to the size of arrays ra and dec to draw the sides of the map should be created
- Then, moment_0_map by summing up intensity by axis=0, i.e. vrad values could be calculated
- Also, Gaussian filter could be used, if necessary
- Using the imshow method in matplotlib.pyplot, we have output the moment map 0, bounded by the values dec (1 to 141) and ra (1 to 168)
- The colorbar method in matplotlib.pyplot allows you to add colour intervals to the map
- Name the axes with the x-label and y-label methods and display it with the show() method

Momentum map 1: (from line 40)
- Get 3D data array and also add dimensions like in the torque map 0
- Create a new 2D array output_data to be filled with the average intensity values
- Iterate through the intensity values relative to velocity vrad and retrieve the average using the numpy mean() function
- We calculate moment_1_map, by the sum of mean intensities versus velocity, and divide by the moment map 0
- Turn on the Gaussian filter, if you like.
- Using the imshow method in matplotlib.pyplot, display the moment map 1, bounded by the values dec (1 to 141) and ra (1 to 168)
- Add colour bars using the colorbar method in matplotlib.pyplot to the map
- Name the axes with the xlabel and ylabel methods and display it with show()

Momentum map 2: (from 66)
- Read moment map 1
- Calculate the moment_2_map by multiplying the intensity by the square of the difference of the sum of the average intensity from speed and moment 1, and divide by the moment map 0
- If desired, apply a Gaussian filter and output the map using the same principle as for moments 0 and 1

To get these maps you need:
- Create an object of class MomentMap and specify path to file fits in arguments.
- Call required methods 

show_zero_moment_map(v_min, v_max, gauss_filter) , where v_min is the lower intensity limit, v_max the upper limit and gauss_filter is the boolean value that determines if smoothing is enabled
show_velocity_map(gauss_filter)
show_two_moment_map()

Everything is already done in the main() function and you only need to install the required packages with a command:
pip install -r requirements.txt

For better understanding, it could be simplified to:

# Compute the moment 0 map
moment0 = np.sum(data_cube, axis=0)

# Compute the moment 1 map
moment1 = np.sum(data_cube * vel_axis[:, None, None], axis=0) / np.sum(data_cube, axis=0)

# Compute the moment 2 map
mean_velocity = np.sum(data_cube * vel_axis[:, None, None], axis=0) / np.sum(data_cube, axis=0)
moment2 = np.sqrt(np.sum(data_cube * (vel_axis[:, None, None] - mean_velocity)**2, axis=0) / np.sum(data_cube, axis=0))




