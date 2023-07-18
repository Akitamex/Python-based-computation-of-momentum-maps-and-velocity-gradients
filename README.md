# Python-based computation of momentum maps and velocity gradients.pdf

	Context. 
Students, professors, and scientists who work in astrophysics actively use .fits files, which show maps of the star clouds,
for their own purposes. However, the calculation of 0-2 moment maps from these files is a difficult procedure which requires a lot
of time to calculate and an understanding of the particular features of the .fits file they use. Some .fits files use non-SI units, which
makes it difficult to work with them, especially when calculating 0-2 momentum maps. The calculation of velocity gradients is a
complex task with a large number of calculations that require the work of all the factors described. At the moment there are good
solutions for calculating velocity gradients in Julia programming language, but the task was set to calculate momentum maps 0-2
and velocity gradients in Python.
	Aims. 
Develop code in Python that will adapt or will ask to help with adapting the .fits file to the SI system, correctly and quickly
calculate moment maps 0-2 and velocity gradients in a given particular dataset and investigate their consequences for the observed
system’s dynamics and magnetic field configuration.
	Methods. 
This work uses Python 3, astropy.io.fits, numpy, and matplotlib.pyplot libraries, which are used to read .fits files,
process information in 3D ranges, and calculate totals to produce 0-2 moment maps and velocity gradients map. Also the study made
use of a 3D PPV cube generated from FITS data to illustrate the Position-Position-Velocity distribution. The VGT QU mv error
function, which takes the 3D PPV cube as input and generates the pseudo Stokes parameters for a particular velocity range, was
used to derive the pseudo Stokes parameters (Qi, Ui).Using the Proj function, the Qi and Ui cubes were projected onto 2D maps (Qa,
Ua). To improve visualization, the imfilter gaussian function was used to apply Gaussian smoothing to the maps. The polarization
angle (psi) was calculated by applying the atan function to the ratio of Ui and Qi, followed by a 90-degree rotation to identify
the direction of the magnetic field. Using matplotlib’s quiver function, the intensity map (Ii) from another FITS file was read and
overlaid with gradient vectors denoting the polarization angle.
	Results. 
The result of this work is the finished code, which can be used on a web site or desktop application to calculate moments
0-2 and velocity gradients, entering data about the structure of the file and its features, and getting pictures of moment maps 0-
2. The ”Velocity Gradients” that plot shows the observed system’s intensity map with superimposed gradient vectors indicating
the direction of velocity gradients. The figure provides useful visual insights into the dataset’s velocity structure and polarization
features.
	Conclusion. 
A unified calculation of 0-2 moment maps is not an easy task based on the complexity of adapting the code to all
variants of .fits files, because these files are used with different structures for different needs of researchers. The different structure
refers to additional information in the file header, file data type, file size and dimensions, which can also be confusing to the person
reading the file and calculating 0-2 momentum maps and velocity gradients. The analysis of velocity gradients and polarization
qualities contributes to a better understanding of the observed system’s dynamics and magnetic field structure. The findings of this
study show the importance of velocity gradients in investigating the underlying physical processes and their impact on the observed
phenomena.

# Calculations and formulas are used to obtain Momentum Maps 0,1 and 2.

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

Compute the moment 0 map
moment0 = np.sum(data_cube, axis=0)

Compute the moment 1 map
moment1 = np.sum(data_cube * vel_axis[:, None, None], axis=0) / np.sum(data_cube, axis=0)

Compute the moment 2 map
mean_velocity = np.sum(data_cube * vel_axis[:, None, None], axis=0) / np.sum(data_cube, axis=0)
moment2 = np.sqrt(np.sum(data_cube * (vel_axis[:, None, None] - mean_velocity)**2, axis=0) / np.sum(data_cube, axis=0))

# Calculations and formulas are used to obtain Velocity Gradients.

VGT\_QU\_mv\_error: Although not explicitly visible in the code, this function is used to calculate the pseudo Stokes parameters (Qi and Ui) as well as their uncertainty (Qie and Uie). These factors represent the data's polarization qualities. The particular implementation of VGT\_QU\_mv\_error is not accessible, although it is likely to entail some type of statistical analysis and cube manipulation (d).

The Proj function converts a 3D data cube into a 2D map along a specified axis. It is utilized in the code to project the pseudo Stokes parameters (Qi and Ui) into their corresponding maps (Qa and Ua). This projection enables for the visualization and investigation of polarization properties in 2D.

imfilter\_gaussian: This function smooths the input map using a Gaussian smoothing filter. It is used in the code to smooth the projected maps (Qa and Ua) with a smoothing width (ker). Smoothing the maps reduces noise and improves feature visibility.

The atan. function computes the arctangent of two arrays element by element. It is utilized in the code to determine the polarization angle (psi) by taking the arctangent of the smoothed map ratios (Ub and Qb). The atan. function aids in determining the polarization vectors' direction.

The programme computes coordinate values (RA, DEC) and line-of-sight velocities (vz). Processing the FITS header information (header) to determine the spatial and velocity dimensions of the data is involved in these calculations. These position and velocity estimates are used to contextualize and label the displays.

imshow and quiver are matplotlib.pyplot library functions for viewing the intensity map and superimposing the gradient vectors, respectively. They allow the "Velocity Gradients" plot to be created by presenting the data and overlaying the gradient vectors on top of it.

# Computation of velocity gradients

To open the FITS file containing the intensity map data, use a Python FITS package such as astropy.io.fits.

Extract the intensity map from the FITS file and carry out any necessary preprocessing. These may involve resolving missing data, noise removal, or any other data adjustments that are required.

To compute the gradient vectors, use an appropriate gradient computation method on the intensity map. The avB2dx function is used in the given code to calculate the gradient vectors, which are stored in the phi variable.

To construct a coordinate grid with the dimensions and spacing appropriate to the resolution of the gradient vectors, use the np.meshgrid function.

To view the intensity map, use the imshow function from matplotlib.pyplot. Customize the color map, colorbar, and axis names according to your study needs.

To superimpose the gradient vectors onto the intensity map, use the quiver function from matplotlib.pyplot. Using the coordinate grid and the derived gradient vectors, specify the vector coordinates and directions. Adjust the appearance of the vectors to your liking, including variables such as headwidth, scale, and color.
\par
Using the title function from matplotlib.pyplot, add a descriptive title to the plot.

To display the plot on your screen for visualization and analysis, use the show function from matplotlib.pyplot.

# Improved description of how Velocity Gradients are calculated

VGT\_QU\_mv\_error(d, dn, Ni, Nf, noise): The Pseudo Stokes parameters, Q and U, as well as their uncertainties, Qie and Uie, are calculated using this function. The 3D PPV cube d, sub-block size n, initial velocity index Ni, final velocity index Nf, and noise level noise are all input parameters. The calculations required to determine the pseudo Stokes parameters and associated uncertainties are completed.

The function Proj(Qi, 3) transforms the Q cube created in the previous step into a Q map. It conducts the projection to create the Q map using the Q cube Qi as input.

Imfilter\_gaussian(Qa, ker): The input map Qa is smoothed using a Gaussian algorithm using this function. For the smoothing procedure, a kernel with a certain FWHM specified by ker is used. The final smoothed map is given back.

avB2dx(cos.(psi), sin.(psi), and dnn) The gradient vectors are calculated using the input polarization angles as psi by this function. The gradient vectors are computed using the cosine and sine components of psi and the relevant mathematical calculations. The resolution of the generated gradient vectors is reduced for display reasons by downsampling them depending on the factor dnn.


