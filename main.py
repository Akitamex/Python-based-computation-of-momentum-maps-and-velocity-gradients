import numpy as np
import glob
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import gaussian_filter

from LazVGT import VGT_QU_mv_error, avB2dx, getcoor, imfilter_gaussian
from LazCore import Gauss, Gauss_x, show_zero_moment_map
from LazCore import Proj

# Open FITS data set
hdul = fits.open("G111_CO.fits")
d2d = hdul[0].data
d,nx, ny, nv = d2d.shape
d = d2d[:, :, :, 0]

# Read header
header = hdul[0].header
# Calculate coordinate and line-of-sight velocity
RA, DEC, vz = getcoor(header, nx, ny, nv)
# Velocity resolution
dv = vz[1] - vz[0]  # unit: m/s
# Calculate noise level
noise = np.nanstd(d[:, :, 0])

# Constructing Pseudo Stokes parameters
dn = 32  # sub-block size
Ni = int(round(abs((vz[0] - vz[0]) / dv))) + 1  # initial velocity
Nf = int(round(abs((vz[-1] - vz[0]) / dv))) + 1  # end velocity
# Calculating Pseudo-Stokes parameters
Qi, Ui, Qie, Uie = VGT_QU_mv_error(d, dn, Ni, Nf, noise)
# Qie and Uie are the uncertainties.
# This function outputs Q, U cubes instead of maps.
# Pixels with noise larger than three sigma levels are blanked out.

Qa = Proj(Qi, 3)  # Project the Q cube into Q map
Ua = Proj(Ui, 3)  # Project the U cube into U map
Qa[Qa == 0] = np.nan
Ua[Ua == 0] = np.nan  # Qa = 0, Ua = 0 means intensity is less than 3 times signal to noise ratio.
ker = 4
# The smoothing width is FWHM = ker * 2.355
Qb = imfilter_gaussian(Qa, ker)  # Gaussian smoothing the map: FWHM = ker * 2.355 pixels
Ub = imfilter_gaussian(Ua, ker)
psi = 0.5 * np.arctan2(Ub, Qb)
# Rotate the polarization vectors by 90 degrees to indicate the magnetic field direction
psi += np.pi / 2

# Decreasing the resolution of psi for visualization purpose
dnn = 4
phi = avB2dx(np.cos(psi), np.sin(psi), dnn)

# Read mom0 file
hdul_map = show_zero_moment_map(hdul)
Ii = hdul_map.data

nx, ny = Ii.shape
Xd, Yd = np.meshgrid(np.arange(dnn // 2, dnn // 2 + ny, dnn),
                     np.arange(dnn // 2, dnn // 2 + nx, dnn))


phi = np.where(np.isnan(phi), 0, phi)

plt.figure()
plt.imshow(Ii, origin="lower", cmap="Greys")
cb = plt.colorbar(pad=0)
cb.ax.tick_params(labelsize=20, direction="in")
plt.clim(0, 80)
cb.set_label(label="Intensity [km/s]", size=20)
plt.xlabel("R.A.(J2000) [degree]", size=20)
plt.ylabel("Dec.(J2000) [degree]", size=20)
plt.tick_params(direction="in", labelsize=16)
plt.xticks([0, 15, 30, 45], ["349.183", "349.085", "349.987", "349.889"], size=20)
plt.yticks([0, 15, 30, 45], ["61.351", "61.448", "61.546", "61.644"], size=20)

# Check if phi has non-zero values before plotting quiver
if np.any(phi):
    plt.quiver(Yd[:-1, :-1], Xd[:-1, :-1], np.cos(phi), np.sin(phi), headwidth=0, scale=20, color="r")
    plt.quiver(Yd[:-1, :-1], Xd[:-1, :-1], -np.cos(phi), -np.sin(phi), headwidth=0, scale=20, color="r")

plt.title("vgt_final")
plt.show()

# Superimpose the gradient and magnetic field vectors on the same intensity map
# Note: Julia uses the Cartesian convention for defining angle, i.e., the angle is from right to top anticlockwise
# plt.figure()
# plt.imshow(Ii, origin="lower", cmap="Greys")
# cb = plt.colorbar(pad=0)
# cb.ax.tick_params(labelsize=20, direction="in")
# plt.clim(0, 80)
# cb.set_label(label="Intensity [km/s]", size=20)
# plt.xlabel("R.A.(J2000) [degree]", size=20)
# plt.ylabel("Dec.(J2000) [degree]", size=20)
# plt.tick_params(direction="in", labelsize=16)
# plt.xticks([0, 15, 30, 45], ["349.183", "349.085", "349.987", "349.889"], size=20)
# plt.yticks([0, 15, 30, 45], ["61.351", "61.448", "61.546", "61.644"], size=20)
# plt.quiver(Yd[:-1, :-1], Xd[:-1, :-1], np.cos(phi), np.sin(phi), headwidth=0, scale=20, color="r")
# plt.quiver(Yd[:-1, :-1], Xd[:-1, :-1], -np.cos(phi), -np.sin(phi), headwidth=0, scale=20, color="r")
# plt.title("VGT_demo")
# plt.show()

# Calculating uncertainties psie
# Qe[np.isnan(Qe)] = 0
# Ue[np.isnan(Ue)] = 0
# Qex = np.sqrt(Proj(Qe**2, 3))
# Uex = np.sqrt(Proj(Ue**2, 3))
# UQ = Ub / Qb
# UQe = np.abs(UQ) * np.sqrt((Qex / Qb)**2.0 + (Uex / Ub)**2.0)
# psie = 0.5 * UQe / (1.0 + UQ**2.0)