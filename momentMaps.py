import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
from numpy import ma
from copy import deepcopy

class MomentMaps:
    def __init__(self, file_path, left_limit, right_limit):
        self.file_path = file_path
        self.left = left_limit
        self.right = right_limit

    def get_fits_file(self):
        return fits.open(self.file_path)

    def get_fits_data(self):
        return fits.getdata(self.file_path, header=True)

    def show_zero_moment_map(self):
        data, header = self.get_fits_data()
        data = data[0, self.left:self.right, :, :]

        fits_file = self.get_fits_file()

        vrad_dim, dec_dim, ra_dim = data.shape
        ra_values = np.arange(ra_dim)
        dec_values = np.arange(dec_dim)

        vrad_reference = 5.62963
        vrad_increment = fits_file[0].header['CDELT3']  # Increment value for VRAD

        vrad_values = [[np.zeros(ra_dim) for _ in range(ra_dim)] for _ in range(dec_dim)]
        vrad_last = []

        for vrad in range(13):
            for dec in range(dec_dim):
                for ra in range(ra_dim):
                    vrad_values[dec][ra] = vrad_reference
            vrad_last.append(deepcopy(vrad_values))
            vrad_reference += vrad_increment/1000


        moment_0_map = np.sum(data*vrad_last, axis=0)/13

        plt.imshow(moment_0_map, origin='lower', extent=(ra_values[0], ra_values[-1], dec_values[0], dec_values[-1]),
                   vmin=-1, vmax=6, cmap="jet")

        plt.colorbar(label='Spectral line intensity, km/s')
        plt.xlabel('x-label')
        plt.ylabel('y-label')
        plt.title('0 Moment Map')
        plt.show()
        return moment_0_map

    def show_velocity_map(self):
        data, header = self.get_fits_data()
        data = data[0, self.left:self.right, :, :]

        fits_file = self.get_fits_file()

        vrad_dim, dec_dim, ra_dim = data.shape

        data = ma.masked_where(data < 0.45, data)

        vrad_increment = fits_file[0].header['CDELT3']
        vrad_reference = 5.62963

        vrad_values = [[np.zeros(ra_dim) for _ in range(ra_dim)] for _ in range(dec_dim)]
        vrad_last = []
        ra_values = np.arange(ra_dim)
        dec_values = np.arange(dec_dim)

        for vrad in range(13):
            for dec in range(dec_dim):
                for ra in range(ra_dim):
                    vrad_values[dec][ra] = vrad_reference
            vrad_last.append(deepcopy(vrad_values))
            vrad_reference+=vrad_increment/1000

        moment_0_map = np.sum(data * vrad_last, axis=0) / 13
        moment_1_map = np.sum(data*vrad_last * vrad_last, axis=0) / (moment_0_map*13)

        plt.imshow(moment_1_map, origin='lower', extent=(ra_values[0], ra_values[-1], dec_values[0], dec_values[-1]),
                    cmap='gist_ncar', vmin=-1, vmax=6)

        plt.colorbar(label='Velocity, km/s')
        plt.xlabel('x-label')
        plt.ylabel('y-label')
        plt.title('1 Moment Map')
        plt.show()
        return moment_1_map

    def show_second_moment_map(self):
        data, header = self.get_fits_data()
        data = data[0, self.left:self.right, :, :]

        fits_file = self.get_fits_file()

        vrad_dim, dec_dim, ra_dim = data.shape

        ra_values = np.arange(ra_dim)
        dec_values = np.arange(dec_dim)

        data = ma.masked_where(data < 0.45, data)

        vrad_reference = 5.62963
        vrad_increment = fits_file[0].header['CDELT3']  # Increment value for VRAD

        vrad_values = [[np.zeros(ra_dim) for _ in range(ra_dim)] for _ in range(dec_dim)]
        vrad_last = []

        for vrad in range(13):
            for dec in range(dec_dim):
                for ra in range(ra_dim):
                    vrad_values[dec][ra] = vrad_reference
            vrad_last.append(deepcopy(vrad_values))
            vrad_reference += vrad_increment/1000

        moment_0_map = np.sum(data * vrad_last, axis=0) / 13
        moment_1_map = np.sum(data * vrad_last * vrad_last, axis=0) / (moment_0_map * 13)
        moment_2_map = np.sum(data * vrad_last * pow(vrad_last-moment_1_map, 2), axis=0) / (moment_0_map*13)

        plt.imshow(moment_2_map, origin='lower', extent=(ra_values[0], ra_values[-1], dec_values[0], dec_values[-1]),
                   )
        plt.colorbar(label='Velocity dispersion, km²/s²')
        plt.xlabel('x-label')
        plt.ylabel('y-label')
        plt.title('2 Moment Map')
        plt.show()


def main():
    file_path = '/Users/stormwind/Desktop/Fits/G82.65-13CO.fits'
    maps = MomentMaps(file_path, 104, 117)
    maps.show_zero_moment_map()
    maps.show_velocity_map()
    maps.show_second_moment_map()


if __name__ == "__main__":
    main()