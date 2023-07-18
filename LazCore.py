# LazCore.py

from copy import deepcopy
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import fits


def show_zero_moment_map(hdul):
        
        left_limit = 104
        right_limit = 117

        data, header = fits.getdata("G111_CO.fits",header = True)
        data = data[0, left_limit:right_limit, :, :]

        fits_file = hdul

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

        return moment_0_map

def Gauss(x, p):
    # Gaussian function
    return p[0] * np.exp(-(x - p[1])**2 * p[2])

def Gauss_x(x, p):
    # Gaussian function with additional offset
    return p[0] * np.exp(-(x - p[1])**2 * p[2]) + p[3]

def dotproductangle(a, b):
    a1 = np.cos(a)
    a2 = np.sin(a)
    b1 = np.cos(b)
    b2 = np.sin(b)
    cab = (a1 * b1 + a2 * b2) / (np.linalg.norm([a1, a2]) * np.linalg.norm([b1, b2]))
    cab[cab > 1] = 1
    cab[cab < -1] = -1
    ab = np.arccos(cab)
    return ab

def bitwise_filter(A, threshold):
    A = A - threshold
    B = (A + np.abs(A)) / (2.0 * np.abs(A))
    B[np.isnan(B)] = 0
    return B

def dot_product_3d(A, B):
    Am = np.mean(A)
    Bm = np.mean(B)
    ratio = Am / Bm
    Ax = bitwise_filter(A, Am)
    Bx = bitwise_filter(B, Bm)
    return np.sum(Ax * Bx) / len(Am)

def hist_new(data, range):
    hist, edges = np.histogram(data, bins=range)
    return edges[:-1], hist

def maxid(ax):
    return np.where(ax == np.max(ax))

def meshgrid(X, Y):
    return np.meshgrid(X, Y)

def atan2(X, Y):
    return np.arctan2(Y, X)

def fit_gaussian_2d(Ax, Ay, binsize):
    phi = np.arctan(Ay / Ax)
    phix = phi[~np.isnan(phi)]
    ax, ac = hist_new(phi[:], np.linspace(-np.pi/2, np.pi/2, binsize+1))
    ax = 0.5 * (ax[:-1] + ax[1:])
    if np.abs(ax[maxid(ac)][0]) < np.pi/4:
        y = ac / np.sum(ac)
        y[np.isnan(y)] = 0
        try:
            fit1, _ = curve_fit(Gauss_x, ax, y, [np.max(y), 0.0, 1.0, 0])
            if np.max(y) == 0:
                sigma = np.nan
                return fit1[1], sigma
            else:
                sigma = margin_error(fit1, 0.05)
                return fit1[1], sigma[1]
        except:
            fit1, _ = curve_fit(Gauss, ax, y, [np.max(y), 0.0, 1.0])
            if np.max(y) == 0:
                sigma = np.nan
                return fit1[1], sigma
            else:
                sigma = margin_error(fit1, 0.05)
                return fit1[1], sigma[1]
    else:
        ax = ax - np.pi/2
        ac = np.fft.fftshift(ac)
        y = ac / np.sum(ac)
        y[np.isnan(y)] = 0
        try:
            fit1, _ = curve_fit(Gauss_x, ax, y, [np.max(y), -np.pi/2, 1.0, 0])
            if np.max(y) == 0:
                sigma = np.nan
                return fit1[1], sigma
            else:
                sigma = margin_error(fit1, 0.05)
                return fit1[1], sigma[1]
        except:
            fit1, _ = curve_fit(Gauss, ax, y, [np.max(y), -np.pi/2, 1.0])
            if np.max(y) == 0:
                sigma = np.nan
                return fit1[1], sigma
            else:
                sigma = margin_error(fit1, 0.05)
                return fit1[1], sigma[1]

def margin_error(fit, alpha):
    alpha /= 2.0
    hessian = np.linalg.inv(fit[1])
    error = np.sqrt(np.diag(hessian)) * np.sqrt(fit[1].shape[0]) * alpha
    return error

def sban2d(Ax, Ay, dn):
    nx, ny = Ax.shape
    Ana = np.zeros((nx // dn, ny // dn))
    Ans = np.zeros((nx // dn, ny // dn))
    for j in range(0, ny, dn):
        for i in range(0, nx, dn):
            is_ = i
            ie = i + dn
            js_ = j
            je = j + dn
            Axx = Ax[is_:ie, js_:je]
            Ayy = Ay[is_:ie, js_:je]
            binsize = 100
            Apeak, Adisp = fit_gaussian_2d(Axx, Ayy, binsize)
            Ana[i // dn, j // dn] = Apeak
            Ans[i // dn, j // dn] = Adisp
    return Ana, Ans

def sobel_conv_2d(A):
    Kx, Ky = sobel_kernel_2d(A)
    Ax = convoluting_kernel(A, Kx)
    Ay = convoluting_kernel(A, Ky)
    return Ax, Ay

def sobel_conv_3d(A):
    Kx, Ky, Kz = sobel_kernel_3d(A)
    Ax = convoluting_kernel(A, Kx)
    Ay = convoluting_kernel(A, Ky)
    Az = convoluting_kernel(A, Kz)
    return Ax, Ay, Az

def sobel_kernel_2d(A):
    Ax = np.zeros_like(A)
    Ay = np.zeros_like(A)
    vp = sobel_parallel(3)
    vl = sobel_perpendicular(3)
    Axx = np.zeros((3, 3))
    Ayy = np.zeros((3, 3))
    for j in range(3):
        for i in range(3):
            Axx[i, j] = vp[i] * vl[j]
            Ayy[i, j] = vl[i] * vp[j]
    Ax[:3, :3] = np.roll(Axx, (1, 1), axis=(0, 1))
    Ay[:3, :3] = np.roll(Ayy, (1, 1), axis=(0, 1))
    Ax = np.roll(Ax, (-1, -1), axis=(0, 1))
    Ay = np.roll(Ay, (-1, -1), axis=(0, 1))
    return Ax, Ay

def sobel_kernel_3d(A):
    Ax = np.zeros_like(A)
    Ay = np.zeros_like(A)
    Az = np.zeros_like(A)
    vp = sobel_parallel(3)
    vl = sobel_perpendicular(3)
    Axx = np.zeros((3, 3, 3))
    Ayy = np.zeros((3, 3, 3))
    Azz = np.zeros((3, 3, 3))
    for k in range(3):
        for j in range(3):
            for i in range(3):
                Axx[i, j, k] = vp[i] * vl[j] * vl[k]
                Ayy[i, j, k] = vl[i] * vp[j] * vl[k]
                Azz[i, j, k] = vl[i] * vl[j] * vp[k]
    Ax[:3, :3, :3] = np.roll(Axx, (1, 1, 1), axis=(0, 1, 2))
    Ay[:3, :3, :3] = np.roll(Ayy, (1, 1, 1), axis=(0, 1, 2))
    Az[:3, :3, :3] = np.roll(Azz, (1, 1, 1), axis=(0, 1, 2))
    Ax = np.roll(Ax, (-1, -1, -1), axis=(0, 1, 2))
    Ay = np.roll(Ay, (-1, -1, -1), axis=(0, 1, 2))
    Az = np.roll(Az, (-1, -1, -1), axis=(0, 1, 2))
    return Ax, Ay, Az

def sobel_parallel(size):
    v = np.zeros(size)
    if size >= 3:
        v[1] = -1
        v[-1] = 1
    else:
        raise ValueError("Size < 3 not supported")
    return v

def sobel_perpendicular(size):
    v = np.zeros(size)
    if size >= 3:
        v[1] = 1
        v[0] = 2
        v[-1] = 1
    else:
        raise ValueError("Size < 3 not supported")
    return v

def convoluting_kernel(A, B):
    Af = np.fft.fftn(A)
    Bf = np.fft.fftn(B)
    Cf = Af * Bf
    C = np.real(np.fft.ifftn(Cf))
    return C

def ppv(d, v, binnum):
    nx, ny, nz = d.shape
    offset = 1e-9
    bindiff = (np.max(v) - np.min(v)) / binnum
    minv = np.min(v)
    p = np.zeros((nx, ny, binnum))
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                vb = np.round((v[i, j, k] - minv) / bindiff).astype(int) + 1
                if vb > binnum:
                    vb = binnum
                p[i, j, vb-1] += d[i, j, k]
    return p

def Proj(d, dim):
    nx, ny, nz = d.shape
    if dim == 1:
        dx = np.zeros((ny, nz))
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    dx[j, k] += d[i, j, k]
    elif dim == 2:
        dx = np.zeros((nx, nz))
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    dx[i, k] += d[i, j, k]
    elif dim == 3:
        dx = np.zeros((nx, ny))
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    dx[i, j] += d[i, j, k]
    return dx

def AM(a, b):
    z = a * b
    ab = np.where(~np.isnan(z))
    ca = np.cos(a[ab])
    sa = np.sin(a[ab])
    cb = np.cos(b[ab])
    sb = np.sin(b[ab])
    c = ca * cb + sa * sb
    return np.mean(2 * c**2 - 1)
