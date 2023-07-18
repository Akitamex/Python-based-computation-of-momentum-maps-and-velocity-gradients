import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit

from LazCore import Gauss, Gauss_x
from LazCore import margin_error

def avB2dx(Ax, Ay, dn):
    nx, ny = Ax.shape
    Ana = np.zeros((nx // dn, ny // dn))
    for j in range(ny // dn):
        for i in range(nx // dn):
            is_ = (i - 1) * dn + 1
            ie = i * dn
            js = (j - 1) * dn + 1
            je = j * dn
            Axx = Ax[is_:ie, js:je]
            Ayy = Ay[is_:ie, js:je]
            Ana[i, j] = np.arctan(np.mean(Ayy[~np.isnan(Ayy)]) / np.mean(Axx[~np.isnan(Ayy)]))
    return Ana

def VGT_QU_mv_error(ppvi, dn, Ni, Nf, noise):
    nx, ny, nv = ppvi.shape
    Q = np.zeros((nx, ny, nv))
    U = np.zeros((nx, ny, nv))
    Qe = np.zeros((nx, ny, nv))
    Ue = np.zeros((nx, ny, nv))
    for k in range(Ni, Nf):
        ds = ppvi[:, :, k]
        dsx = np.roll(ds, shift=(1, 0), axis=(0, 1)) - ds
        dsy = np.roll(ds, shift=(0, 1), axis=(0, 1)) - ds
        dna, dns = sban2d_SNR_mv(dsx, dsy, dn, ds, noise)
        for ii in range(nx):
            for jj in range(ny):
                if np.isnan(dna[ii, jj]):
                    Q[ii, jj, k] = 0
                    U[ii, jj, k] = 0
                    Qe[ii, jj, k] = 0
                    Ue[ii, jj, k] = 0
                else:
                    Q[ii, jj, k] = ds[ii, jj] * np.cos(2.0 * dna[ii, jj])
                    U[ii, jj, k] = ds[ii, jj] * np.sin(2.0 * dna[ii, jj])
                    ce = np.abs(2 * np.sin(2.0 * dna[ii, jj]) * dns[ii, jj])
                    se = np.abs(2 * np.cos(2.0 * dna[ii, jj]) * dns[ii, jj])
                    Qe[ii, jj, k] = np.abs(Q[ii, jj, k]) * np.sqrt((noise / ds[ii, jj]) ** 2 + (ce / np.cos(2.0 * dna[ii, jj])) ** 2)
                    Ue[ii, jj, k] = np.abs(U[ii, jj, k]) * np.sqrt((noise / ds[ii, jj]) ** 2 + (se / np.sin(2.0 * dna[ii, jj])) ** 2)
        print(f"{k}/{Nf}")
    return Q, U, Qe, Ue

def hist_new(data, range):
    h, ax = np.histogram(data, bins=range)
    ac = h
    return ax, ac

def sban2d_mv(Ax, Ay, dn):
    nx, ny = Ax.shape
    Ana = np.zeros((nx, ny))
    Ans = np.zeros((nx, ny))
    for j in range(dn // 2, ny - dn // 2):
        for i in range(dn // 2, nx - dn // 2):
            is_ = i - dn // 2 + 1
            ie = i + dn // 2
            js = j - dn // 2 + 1
            je = j + dn // 2
            try:
                Axx = Ax[is_:ie, js:je]
                Ayy = Ay[is_:ie, js:je]
                binsize = 100
                Apeak, Adisp = fit_gaussian_2d(Axx, Ayy, binsize)
                Ana[i, j] = Apeak
                Ans[i, j] = Adisp
            except:
                Ana[i, j] = np.nan
                Ans[i, j] = np.nan
    Ana[1:dn // 2 - 1, :] = np.nan
    Ana[:, 1:dn // 2 - 1] = np.nan
    Ana[nx - dn // 2 + 1:nx, :] = np.nan
    Ana[:, ny - dn // 2 + 1:ny] = np.nan
    Ans[np.isnan(Ana)] = np.nan
    return Ana, Ans

def getmap(ppv, vz, peak):
    nx, ny, nv = ppv.shape
    I = np.zeros((nx, ny))
    C = np.zeros((nx, ny))
    Ch = np.zeros((nx, ny))
    for i in range(nx):
        for j in range(ny):
            for k in range(nv):
                I[i, j] += ppv[i, j, k]
                C[i, j] += ppv[i, j, k] * vz[k]
    C /= I
    C[C > np.max(vz)] = np.nan
    C[C < np.min(vz)] = np.nan
    dC = np.nanstd(C[~np.isnan(C)])
    for i in range(nx):
        for j in range(ny):
            vchp = vz[peak]
            vch = vz[k]
            if abs(vch - vchp) < 0.5 * dC:
                Ch[i, j] += ppv[i, j, k]
    return I, C, Ch

def getcoor(header, nx, ny, nv):
    CRPIX1 = header["CRPIX1"]
    CRVAL1 = header["CRVAL1"]
    CDELT1 = header["CDELT1"]
    CRPIX2 = header["CRPIX2"]
    CRVAL2 = header["CRVAL2"]
    CDELT2 = header["CDELT2"]
    CRPIX3 = header["CRPIX3"]
    CRVAL3 = header["CRVAL3"]
    CDELT3 = header["CDELT3"]
    RA = np.zeros(nx)
    DEC = np.zeros(ny)
    vz = np.zeros(nv)
    for i in range(nx):
        RA[i] = CRVAL1 + (i - CRPIX1) * CDELT1
    for j in range(ny):
        DEC[j] = CRVAL2 + (j - CRPIX2) * CDELT2
    for k in range(nv):
        vz[k] = CRVAL3 + (k - CRPIX3) * CDELT3
    xi = CRVAL1 + (0 - CRPIX1) * CDELT1
    xe = CRVAL1 + (nx - CRPIX1) * CDELT1
    yi = CRVAL2 + (0 - CRPIX2) * CDELT2
    ye = CRVAL2 + (ny - CRPIX2) * CDELT2
    print(f"xi={xi}; xe={xe}; yi={yi}; ye={ye}")
    return RA, DEC, vz

def getpeak(d):
    nx, ny, nv = d.shape
    peak = np.zeros(nv)
    for k in range(nv):
        for j in range(ny):
            for i in range(nx):
                if ~np.isnan(d[i, j, k]):
                    peak[k] += d[i, j, k]
    return peak

def imfilter_gaussian(d, p):
    im = gaussian_filter(d, sigma=p)
    return im

def de_resolution(b, dn):
    nx, ny = b.shape
    bx = np.zeros((nx // dn, ny // dn))
    for i in range(nx // dn):
        for j in range(ny // dn):
            a = b[i * dn - dn + 1:i * dn, j * dn - dn + 1:j * dn]
            bx[i, j] = np.mean(a)
    return bx

def fit_gaussian_2d(Ax, Ay, binsize):
    phi = np.arctan(Ay / Ax)
    phix = phi[~np.isnan(phi)]
    ax, ac = hist_new(phix[:], np.linspace(-np.pi / 2, np.pi / 2, binsize + 1))
    ax = 0.5 * (ax[:-1] + ax[1:])
    if np.abs(ax[np.argmax(ac)])[0] < np.pi / 4:
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
        ax = ax - np.pi / 2
        ac = np.fft.fftshift(ac)
        y = ac / np.sum(ac)
        y[np.isnan(y)] = 0
        try:
            fit1, _ = curve_fit(Gauss_x, ax, y, [np.max(y), -np.pi / 2, 1.0, 0])
            if np.max(y) == 0:
                sigma = np.nan
                return fit1[1], sigma
            else:
                sigma = margin_error(fit1, 0.05)
                return fit1[1], sigma[1]
        except:
            fit1, _ = curve_fit(Gauss, ax, y, [np.max(y), -np.pi / 2, 1.0])
            if np.max(y) == 0:
                sigma = np.nan
                return fit1[1], sigma
            else:
                sigma = margin_error(fit1, 0.05)
                return fit1[1], sigma[1]

def sban2d_SNR(Ax, Ay, dn, Ch, noise):
    nx, ny = Ax.shape
    Ana = np.zeros((nx // dn, ny // dn))
    Ans = np.zeros((nx // dn, ny // dn))
    Ax[Ch < 3 * noise] = np.nan
    Ay[Ch < 3 * noise] = np.nan
    for j in range(1, ny // dn):
        for i in range(1, nx // dn):
            is_ = (i - 1) * dn + 1
            ie = i * dn
            js = (j - 1) * dn + 1
            je = j * dn
            Axx = Ax[is_:ie, js:je]
            Ayy = Ay[is_:ie, js:je]
            Axx = Axx[~np.isnan(Axx)]
            Ayy = Ayy[~np.isnan(Ayy)]
            if len(Axx) > 100:
                binsize = 100
                Apeak, Adisp = fit_gaussian_2d(Axx, Ayy, binsize)
                Ana[i, j] = Apeak
                Ans[i, j] = Adisp
            else:
                Ana[i, j] = np.nan
                Ans[i, j] = np.nan
    return Ana, Ans

def sban2d_SNR_mv(Ax, Ay, dn, Ch, noise):
    nx, ny = Ax.shape
    Ana = np.zeros((nx, ny))
    Ans = np.zeros((nx, ny))
    Ax[Ch < 3 * noise] = np.nan
    Ay[Ch < 3 * noise] = np.nan
    for j in range(dn // 2, ny - dn // 2):
        for i in range(dn // 2, nx - dn // 2):
            is_ = i - dn // 2 + 1
            ie = i + dn // 2
            js = j - dn // 2 + 1
            je = j + dn // 2
            try:
                Axx = Ax[is_:ie, js:je]
                Ayy = Ay[is_:ie, js:je]
                Axx = Axx[~np.isnan(Axx)]
                Ayy = Ayy[~np.isnan(Ayy)]
                if len(Axx) > 100:
                    binsize = 100
                    Apeak, Adisp = fit_gaussian_2d(Axx, Ayy, binsize)
                    Ana[i, j] = Apeak
                    Ans[i, j] = Adisp
                else:
                    Ana[i, j] = np.nan
                    Ans[i, j] = np.nan
            except:
                Ana[i, j] = np.nan
                Ans[i, j] = np.nan
    Ana[1:dn // 2 - 1, :] = np.nan
    Ana[:, 1:dn // 2 - 1] = np.nan
    Ana[nx - dn // 2 + 1:nx, :] = np.nan
    Ana[:, ny - dn // 2 + 1:ny] = np.nan
    Ans[np.isnan(Ana)] = np.nan
    return Ana, Ans

def VGT_QU(ppvi, dn, Ni, Nf, noise):
    nx, ny, nv = ppvi.shape
    nxx, nyy = nx // dn, ny // dn
    Q = np.zeros((nxx, nyy, nv))
    U = np.zeros((nxx, nyy, nv))
    Q2d = np.zeros((nxx, nyy))
    U2d = np.zeros((nxx, nyy))
    for k in range(Ni, Nf):
        ds = ppvi[:, :, k]
        dsx = np.roll(ds, shift=(1, 0), axis=(0, 1)) - ds
        dsy = np.roll(ds, shift=(0, 1), axis=(0, 1)) - ds
        dna, dns = sban2d_SNR(dsx, dsy, dn, ds, noise)
        Ia = de_resolution(ppvi[:, :, k], dn)
        for ii in range(nxx):
            for jj in range(nyy):
                if np.isnan(dna[ii, jj]):
                    Q[ii, jj, k] = 0
                    U[ii, jj, k] = 0
                else:
                    Q[ii, jj, k] = Ia[ii, jj] * np.cos(2 * dna[ii, jj])
                    U[ii, jj, k] = Ia[ii, jj] * np.sin(2 * dna[ii, jj])
        Q2d += Q[:, :, k]
        U2d += U[:, :, k]
    return Q2d, U2d

def VGT_QU_mv(ppvi, dn, Ni, Nf, noise):
    nx, ny, nv = ppvi.shape
    Q = np.zeros((nx, ny, nv))
    U = np.zeros((nx, ny, nv))
    Q2d = np.zeros((nx, ny))
    U2d = np.zeros((nx, ny))
    for k in range(Ni, Nf):
        ds = ppvi[:, :, k]
        dsx = np.roll(ds, shift=(1, 0), axis=(0, 1)) - ds
        dsy = np.roll(ds, shift=(0, 1), axis=(0, 1)) - ds
        dna, dns = sban2d_SNR_mv(dsx, dsy, dn, ds, noise)
        for ii in range(nx):
            for jj in range(ny):
                if np.isnan(dna[ii, jj]):
                    Q[ii, jj, k] = 0
                    U[ii, jj, k] = 0
                else:
                    Q[ii, jj, k] = ds[ii, jj] * np.cos(2.0 * dna[ii, jj])
                    U[ii, jj, k] = ds[ii, jj] * np.sin(2.0 * dna[ii, jj])
        Q2d += Q[:, :, k]
        U2d += U[:, :, k]
    return Q2d, U2d

def linspace(a, b, c):
    width = (b - a) / c
    x = np.zeros(c)
    xx = np.linspace(a, b, num=c)
    for i in range(c):
        x[i] = xx[i]
    return x
