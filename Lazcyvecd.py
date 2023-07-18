# LazCyvecd.py

import numpy as np

def cgrad(v, dx):
    vdx = np.zeros_like(v)
    vdy = np.zeros_like(v)
    vdz = np.zeros_like(v)
    nx, ny, nz = v.shape
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                ip = (i + 1) % nx
                jp = (j + 1) % ny
                kp = (k + 1) % nz
                vdx[i, j, k] = (v[ip, j, k] - v[i, j, k]) / dx
                vdy[i, j, k] = (v[i, jp, k] - v[i, j, k]) / dx
                vdz[i, j, k] = (v[i, j, kp] - v[i, j, k]) / dx
    return vdx, vdy, vdz

def cgrad2d(v, dx):
    vdx = np.zeros_like(v)
    vdy = np.zeros_like(v)
    nx, ny = v.shape
    for i in range(nx):
        for j in range(ny):
            ip = (i + 1) % nx
            jp = (j + 1) % ny
            vdx[i, j] = (v[ip, j] - v[i, j]) / dx
            vdy[i, j] = (v[i, jp] - v[i, j]) / dx
    return vdx, vdy

def cdiv(vx, vy, vz, dx):
    vdiv = np.zeros_like(vx)
    nx, ny, nz = vx.shape
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                ip = (i + 1) % nx
                jp = (j + 1) % ny
                kp = (k + 1) % nz
                vdiv[i, j, k] = (vx[ip, j, k] - vx[i, j, k] + vy[i, jp, k] - vy[i, j, k] + vz[i, j, kp] - vz[i, j, k]) / dx
    return vdiv

def pcurl(ax, ay, az, bx, by, bz):
    cx = ay * bz - az * by
    cy = az * bx - ax * bz
    cz = ax * by - ay * bx
    return cx, cy, cz

def ccurl(vx, vy, vz, dx):
    vdx = np.zeros_like(vx)
    vdy = np.zeros_like(vx)
    vdz = np.zeros_like(vx)
    nx, ny, nz = vx.shape
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                ip = (i + 1) % nx
                jp = (j + 1) % ny
                kp = (k + 1) % nz
                vdx[i, j, k] = (vy[i, j, kp] - vy[i, j, k] - vz[i, jp, k] + vz[i, j, k]) / dx
                vdy[i, j, k] = (vz[ip, j, k] - vz[i, j, k] - vx[i, j, kp] + vx[i, j, k]) / dx
                vdz[i, j, k] = (vx[i, jp, k] - vx[i, j, k] - vy[ip, j, k] + vy[i, j, k]) / dx
    return vdx, vdy, vdz

