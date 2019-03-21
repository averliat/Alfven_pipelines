
import numpy as np
from numpy.fft import ifftn
import numpy.random as npr

import tables as T

nx = 1 << 5
p = -11. / 3.

noise = 1e-3

# Wave vectors
k1_alias = np.arange(nx, dtype=np.float64)
k1 = np.where(k1_alias >= nx / 2, k1_alias - nx, k1_alias)
kx, ky, kz = np.meshgrid(k1, k1, k1)
k = np.sqrt(kx**2 + ky**2 + kz**2)

# Phase
phx = np.exp(2. * np.pi * 1j * npr.random((nx, nx, nx)))
phy = np.exp(2. * np.pi * 1j * npr.random((nx, nx, nx)))
phz = np.exp(2. * np.pi * 1j * npr.random((nx, nx, nx)))

# Noise
nox = noise ** npr.uniform(-1., 1., (nx, nx, nx))
noy = noise ** npr.uniform(-1., 1., (nx, nx, nx))
noz = noise ** npr.uniform(-1., 1., (nx, nx, nx))

# Fourier coefficients
vxf = phx * k**(p / 2.) * nox
vyf = phy * k**(p / 2.) * noy
vzf = phz * k**(p / 2.) * noz

# Remove NaNs at origin
vxf[0, 0, 0] = 0.
vyf[0, 0, 0] = 0.
vzf[0, 0, 0] = 0.

# Inverse transform
vx = ifftn(vxf)
vy = ifftn(vyf)
vz = ifftn(vzf)

# Keep only real part
vx = np.real(vx)
vy = np.real(vy)
vz = np.real(vz)

# Save cube
h5fd = T.open_file("fake_cube.hdf5", mode='w')
h5fd.create_array("/", "nx", nx)
h5fd.create_array("/", "vx", vx)
h5fd.create_array("/", "vy", vy)
h5fd.create_array("/", "vz", vz)
h5fd.close()

