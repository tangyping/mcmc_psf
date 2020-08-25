import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling.models import Gaussian2D
from astropy.io import fits
import scipy.signal as ss
import os

os.system("rm *_med.fits")

p1_dat = np.loadtxt('p1_med.dat')
p2_dat = np.loadtxt('p2_med.dat')
p3_dat = np.loadtxt('p3_med.dat')

hdu = fits.open('signal_1.fits') #just copy head from signal_1.fits
header = hdu[0].header

l_bk=header['NAXIS1']
s_bk=header['NAXIS2']
nf = 5

p1_img = np.zeros([s_bk, l_bk])
p2_img = np.zeros([s_bk, l_bk])
p3_img = np.zeros([s_bk, l_bk])

for ix in range(s_bk):
  for iy in range(l_bk):
    p1_img[ix,iy] = p1_dat[ix*l_bk+iy]
    p2_img[ix,iy] = np.exp(p2_dat[ix*l_bk+iy])
    p3_img[ix,iy] = p3_dat[ix*l_bk+iy]

fits.writeto('p1_med.fits', p1_img, header)
fits.writeto('p2_med.fits', p2_img, header)
fits.writeto('p3_med.fits', p3_img, header)
