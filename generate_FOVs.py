# from astropy.models import Gaussian2D
from numpy.random import uniform, normal

def gaussian2D(yy, xx, amp, yc, xc, ys, xs):
  ychi  = (yy-yc) / ys
  xchi  = (xx-xc) / xs
  
  return amp*exp(-0.5*(ychi**2 + xchi**2))

# set up size of data cube
nTimes    = 1000
imageSize = 1024
nStars    = 100
fwhm      = 3.0

# set up time series variability to inject
tMin      = 0
tMax      = 10
times     = np.linspace(tMin, tMax, nTimes)

minPeriod = -3
maxPeriod = 2.
nPeriods  = 10
vPeriods  = np.logspace(minPeriod, maxPeriod, nPeriods)
vAngFreqs = 2*pi / vPeriods

stdAmp    = 5e-3
vSinAmps  = normal(0, stdAmp, (nStars,nPeriods))
vCosAmps  = normal(0, stdAmp, (nStars,nPeriods))

starModels    = np.zeros((nStars, nTimes))
for star, sAmps, cAmps in zip(starModels, vSinAmps, vCosAmps):
  for sAmp, cAmp, angfreq in zip(sAmps, cAmps, vAngFreqs):
    star += sAmp*sin(angfreq*times/per) + cAmp*cos(angfreq*times)

# Set up FOV -- stellar positions and amplitudes 
sAmplitudes = uniform(10,100, nStars)

ycenters  = uniform(0,imageSize,nStars)
xcenters  = uniform(0,imageSize,nStars)

ywidths   = normal(fwhm, 1e-2*fwhm, nStars)
xwidths   = normal(fwhm, 1e-2*fwhm, nStars)

yy,xx = np.indices((imageSize,imageSize))
image = np.zeros((imageSize,imageSize))

for yc, xc, ys, xs in zip(amplitudes, ycenters, xcenters, ywidths, xwidths):
  image += gaussian2D(yy, xx, amp, yc, xc, ys, xs)
