# from astropy.models import Gaussian2D
from numpy import pi, sin, cos, exp, linspace, logspace, zeros, indices, empty
from numpy.random import uniform, normal

def gaussian2D(yy, xx, amp, yc, xc, ys, xs):
  ychi  = (yy-yc) / ys
  xchi  = (xx-xc) / xs
  
  return amp*exp(-0.5*(ychi**2 + xchi**2))

def sineWave(times, sAmp, cAmp, angfreq):
  return sAmp*sin(angfreq*times) + cAmp*cos(angfreq*times)

def stellarModels(times, nTimes, nStars, minLogPeriod, maxLogPeriod, nPeriods, stdAmp):
  vPeriods  = logspace(minLogPeriod, maxLogPeriod, nPeriods)
  vAngFreqs = 2*pi / vPeriods

  vSinAmps  = normal(0, stdAmp, (nStars,nPeriods))
  vCosAmps  = normal(0, stdAmp, (nStars,nPeriods))

  starModels    = zeros((nStars, nTimes))
  for star, sAmps, cAmps in zip(starModels, vSinAmps, vCosAmps):
    for sAmp, cAmp, angfreq in zip(sAmps, cAmps, vAngFreqs):
      star += sineWave(times, sAmp, cAmp, angfreq)
  
  return starModels

def generate_field_of_view(nTimes=100, imageSize=100, nStars=10, 
                           fwhm=3., tMin=0, tMax=10, nPeriods=5,
                           minLogPeriod=-3, maxLogPeriod=2,
                           stdAmp=5e-3, returnAll=False):
  
  # set up size of data cube
  imageCube = empty((nTimes, imageSize, imageSize))

  # set up time series variability to inject
  times = linspace(tMin, tMax, nTimes)
  
  # Compute stellar variability per star
  starModels  = stellarModels(times, nTimes, nStars, minLogPeriod, maxLogPeriod, nPeriods, stdAmp)
  
  # Set up FOV -- stellar positions and amplitudes 
  sAmplitudes = uniform(10,100, nStars)

  ycenters  = uniform(0,imageSize,nStars)
  xcenters  = uniform(0,imageSize,nStars)

  ywidths   = normal(fwhm, 1e-2*fwhm, nStars)
  xwidths   = normal(fwhm, 1e-2*fwhm, nStars)

  yy,xx = indices((imageSize,imageSize))
  image = zeros((imageSize,imageSize))

  for kt, _ in enumerate(times):
    for ks, (amp, yc, xc, ys, xs) in enumerate(zip(sAmplitudes, ycenters, xcenters, ywidths, xwidths)):
      imageCube[kt] += gaussian2D(yy, xx, amp, yc, xc, ys, xs)*starModels[ks][kt]

  if returnAll:
    return imageCube, starModels
  else:
    return imageCube

if __name__ == '__main__':
  from pylab import imshow, ion
  ion()
  
  imageCube, starModels = generate_field_of_view(returnAll=True)
  
  imshow(imageCube[0])
