from photutils import CircularAperture, Gaussian2D, Const2D
from photutils import DAOStarFinder, find_peaks

from astropy.modeling import models, fitting
from astropy.stats import sigma_clipped_stats

def extract_lightcurves_from_FOV(data, nFWHM=3.0, nPixelsSeeing=3.0, nSig = 5.0, sigma=3.0, iters=5)
  mean, median, std = sigma_clipped_stats(data, sigma=sigma, iters=iters)
  
  aperRad = 1.5*nPixelsSeeing
  psSize  = nFWHM*nPixelSeeing

  useDAO = False:
  if useDAO:
    daofind = DAOStarFinder(fwhm=nPixelsSeeing, threshold=nSig*std)
    sources = daofind(data - median)
  else:
    sources = find_peaks(data, threshold=nSig*std, box_size=2*psSize)

  postageStamps = []
  for source in sources:
    yc0,xc0 = source['ycentroid'], source['xcentroid']
    postageStamps.append(data[yc0-psSize:yc0+psSize, xc0-psSize, xc0+psSize])

  amplitude = 1.0
  offset    = median
  y_mean    = psSize
  x_mean    = psSize
  y_stddev  = nPixelsSeeing
  x_stddev  = nPixelsSeeing
  theta     = 0.0

  yy, xx    = np.indices((psSize, psSize)) - 0.5*psSize # center at zero for later
  initModel = Gaussian2D(amplitude=amplitude, x_mean=x_mean, y_mean=y_mean, x_stddev=x_sstddev, y_stddev=y_stddev, theta=theta)
  initModel = initModel + Const2D(amplitude=offset)

  fit_lvmq  = fitting.LevMarLSQFitter()

  positions   = []
  apertures   = []
  phot_tables = []
  gaussFits   = []
  for image in images:
    gaussFits.append([])
    for postageStamp in postageStamps:
      nowModel = initModel
      nowModel.amplitude = postageStamp.max()
      nowModel.y_mean = yc
      nowModel.x_mean = xc
      
      gaussFits[-1].append(fit_lvmq(nowModel, xx+xc,yy+yc, postageStamp))
    
    positions.append([gfit['xcentroid'], gfit['ycentroid'] for gfits in gaussFits[-1]])
    
    apertures.append(CircularAperture(positions[-1], r=aperRad))
    phot_tables.append(aperture_photometry(data, apertures))
  
  return phot_tables
