#! /usr/bin/env python

import numpy as np
from numpy.polynomial.polynomial import polyfit
from numpy.polynomial import chebyshev
from numpy.polynomial import legendre  

from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline as USplineFit
from scipy.optimize import curve_fit
from scipy.signal import gaussian
from scipy.optimize import fmin

import matplotlib.pyplot as plt

'''
fit_me:  y-values to be fit
xi:      initial x-values of fit_me
xf:      final x-values at which fit is to be calculated

fitted:  outputs chosen function calculated at xf
'''

def poly_fit(fit_me, xi, xf, weights,order=10, *args, **kwargs):
  # Fit returns coefficients, highest-order first
  pfit = polyfit(xi, fit_me, deg=order,w=weights)

  # Generate poly function
  pfunc = np.poly1d(pfit)
  
  # Evaluate at 'xf'
  fitted = pfunc(xf)
  residual = pfunc(xi) - fit_me
  return (fitted, residual)

def exp_fit(fit_me, xi, xf, weights=None,*args, **kwargs):
  # Make an exponential function
  expfunc = lambda x,a,b,c: a*np.exp(-b*x) + c

  # Initial guess
  a0 = fit_me[0]
  b0 = 1e-3     # GAH, hardcoding
  c0 = 0
  
  # Fit curve to fit_me
  if weights:
    weights = 1/weights
  popt, pcov = curve_fit(expfunc, xi, fit_me, p0=[a0,b0,c0],sigma=weights)

  # Interpolate for 'xf'
  fitted = expfunc(xf, *popt)
  residual = expfunc(x, *popt) - fit_me
  return (fitted, residual)

def gaussian_fit(fit_me, xi, xf, *args, **kwargs):
  pass  # Shit be broke
  '''
  gaussfunc = lambda x,a,mu,s2: a * np.exp(-((x-mu)**2)/(s2))
  popt, pcov = curve_fit(gaussfunc,xi,fit_me,p0=p0)

  fitted = gaussfunc(xf,*popt)

  '''
  
def lorentz_fit(fit_me, xi, xf, weights,*args, **kwargs):
  lorentzfunc = lambda x,x0,g: g/np.pi/((x-x0)**2+g**2)
  
  p0 = [xi[0],fit_me[0]]
  popt, pcov = curve_fit(lorentzfunc,xi,fit_me,p0=p0,sigma=1/weights)
  fitted = lorentzfunc(xf,*popt)
  
  residual = lorentzfunc(xi, *popt) - fit_me

  return (fitted, residual)
  
def interp1d_fit(fit_me, xi, xf, *args):
  f = interp1d(xi, fit_me, kind='cubic')
  
  fitted = np.zeros(xf.shape)
  fitted[:-1] = f(xf[:-1])[:]
  fitted[-1] = (fitted[-1])

  # Should return residual as second thingie in the future  
  return (fitted, fitted)
  
def chebyshev_fit(fit_me, xi, xf, weights,order=7, *args, **kwargs):
  chebfunc = chebyshev.chebfit(xi,fit_me,deg=order,w=weights)
  fitted = chebyshev.chebval(xf,chebfunc)
  residual = chebyshev.chebval(xi, chebfunc) - fit_me
  return (fitted, residual)

def legendre_fit(fit_me, xi, xf, weights,order=7,*args, **kwargs):
  legfunc = legendre.legfit(xi,fit_me,deg=order,w=weights)
  fitted = legendre.legval(xf,legfunc)
  residual = legendre.legval(xi, legfunc) - fit_me
  return (fitted, residual)
  
def spline_fit(fit_me, xi, xf,weights,order=3,s=1,*args,**kwargs):
  if order > 5:
    order = 5
    print "WARNING: spline_fit: 'k' order cannot exceed '5'"
    
  splinefunc = USplineFit(xi,fit_me, k=order,s=s,w=weights)
  fitted = splinefunc(xf)
  residual = splinefunc(xi) - fit_me

  return (fitted, residual)

def blackbody_fit(fit_me, xi, xf, *args, **kwargs):
  bbfunc = lambda xi, a, b, c: (a/xi**5) * 1.0 / (np.exp(b/xi)-1.) + c
    
  p0 = [10., 1e-2, 0.]

  weights = [1 if x > 5000 else 10000 for x in xi ]
  popt, pconv = curve_fit(bbfunc, xi, fit_me, p0=p0, sigma=weights)

  fitted = bbfunc(xf, *popt)
  residual = bbfunc(xi, *popt) - fit_me

  #plt.plot(xi, fit_me, 'k')
  #plt.plot(xi, fitted, 'r')
  #plt.show()

  return (fitted, residual)

