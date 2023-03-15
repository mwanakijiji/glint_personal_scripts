#!/usr/bin/env python
# coding: utf-8

import numpy as np


# ### Find thickness of MLA substrate
# Made from notebook 2023 03 15 by E.S.

# user inputs
NA = 0.125 # numerical aperture 
D = 66. # lenslet diameter (um)
n_substrate = 1.5255 # index of refraction of substrate
sig_NA = 0.005 # variation of NA (based on Simon email, conserv. sig_NA=0.005, more leeway sig_NA=0.015)
sig_D = 1 # variation of D
sig_n_substrate = 0.001 # variation of n

# print info
#   note based on plot from Simon, NA could be
#   conservative: 0.125 pm 0.005
#   more leeway: 0.125 pm 0.015
print("------ Input values: ------")
print("NA =", NA)
print("D (um) =", D)
print("n_substrate =",n_substrate)
print("--- Tolerances: ---")
print("sig_NA =",sig_NA)
print("sig_D =",sig_D)
print("sig_n_substrate =",sig_n_substrate)
print("----------------------")
print("------ Results: ------")

def foc_length(D_pass,NA_pass,n_pass):
    # thickness of the substrate
    
    return D_pass/(2.*np.tan(np.arcsin(NA_pass/n_pass)))


f_test = foc_length(D_pass=D,NA_pass=NA,n_pass=n_substrate)

print("Substrate thickness (um):")
print("f=",f_test)

# ### Find tolerance in thickness, based on propagated error
# assumes no correlation between parameters (d'=partial deriv):
# sig_f**2 = ( (d'(f)/dD)*sig_D )**2 + ( (d'(f)/dNA)*sig_NA )**2 + ( (d'(f)/dn)*sig_n )**2

def dfdx(x_pass, y_pass, z_pass):
    
    return 1./( 2.* np.tan(np.arcsin(y_pass/z_pass)) )

def dfdy(x_pass, y_pass, z_pass):
    
    piece1 = (-x_pass/2.)*np.power(np.tan(np.arcsin(y_pass/z_pass)),-2.)
    piece2 = np.power(np.cos(np.arcsin(y_pass/z_pass)),-2.) # sec = 1./cos
    piece3 = ( 1./np.sqrt(1-(y_pass/z_pass)**2) ) * ( 1./z_pass )
    
    d2f_dxdy = piece1*piece2*piece3
    
    return d2f_dxdy

def dfdz(x_pass, y_pass, z_pass):
    
    piece1 = (-x_pass/2.)*np.power(np.tan(np.arcsin(y_pass/z_pass)),-2.)
    piece2 = np.power(np.cos(np.arcsin(y_pass/z_pass)),-2.) # sec = 1./cos
    piece3 = ( 1./np.sqrt(1-np.power(y_pass/z_pass,2)) ) * ( -y_pass/z_pass**2. )
    
    d2f_dxdy = piece1*piece2*piece3
    
    return d2f_dxdy

def prop_error(dfdx_pass,dfdy_pass,dfdz_pass,sigx_pass,sigy_pass,sigz_pass):
    # assumes no correlation in D and NA
    
    sigf_2 = np.power(dfdx_pass*sigx_pass,2.) + np.power(dfdy_pass*sigy_pass,2.) + np.power(dfdz_pass*sigz_pass,2.)
    
    return np.sqrt(sigf_2)

# example

# cosmetic renaming
x = D
y = NA
z = n_substrate
sigx = sig_D
sigy = sig_NA
sigz = sig_n_substrate

dfdx_example = dfdx(x,y,z)
dfdy_example = dfdy(x,y,z)
dfdz_example = dfdz(x,y,z)

tol_test = prop_error(dfdx_example,dfdy_example,dfdz_example,sigx,sigy,sigz)

print("Tolerance (um):")
print("sig_f=",tol_test)
