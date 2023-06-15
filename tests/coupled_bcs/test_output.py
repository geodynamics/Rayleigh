#!/usr/bin/env python3

from rayleigh_diagnostics import Shell_Slices, main_input
import numpy as np
import sys, os

error = False

def test_bcs(case):
  """
  For a given case check if the values contained in the shell slices match the expected values based on the coefficients in
  main_input
  """
  error = False

  mi = main_input(os.path.join(case, 'main_input'))
  bcs = mi.vals['boundary_conditions']
  slices = Shell_Slices('00000002', path=os.path.join(case, 'Shell_Slices/'))
  
  # values produced by the code
  chi1_rmin = slices.vals[:,:,1,slices.lut[10001],0]
  chi1_rmax = slices.vals[:,:,0,slices.lut[10001],0]
  T_rmin    = slices.vals[:,:,1,slices.lut[501],0]
  T_rmax    = slices.vals[:,:,0,slices.lut[501],0]
  chi2_rmin = slices.vals[:,:,1,slices.lut[10201],0]
  chi2_rmax = slices.vals[:,:,0,slices.lut[10201],0]

  # derivatives produced by the code
  dchi1dr_rmin = slices.vals[:,:,1,slices.lut[10004],0]
  dchi1dr_rmax = slices.vals[:,:,0,slices.lut[10004],0]
  dTdr_rmin    = slices.vals[:,:,1,slices.lut[507],0]
  dTdr_rmax    = slices.vals[:,:,0,slices.lut[507],0]
  dchi2dr_rmin = slices.vals[:,:,1,slices.lut[10204],0]
  dchi2dr_rmax = slices.vals[:,:,0,slices.lut[10204],0]

  # coefficients input
  chi1bc_rmax = float(bcs['chi_a_top(1)'].replace('d', 'e'))
  chi1bc_rmin = float(bcs['chi_a_bottom(1)'].replace('d', 'e'))

  dTdrbc_rmax          = float(bcs['dtdr_top'].replace('d', 'e')) 
  dtdr_chi1_coeff_rmax = float(bcs['dtdr_chi_coeff_top(1)'].replace('d', 'e'))

  Tbc_rmin             = float(bcs['t_bottom'].replace('d', 'e'))
  t_chi1_coeff_rmin    = float(bcs['t_chi_coeff_bottom(1)'].replace('d', 'e'))
  t_dchi1dr_coeff_rmin = float(bcs['t_dchidr_coeff_bottom(1)'].replace('d', 'e'))

  chi2bc_rmax             = float(bcs['chi_a_top(2)'].replace('d', 'e'))
  chi2_chi1_coeff_rmax    = float(bcs['chi_chi_coeff_top(2,1)'].replace('d', 'e'))
  chi2_dchi1dr_coeff_rmax = float(bcs['chi_dchidr_coeff_top(2,1)'].replace('d', 'e'))

  dchi2drbc_rmin          = float(bcs['dchidr_a_bottom(2)'].replace('d', 'e'))
  dchi2dr_chi1_coeff_rmin = float(bcs['dchidr_chi_coeff_bottom(2,1)'].replace('d', 'e'))

  # check if values produced match expected values based on input coefficients
  def failed(values, expected, tol=1.e-10):
    return np.any(np.abs(values - expected) > tol)

  error = error or failed(chi1_rmax, chi1bc_rmax)
  error = error or failed(chi1_rmin, chi1bc_rmin)

  error = error or failed(dTdr_rmax, dtdr_chi1_coeff_rmax*chi1_rmax + dTdrbc_rmax)
  error = error or failed(T_rmin, t_chi1_coeff_rmin*chi1_rmin + t_dchi1dr_coeff_rmin*dchi1dr_rmin + Tbc_rmin)

  error = error or failed(chi2_rmax, chi2_chi1_coeff_rmax*chi1_rmax + chi2_dchi1dr_coeff_rmax*dchi1dr_rmax + chi2bc_rmax)
  error = error or failed(dchi2dr_rmin, dchi2dr_chi1_coeff_rmin*chi1_rmin + dchi2drbc_rmin)

  return error

lerror = test_bcs('const')
if lerror:
  print("ERROR: const case produced unexpected output relative to input coefficients!")
error = error or lerror

if error:
  sys.exit(1)
sys.exit(0)

