# -*- coding: utf-8 -*-

"""Collection of functions to produce RTTOV-type SSP."""

try:
  import numpy as np
except ImportError:
  raise Exception( 'Module *numpy* required but not found.' )

try:
    import typhon.arts.scattering as tas
    import typhon.arts.xml as tax
    use_typhon = True
except ImportError:
    #from . import typhontools as tt
    #use_typhon = False
    raise Exception( 'Non-typhon *SingleScatteringData* not yet implemented' )

try:
  from . import utils
except Exception as e:
  print('Module *rttov* requires module *utils*, but failed to import it:\n%s' %str(e))

try:
  from . import assp
except Exception as e:
  print('Module *rttov* requires module *assp*, but failed to import it:\n%s' %str(e))

# for debugging
#import pdb


def get_assp(# assp_import_ssdb parameters
             habit_id, orientation,
             size_range=None, size_type='dmax',
             freq_range=None, temp_range=None,
             allow_nodata=False,
             # assp_interp_t parameters
             t_grid=None, t_interpm='pchip', t_allow_extrap=False,min_tpoints=3,
             # assp_interp_f parameters
             f_grid=None, f_interpm='linear', f_allow_extrap=False, f_extend=True,
             # assp_interp_size parameters
             size_grid=None, s_interpm='pchip', s_allow_extrap=False):
  #2017-05-14 Jana Mendrok
  """Get data (in ARTS SSD and SMD format) over f, T, D for one habit.
  
  Imports single scattering data within given frequency, temperature, and size
  ranges from ARTS single scattering database, converts them to ARTS format, and
  interpolates to given grids.
  
  Parameters
  ----------
  A: variable type
    variable description
  
  Returns
  -------
  S: (flat) list of objects
    ARTS-type single scattering data on common f and T grids for multiple
    (sorted) sizes of one habit.
  M: (flat) list of objects
    ARTS-type scattering meta data corresponding to S.
  md: list
    mass-dimension relationship parameters alpha and beta.
  """
  #FIXME: should we keep track somehow of which data is original (directly from
  # the database) and which is interpolated/extrpolated/extend-filled?
  
  # Derive all SSP data (within D, f, T ranges) for one habit. Stored in 
  # ARTS-SSP type structure single scattering and meta data structures.
  S,M = assp.assp_import_ssdb(habit_id, orientation,
                              size_range=size_range, size_type=size_type,
                              freq_range=freq_range, temp_range=temp_range,
                              allow_nodata=allow_nodata)
  
  # Interpolate the SSP data to a common T-grid
  assp.assp_interp_t(S, new_t_grid=t_grid,
                     interpm=t_interpm, allow_extrap=t_allow_extrap,
                     min_tpoints=min_tpoints)

  # Interpolate the SSP data to a common f-grid.
  # If new f_grid is not given, then determine it here (this is part of
  # assp_interp_f, but we might need it before we execute assp_interp_f to 
  # determine up to which freq we have to fill data (and we don't want this
  # fill f_limit to be included in a new, auto-determined common f_grid).
  if (f_grid is None):
    all_f = np.array([])
    all_f,count = assp.flatappend_from_recursive_S(S,all_f,'f_grid',count=0)
    new_f_grid = np.unique( all_f )
  # Since we only calculate SSP for x up to 10, we might first want to make
  # sure that the data covers the full requested f_grid (particularly when we
  # don't allow extrapolation).
  if f_extend:
    assp.assp_extend_f(S,f_limit=new_f_grid[-1],zerofill=False)
  # Now do the interpolation. Onto the previously derived common f_grid.
  assp.assp_interp_f(S, new_f_grid=new_f_grid,
                     interpm=f_interpm,allow_extrap=f_allow_extrap)

  # Sort in size (and interpolate the SSP and meta data to a new size grid if requested)
  S,M = assp.assp_interp_size(S, M, new_size_grid=size_grid,
                              size_type=size_type,
                              interpm=s_interpm, allow_extrap=s_allow_extrap)
  
  habsum = utils.ssdb_summary( habit_id, printinfo=False )
  mD = [ habsum['ALPHA'], habsum['BETA'] ] # we keep them as string here
  
  return S,M,mD


def calc_rssp(S,M,mD):
  #2017-05-15 Jana Mendrok
  """Derives single scattering properties as required for RTTOV.
  
  description
  
  Parameters
  ----------
  S: (flat) list of objects
    ARTS-type single scattering data on common f and T grids for multiple
    (sorted) sizes of one habit.
  M: (flat) list of objects
    ARTS-type scattering meta data corresponding to S.
  md: list
    mass-dimension relationship parameters alpha and beta.
  
  Returns
  -------
  rssp: dict
    Dictionary holding RTTOV-type single scattering data (Cext, Csca, g, Cbsc)
    along with corresponding grids (f, T, D) and habit characterising
    parameters (m-D relation's a and b) for one habit.
  """
  nf = S[0].f_grid.size
  nT = S[0].T_grid.size
  nD = len(S)
  
  data = -np.ones((nf,nT,nD,4)) # 4 data parameters: Cext, Csca, g, Cbsc
  freq = S[0].f_grid
  temp = S[0].T_grid
  dmax = -np.ones(nD)
  dveq = -np.ones(nD)
  mass = -np.ones(nD)

  for i in np.arange(nD):
    assert( assp.db_ptype[S[i].ptype]==20 ), \
      'RTTOV can only handle totally randomly oriented particles, but element #%i is %s.' %(i,S[i].ptype)
    assert( (S[i].f_grid==freq).all() ), \
      'Frequency grids of all scattering elements must be the same, but element #%i deviates.' %i
    assert( (S[i].T_grid==temp).all() ), \
      'Temperature grids of all scattering elements must be the same, but element #%i deviates.' %i
    
    # FIXME: assert non-negative values for dmax, dveq, mass; Cext, Csca, Cbsc

    dmax[i] = M[i].diameter_max
    dveq[i] = M[i].diameter_volume_equ
    mass[i] = M[i].mass
    
    # Cext
    data[:,:,i,0] = np.squeeze(S[i].ext_mat_data)
    # Csca
    data[:,:,i,1] = np.squeeze(S[i].ext_mat_data-S[i].abs_vec_data)
    # g
    data[:,:,i,2] = assp2g(S[i])
    # Cbsc
    assert( S[i].za_grid[-1]==180. ), \
      'Last entry in zenith angle needs to be 180deg, but is %.3fdeg for element #%i.' %(S[i].za_grid[-1],i)
    data[:,:,i,3] = S[i].pha_mat_data[:,:,-1,0,0,0,0]
    
  # Cbsc needs rescaling by 4Pi (see atmlab.scattering.assp2backcoef)
  data[...,3] = 4*np.pi*data[...,3]

  rssp = {}
  rssp['data'] = data
  rssp['freq'] = freq
  rssp['temp'] = temp
  rssp['dmax'] = dmax
  rssp['dveq'] = dveq
  rssp['mass'] = mass
  rssp['mDparams'] = np.array(mD)
  return rssp


def write_rssp(rssp,filename,dirname=None):
  #2017-05-15 Jana Mendrok
  """Writes RTTOV-type single scattering data to Mie-table input file.
  
  description
  
  Parameters
  ----------
  rssp: dict
    Dictionary holding RTTOV-type single scattering data, corresponding grids,
    and habit characterising parameters for one habit.
  filename: str
    Name of file where to write data for this habit.
  dirname: str
    Name of folder where to write data.
  """
  if dirname is not None:
    try:
      import os.path
    except ImportError:
      raise Exception( 'Module *os.path* required but not found.' )
    assert( os.path.isdir(dirname) ), \
      '%s is not a directory.' %dirname
    ffname = os.path.join(dirname,filename)
  else:
    ffname = filename
  
  f = open(ffname,'wb')
  # header (dimensions, grids, jhabit params)
  f.write( b'#  nf,  nT,  nD\n' )
  f.write( b'%5i%5i%5i\n' %(rssp['freq'].size,rssp['temp'].size,rssp['dmax'].size) )
  f.write( b'# f-grid [Hz]\n' )
  np.savetxt( f, rssp['freq'].reshape(1,-1), fmt='%.6e' )
  f.write( b'# T-grid [K]\n' )
  np.savetxt( f, rssp['temp'].reshape(1,-1), fmt='%7.3f' )
  f.write( b'# Dmax-grid [m]\n' )
  np.savetxt( f, rssp['dmax'].reshape(1,-1), fmt='%.6e' )
  f.write( b'# Dveq-grid [m]\n' )
  np.savetxt( f, rssp['dveq'].reshape(1,-1), fmt='%.6e' )
  f.write( b'# mass-grid [kg]\n' )
  np.savetxt( f, rssp['mass'].reshape(1,-1), fmt='%.6e' )
  f.write( b'# a,b of m = a * Dmax^b\n' )
  np.savetxt( f, rssp['mDparams'].reshape(1,-1), fmt='%s' )
  # SSP data
  f.write( b'# Cext [m2], Csca [m2], g [-], Cbsc [m2] over (f,T,D) (f=outer, D=inner loop)\n' )
  np.savetxt( f, rssp['data'].reshape(-1,rssp['data'].shape[-1]), fmt='%.8e' )
  f.close()


# NOTE: This is an equivalent to atmlab's scattering.assp2g.
# FIXME: It should rather/as well become part of typhon. And maybe reside
#        somewhere else here (in a tools collection for non-typhon users?)
def assp2g(s):
  #2017-05-15 Jana Mendrok
  """Derives asymetry parameter g from ARTS-type single scattering data.
  
  For a normalised phase function (p), g equals the 4pi integral of p*cos(th),
  where th is the scattering angle. For pure isotropic scattering g = 0, while
  pure forward scattering has g=1.

  WARNING: this function does not handle the extreme cases of delta-function
  type of forward or backward scattering lobes. A g of zero is returned for
  these cases.
 
  Parameters
  ----------
  s: object
    ARTS-type single scattering data.
  
  Returns
  -------
  g: 2-D numpy array
    Asymmetry parameters. One value for per frequency and temperature in s.
  """
  # Check input  
  # FIXME: adapt for non-typhon use
  if use_typhon:
    assert( isinstance(s,tas.SingleScatteringData)==1 ), \
      'Only single element *s* is handled (i.e. length(s) must be 1).'
  else:
    raise Exception( 'Non-typhon *SingleScatteringData* handling not yet implemented.' )
  assert( s.ptype=='totally_random' ), \
    'So far just totally random orientation is handled.'
  assert( s.za_grid[0]==0. ), \
    'First value of s.za_grid must be 0.'
  assert( s.za_grid[-1]==180. ), \
    'Last value of S.za_grid must be 180.'
  
  nf = s.f_grid.size
  nT = s.T_grid.size
  g = np.zeros( (nf, nT) )

  # ARTS uses pure phase matrix values, and not a normalised phase function,
  # hence we need to include a normalisation.

  za_rad = np.pi/180. * s.za_grid
  azi_w = np.abs(np.sin(za_rad))   # Azimuthal weighting
  cterm = np.cos(za_rad)          # Avoid recalculate cos term

  for j in np.arange(nf): 
    for k in np.arange(nT):  
      p = np.squeeze( s.pha_mat_data[j,k,:,0,0,0,0] )
    
      # All pi factors disappear as they are part of both integrals
      normfac = np.trapz( p*azi_w, x=s.za_grid )
      # If zero, this means that p==0 and should indicate very small particles
      # that have g=0, ie for which we leave g at the inititalised value of 0.
      if (normfac!=0.):
        g[j,k] = np.trapz( cterm*p*azi_w, s.za_grid ) / normfac
  
  return g



#def myfunc(A):
#  #YYYY-MM-DD Author Name
#  """one line summary
#  
#  description
#  
#  Parameters
#  ----------
#  A: variable type
#    variable description
#  
#  Returns
#  -------
#  B: variable type
#    variable description
#  """
#  codeblock
#  return B
