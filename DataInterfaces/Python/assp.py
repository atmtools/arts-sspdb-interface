# -*- coding: utf-8 -*-

"""Collection of functions to produce ARTS-type SSP."""


try:
  import numpy as np
except ImportError:
  raise Exception( 'Module *numpy* required but not found.' )

try:
  import scipy.interpolate as spi
except ImportError:
  raise Exception( 'Module *scipy.interpolate* required but not found.' )

try:
  import typhon.arts.scattering as tas
  #import typhon.arts.xml as tax
  use_typhon = True
except ImportError:
  #from . import typhontools as tt
  use_typhon = False
  #raise Exception( 'Non-typhon *SingleScatteringData* not yet implemented' )

try:
  from . import utils
except ImportError:
  try: #relative import does not work from scripts in the package folder
    import utils
#  from Python import utils
#except ImportError:
#  try:
#    import utils
  except ImportError:
    raise Exception( 'Module *utils* of database interface required but not found.' )
except Exception as e:
  print('Module *assp* requires module *utils*, but failed to import it:\n%s' %str(e))

# for debugging
#import pdb



# ARTS ptype definitions
a_ptype = {
  20: "totally_random",
  30: "azimuthally_random",
  }
# inverted ARTS ptype dict
ai_ptype = dict(zip(a_ptype.values(),a_ptype.keys()))

# SSP-DB ptype definitions (in case ARTS and SSP deviate in future)
db_ptype = {
  "totally_random": 20,
  "azimuthally_random": 30,
  "random": 20, # as in order instance of SSP-DB
  }



def assp_import_ssdb(habit_id, orientation,
                     habit_folder=None,
                     size_range=None, size_type='dveq',
                     freq_range=None, temp_range=None,
                     allow_nodata=False):
  #2017-03-20 Jana Mendrok
  #2017-10-16 Jana Mendrok: Added functionality for azimuthally random orientation
  #                         following M.Brath's draft.
  """Reads data from SSDB database for one habit and compiles it into ARTS format data.
  
  Default is to read all data for a specified habit and orientation
  combination, but the reading can be restricted in terms of size, frequency
  and temperature. For frequency and temperature, data are selected as: 
    limit1 <= data <= limit2
  while for frequency data at the lower limit is excluded:
    limit1 < data <= limit2.

  Default is to issue an error as soon as data are missing for a frequency
  and temperature combination. Allow missing data by setting the optional
  argument *allow_nodata* to true. Note that sets of frequencies and
  temperatures are allowed to differ between sizes, independently of how
  *allow_nodata* is set. For example, for one size there can be a single
  temperature, while other sizes have a temperature grid with several
  elements.
  
  Parameters
  ----------
  habit_id: int
    Habit id number.
  orientation: str
    Descriptor of orientation to explore for the chosen habit.
  habit_folder: str
    Full path to a habit folder.
    Temporary add for testing az.random not-yet-in-DB-structure. If used,
    provide a dummy habit_id.
  allow_nodata: bool
    See above. Default is false.
  size_range: 2-element list
    Particle size limits [unit: m]. Ignore data outside these limits. If not
    given, no limits are applied, i.e. all available data is considered.
  size_type: str
    Quantity used for size cropping. Allowed options are 'dveq' (default),
    'dmax' and 'mass'.
  freq_range: 2-element list
    Frequency limits [unit: Hz]. Ignore data outside these limits. If not
    given, no limits are applied, i.e. all available data is considered.
  temp range: 2-element list
    Temperature limits [unit: K]. Ignore data outside these limits. If not
    given, no limits are applied, i.e. all available data is considered.

  Returns
  -------
  S: list of objects (of SingleScatteringData)
    Array of ARTS-type single scattering data, one entry per particle size.
  M: list of objects (of ScatteringMetaData)
    Array of ARTS-type scattering meta data, one entry per particle size.
  """
  # import data
  if (habit_folder is None):
    data = utils.ssdb_import_habit( habit_id, orientation,
                                  allow_nodata=allow_nodata,
                                  size_range=size_range, size_type=size_type,
                                  freq_range=freq_range, temp_range=temp_range )
  else:
    data = utils.ssdb_import_habit( habit_id, orientation,
                                    habit_folder=habit_folder,
                                  allow_nodata=allow_nodata,
                                  size_range=size_range, size_type=size_type,
                                  freq_range=freq_range, temp_range=temp_range )

  # Loop particle sizes and convert data
  S = []
  M = []
  print( 'Converting extracted SSDB data to ARTS format data.' )
  print( '%i scattering elements to process' %len(data) )
  for i in np.arange(len(data)):
    if ((i+1)%5==0):
      print( '  processing element %i' %(i+1) )
    try:
      if (data[i]['SSD'].size>0):
        S1,M1 = ssdb2assp( data[i]['SSD'], data[i]['freq'],
                           data[i]['temp'], data[i]['nodata'] )
        S.append( S1 )
        M.append( M1 )
      else:
        try:
          print('For this particle (%s=%.0f um), no SSD available -'
                ' skipping this particle (ie. no equivalent entry in S & M).\n'
                %(size_type,data[i][size_type]*1e6))
        except:
          print('For particle #%i, no data available -'
                ' skipping this particle (ie. no equivalent entry in S & M).\n'
                %(i))
    except KeyError:
        print('For particle #%i, no data available -'
              ' skipping this particle (ie. no equivalent entry in S & M).\n'
              %(i))

  return S, M


def ssdb2assp(SSD, freq, temp, nodata,interpm='linear'):
  #2017-03-22 Jana Mendrok
  #2017-11-28 Robin Ekelund: Fixed sph importing error
  #2021-02-19 Vasileios Barlakas: Extended towards liquid oriented hydrometeors
  """Converts internal SSD to ARTS single scattering properties.

  The input SSD data should be imported through utils.ssdb_import_data using
  the grid_sort option.

  Parameters
  ----------
  SSD: 2-D numpy array of dictionary (dim: [freq,temp])
    Single scattering data over frequencies and temperatures.
  freq: 1-D numpy array
    Frequency grid of the SSD data.
  temp: 1-D numpy array
    Temperature grid of the SSD data.
  nodata: numpy array
    Boolean matrix, flagging empty SSD elements.
  interpm: str
    Method for zenith angle interpolation.

    Allowed are:
     TRO: 'linear', 'nearest', 'zero', 'slinear', 'quadratic','cubic'
     from scipy's interp1d (where the latter four use spline interpolation of
     zeroth to third order, respectively) as well as 'pchip' applying scipy's
     PchipInterpolator (supposed to provide similar, but not necessarily
     identical results to Matlab's pchip interpolation from interp1).
     ARO (liquid only): 'linear' and 'nearest' on the basis of RegularGridInterpolator
    Default is 'linear'.

  Returns
  -------
  S: object (typhon.arts.scattering.SingleScatteringData
    ARTS-type single scattering data of a specific scattering element (SingleScatteringData).
  M: object (typhon.arts.scattering.ScatteringMetaData)
    ARTS-type scattering meta data of a specific scattering element (ScatteringMetaData).
  """
  #####
  #FIXME: implement non-typhon version!
  #####

  # Basic sanity check of input
  if not( SSD.shape[0]==freq.size ):
    raise Exception( 'Mismatch in size between *SSD* and *freq*.' )
  if not( SSD.shape[1]==temp.size ):
    raise Exception( 'Mismatch in size between *SSD* and *temp*.' )
  if not( all('SingleScatteringData' in d for d in SSD.reshape(-1)) and
          all('ShapeData' in d for d in SSD.reshape(-1)) ):
    raise Exception( \
      'At least one element of *SSD* does not seem to hold the expected data format.' )

  # Set M
  #  (here from the first element in SSD; later we check that the crucial meta
  #   data is consistent between all SSD)
  if use_typhon:
    M = tas.ScatteringMetaData()
    if (SSD.size>0):
      M.description         = 'Meta data for '+SSD[0,0]['ShapeData']['description']
      M.source              = SSD[0,0]['ShapeData']['source']
      M.refr_index          = SSD[0,0]['ShapeData']['refrIndex_model']
      assert( (len(SSD[0,0]['ShapeData']['mass']['dim'])==0) and
              (SSD[0,0]['ShapeData']['mass']['data'].size==1) ), \
        'Mass data has wrong dimensions.'
      assert( (len(SSD[0,0]['ShapeData']['diameter_max']['dim'])==0) and
              (SSD[0,0]['ShapeData']['diameter_max']['data'].size==1) ), \
        'Max. diameter data has wrong dimensions.'
      assert( (len(SSD[0,0]['ShapeData']['diameter_vol_eq']['dim'])==0) and
              (SSD[0,0]['ShapeData']['diameter_vol_eq']['data'].size==1) ), \
        'Vol. equ. diameter data has wrong dimensions.'
      assert( (len(SSD[0,0]['ShapeData']['diameter_area_eq_aerodynamical']['dim'])==0) and
              (SSD[0,0]['ShapeData']['diameter_area_eq_aerodynamical']['data'].size==1) ), \
        'Aerodyn. area equ. diameter data has wrong dimensions.'
      M.mass                = np.float(SSD[0,0]['ShapeData']['mass']['data'])
      M.diameter_max        = np.float(SSD[0,0]['ShapeData']['diameter_max']['data'])
      M.diameter_volume_equ = np.float(SSD[0,0]['ShapeData']['diameter_vol_eq']['data'])
      M.diameter_area_equ_aerodynamical \
                            = np.float(SSD[0,0]['ShapeData']['diameter_area_eq_aerodynamical']['data'])
  else:
    raise Exception( 'Non-typhon *ScatteringMetaData* not yet implemented.' )

  phase = SSD[0,0]['ShapeData']['phase']

  # Basic data of S
  if use_typhon:
    S = tas.SingleScatteringData()
    S.version     = 3
    if (SSD.size>0):
      try:
        ptypeID = db_ptype[SSD[0,0]['SingleScatteringData']['orient_type']]
      except KeyError:
        raise Exception( \
          "Database ptype '%s' is unknown." %SSD[0,0]['SingleScatteringData']['orient_type'] )
      try:
        S.ptype       = a_ptype[ptypeID]
      except KeyError:
        raise Exception( \
          "No ARTS ptype equivalent defined for database ptype '%s' (internal ID=%i)." \
          %(SSD[0,0]['SingleScatteringData']['orient_type'],ptypeID) )
      S.description = SSD[0,0]['ShapeData']['description']
      S.f_grid  = freq
      S.T_grid  = temp
  else:
    raise Exception( 'Non-typhon *SingleScatteringData* not yet implemented' )

  # Fill S according to ptype
  nf = freq.size
  nt = temp.size

  if (SSD.size>0):

    #####
    # totally random orientation
    #####
    if ( S.ptype==a_ptype[20] ):

      pha_inds = [ [1,2,0,0], [2,3,0,0], [0,0,4,-5], [0,0,5,6] ]

      # Empty aa_grid
      S.aa_grid = np.array([])

      # Create a common za_grid
      # (instead of using union1d over and over, we concatenate all za_scat data
      #  and derive the unique value from complete array. to be tested whether
      #  this is indeed faster than union1d use.)
      za_grid = np.array([])
      for f in np.arange(freq.size):
        for t in np.arange(temp.size):
          if (not nodata[f,t]):
            # A number of asserts to double-check the input
            assert( abs(M.mass-SSD[f,t]['ShapeData']['mass']['data'])<1e-12 ), \
              'Mass info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
              %(freq[f]*1e-9,temp[t])
            assert( abs(M.diameter_max-SSD[f,t]['ShapeData']['diameter_max']['data'])<1e-12 ), \
              'Dmax info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
              %(freq[f]*1e-9,temp[t])
            assert( abs(M.diameter_volume_equ-SSD[f,t]['ShapeData']['diameter_vol_eq']['data'])<1e-12 ), \
              'Dveq info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
              %(freq[f]*1e-9,temp[t])

            try:
              ptypeID = db_ptype[SSD[f,t]['SingleScatteringData']['orient_type']]
            except KeyError:
              raise Exception( \
                "Database ptype for SSD at f=%.1fGHz and T=%.1fK ('%s') is unknown." \
                %(freq[f]*1e-9,temp[t],SSD[f,t]['SingleScatteringData']['orient_type']) )
            try:
              ptypeName = a_ptype[ptypeID]
            except KeyError:
              raise Exception( \
                "No ARTS ptype equivalent defined for database ptype '%s' (internal ID=%i) for SSD at f=%.1fGHz and T=%.1fK." \
                %(SSD[f,t]['SingleScatteringData']['orient_type'],ptypeID,freq[f]*1e-9,temp[t]) )
            assert( S.ptype==ptypeName ), \
              'ptype for SSD at f=%.1fGHz and T=%.1fK inconsistent with initial value.' \
              %(freq[f]*1e-9,temp[t])

            if not( abs(freq[f]-SSD[f,t]['SingleScatteringData']['frequency']['data'])<1e3 ):
              raise Exception( \
                'Freq info in SSD at f=%.1fGHz and T=%.1fK inconsistent with freq grid value.' \
                %(freq[f]*1e-9,temp[t]) )
            if not( abs(temp[t]-SSD[f,t]['SingleScatteringData']['temperature']['data'])<0.001 ):
              raise Exception( \
                'Temp info in SSD at f=%.1fGHz and T=%.1fK inconsistent with temp grid value.' \
                %(freq[f]*1e-9,temp[t]) )
            assert( pha_inds==(SSD[f,t]['SingleScatteringData']['phaMat_index']['data']).tolist() ), \
              'phamat indexing info in SSD at f=%.1fGHz and T=%.1fK inconsistent ptype.' \
              %(freq[f]*1e-9,temp[t])

            assert( (SSD[f,t]['SingleScatteringData']['za_scat']['data'][0]==0.) and
                    (SSD[f,t]['SingleScatteringData']['za_scat']['data'][-1]==180.) ), \
              'za_scat at f=%.1fGHz and T=%.1fK does not fulfill basic requirement\n' + \
              'of first value being at 0deg and/or last value at 180deg.' %(freq[f]*1e-9,temp[t])
            za_grid = np.append(za_grid,SSD[f,t]['SingleScatteringData']['za_scat']['data'])
      S.za_grid = np.unique(za_grid)
      nza = S.za_grid.size
      S.pha_mat_data = np.empty((nf,nt,nza,1,1,1,6))*np.NAN
      S.ext_mat_data = np.empty((nf,nt,1,1,1))*np.NAN
      S.abs_vec_data = np.empty((nf,nt,1,1,1))*np.NAN

      for f in np.arange(freq.size):
        for t in np.arange(temp.size):
          if (not nodata[f,t]):
            # Fill data fields

            # Old, linear za-interpolation only version:
            #for i in np.arange(S.pha_mat_data.shape[-1]):
            #  S.pha_mat_data[f,t,:,0,0,0,i] = \
            #    np.interp( S.za_grid,
            #              SSD[f,t]['SingleScatteringData']['za_scat']['data'],
            #              SSD[f,t]['SingleScatteringData']['phaMat_data']['data'][i,0,0,0,:] )

            # selectable za-interpolation method:
            if (interpm=='pchip'):
              if not( SSD[f,t]['SingleScatteringData']['za_scat'].size>2 ):
                raise Exception( \
                  "Interpolation method '%s' requires at least 3 sample points" + \
                  " for evaluating the function,\n" + \
                  "but SSD at f=%.1fGHz and T=%.1fK provides only %i" + \
                  " zenith angle sample points." \
                  %(interpm,freq[f]*1e-9,temp[t],i,
                    SSD[f,t]['SingleScatteringData']['za_scat'].size) )
              fi = spi.PchipInterpolator(SSD[f,t]['SingleScatteringData']['za_scat']['data'],
                                         SSD[f,t]['SingleScatteringData']['phaMat_data']['data'],
                                         axis=4)
              S.pha_mat_data[f,t,...] = fi(S.za_grid).T
            else:
              fi = spi.interp1d(SSD[f,t]['SingleScatteringData']['za_scat']['data'],
                                SSD[f,t]['SingleScatteringData']['phaMat_data']['data'],
                                axis=4,kind=interpm,assume_sorted=True)
              S.pha_mat_data[f,t,...] = fi(S.za_grid).T

            S.ext_mat_data[f,t,...] = SSD[f,t]['SingleScatteringData']['extMat_data']['data']
            S.abs_vec_data[f,t,...] = SSD[f,t]['SingleScatteringData']['absVec_data']['data']

    #####
    # azimuthally random orientation
    #####
    elif ( S.ptype==a_ptype[30] ):
      #####
      # Separate according to phase: liquid (gridded data) vs ice (spherical harmonics)
      #####
      if phase=='ice':

          try:
            from . import sph
          except ImportError:
            try: #relative import does not work from scripts in the package folder
              import sph
              # from Python import sph
            except ImportError:
              raise Exception( 'Module *sph* of database interface required but not found.' )
          except Exception as e:
             print('Module *assp* requires module *sph*, but failed to import it:\n%s' %str(e))

          S.za_grid=SSD[0,0]['SingleScatteringData']['za_inc']['data']

          for f in np.arange(freq.size):
            for t in np.arange(temp.size):
              if (not nodata[f,t]):
                # A number of asserts to double-check the input
                assert( abs(M.mass-SSD[f,t]['ShapeData']['mass']['data'])<1e-12 ), \
                  'Mass info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
                  %(freq[f]*1e-9,temp[t])
                assert( abs(M.diameter_max-SSD[f,t]['ShapeData']['diameter_max']['data'])<1e-12 ), \
                  'Dmax info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
                  %(freq[f]*1e-9,temp[t])
                assert( abs(M.diameter_volume_equ-SSD[f,t]['ShapeData']['diameter_vol_eq']['data'])<1e-12 ), \
                  'Dveq info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
                  %(freq[f]*1e-9,temp[t])

                try:
                  ptypeID = db_ptype[SSD[f,t]['SingleScatteringData']['orient_type']]
                except KeyError:
                  raise Exception( \
                    "Database ptype for SSD at f=%.1fGHz and T=%.1fK ('%s') is unknown." \
                    %(freq[f]*1e-9,temp[t],SSD[f,t]['SingleScatteringData']['orient_type']) )
                try:
                  ptypeName = a_ptype[ptypeID]
                except KeyError:
                  raise Exception( \
                    "No ARTS ptype equivalent defined for database ptype '%s' (internal ID=%i) for SSD at f=%.1fGHz and T=%.1fK." \
                    %(SSD[f,t]['SingleScatteringData']['orient_type'],ptypeID,freq[f]*1e-9,temp[t]) )
                assert( S.ptype==ptypeName ), \
                  'ptype for SSD at f=%.1fGHz and T=%.1fK inconsistent with initial value.' \
                  %(freq[f]*1e-9,temp[t])

                if not( abs(freq[f]-SSD[f,t]['SingleScatteringData']['frequency']['data'])<1e3 ):
                  raise Exception( \
                    'Freq info in SSD at f=%.1fGHz and T=%.1fK inconsistent with freq grid value.' \
                    %(freq[f]*1e-9,temp[t]) )
                if not( abs(temp[t]-SSD[f,t]['SingleScatteringData']['temperature']['data'])<0.001 ):
                  raise Exception( \
                    'Temp info in SSD at f=%.1fGHz and T=%.1fK inconsistent with temp grid value.' \
                    %(freq[f]*1e-9,temp[t]) )
                assert( (SSD[f,t]['SingleScatteringData']['za_inc']['data'][0]==0.) and
                        (SSD[f,t]['SingleScatteringData']['za_inc']['data'][-1]==180.) ), \
                  'za_inc at f=%.1fGHz and T=%.1fK does not fulfill basic requirement\n' + \
                  'of first value being at 0deg and/or last value at 180deg.' %(freq[f]*1e-9,temp[t])

                # for now we require all data to be on identical za_inc grids
                assert( (S.za_grid==SSD[f,t]['SingleScatteringData']['za_inc']['data']).all() ), \
                  'za_inc grid f=%.1fGHz and T=%.1fK inconsistent with global za_grid.' \
                  %(freq[f]*1e-9,temp[t])

          pha_inds = SSD[0,0]['SingleScatteringData']['phaMat_index']['data']
          ext_inds = SSD[0,0]['SingleScatteringData']['extMat_index']['data']
          abs_inds = SSD[0,0]['SingleScatteringData']['absVec_index']['data']

          # Set scattering grid sizes
          #   ARTS requires za_inc==za_sca, sph requires (as implemented)
          #    aa_sca[0-180]=za_sca. Hence here we fix everything to the za_inc, which
          #    is the coarser grid in our ADDA calcs.
          #   Basically, that could be replaced by two (independent) input parameters,
          #    one for za, one for aa. However, that requires more interpolations
          #    (and rotations?) in sph.
          nza = len(SSD[0,0]['SingleScatteringData']['za_inc']['data'])

          got_grid=False
          for f in np.arange(freq.size):
            for t in np.arange(temp.size):
              if (not nodata[f,t]):

                # extract optical properties
                ext_mat = SSD[f,t]['SingleScatteringData']['extMat_data']['data']
                abs_vec = SSD[f,t]['SingleScatteringData']['absVec_data']['data']

                phamat_real = SSD[f,t]['SingleScatteringData']['phaMat_data_real']['data']
                phamat_imag = SSD[f,t]['SingleScatteringData']['phaMat_data_imag']['data']
                phamat_sph = phamat_real+1j*phamat_imag

                # get truncation level
                lmax = sph.get_lmax(SSD[f,t]['SingleScatteringData']['sph_coeffs_scat']['data'])
                # transform spherical harmonics to grid
                [phase_matrix_reg,theta_s,phi_s] = \
                  sph.create_regular_representation(phamat_sph, lmax,
                                                    grid_size_theta=nza, gridtype='regular')
                assert( (S.za_grid==theta_s).all() ), \
                  'Spherical harmonics output polar grid at f=%.1fGHz and T=%.1fK' + \
                  ' inconsistent with global za_grid.' %(freq[f]*1e-9,temp[t])

                # allocate and set angle grids
                if got_grid==False:
                  S.aa_grid=phi_s
                  S.ext_mat_data = np.empty((nf,nt,len(S.za_grid),1,3))*np.NAN
                  S.abs_vec_data = np.empty((nf,nt,len(S.za_grid),1,2))*np.NAN
                  S.pha_mat_data = np.empty((nf,nt,len(S.za_grid),len(S.aa_grid),
                                            len(S.za_grid),1,16))*np.NAN
                  got_grid=True

                # Fill data fields
                S.ext_mat_data[f,t,...] = np.transpose(ext_mat,(2,1,0))
                S.abs_vec_data[f,t,...] = np.transpose(abs_vec,(2,1,0))
                S.pha_mat_data[f,t,...] = np.transpose(phase_matrix_reg,(3,4,2,1,0))

      elif phase=='liquid':
          S.aa_grid = np.array([])
          S.za_grid = np.array([])

          pha_inds = [ [1,5,9,13], [2,6,10,14], [3,7,11,15], [4,8,12,16] ]

          for f in np.arange(freq.size):
            for t in np.arange(temp.size):
              if (not nodata[f,t]):
                # A number of asserts to double-check the input
                assert( abs(M.mass-SSD[f,t]['ShapeData']['mass']['data'])<1e-12 ), \
                  'Mass info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
                  %(freq[f]*1e-9,temp[t])
                assert( abs(M.diameter_max-SSD[f,t]['ShapeData']['diameter_max']['data'])<1e-12 ), \
                  'Dmax info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
                  %(freq[f]*1e-9,temp[t])
                assert( abs(M.diameter_volume_equ-SSD[f,t]['ShapeData']['diameter_vol_eq']['data'])<1e-12 ), \
                  'Dveq info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
                  %(freq[f]*1e-9,temp[t])
                try:
                  ptypeID = db_ptype[SSD[f,t]['SingleScatteringData']['orient_type']]
                except KeyError:
                  raise Exception( \
                    "Database ptype for SSD at f=%.1fGHz and T=%.1fK ('%s') is unknown." \
                    %(freq[f]*1e-9,temp[t],SSD[f,t]['SingleScatteringData']['orient_type']) )
                try:
                  ptypeName = a_ptype[ptypeID]
                except KeyError:
                  raise Exception( \
                    "No ARTS ptype equivalent defined for database ptype '%s' (internal ID=%i) for SSD at f=%.1fGHz and T=%.1fK." \
                    %(SSD[f,t]['SingleScatteringData']['orient_type'],ptypeID,freq[f]*1e-9,temp[t]) )
                assert( S.ptype==ptypeName ), \
                  'ptype for SSD at f=%.1fGHz and T=%.1fK inconsistent with initial value.' \
                  %(freq[f]*1e-9,temp[t])

                if not( abs(freq[f]-SSD[f,t]['SingleScatteringData']['frequency']['data'])<1e3 ):
                  raise Exception( \
                    'Freq info in SSD at f=%.1fGHz and T=%.1fK inconsistent with freq grid value.' \
                    %(freq[f]*1e-9,temp[t]) )
                if not( abs(temp[t]-SSD[f,t]['SingleScatteringData']['temperature']['data'])<0.001 ):
                  raise Exception( \
                    'Temp info in SSD at f=%.1fGHz and T=%.1fK inconsistent with temp grid value.' \
                    %(freq[f]*1e-9,temp[t]) )
                assert( pha_inds==(SSD[f,t]['SingleScatteringData']['phaMat_index']['data']).tolist() ), \
                  'phamat indexing info in SSD at f=%.1fGHz and T=%.1fK inconsistent ptype.' \
                  %(freq[f]*1e-9,temp[t])

                assert( (SSD[f,t]['SingleScatteringData']['za_scat']['data'][0]==0.) and
                        (SSD[f,t]['SingleScatteringData']['za_scat']['data'][-1]==180.) ), \
                  'za_scat at f=%.1fGHz and T=%.1fK does not fulfill basic requirement\n' + \
                  'of first value being at 0deg and/or last value at 180deg.' %(freq[f]*1e-9,temp[t])
                S.za_grid = np.append(S.za_grid,SSD[f,t]['SingleScatteringData']['za_scat']['data'])
                S.aa_grid = np.append(S.aa_grid,SSD[f,t]['SingleScatteringData']['aa_scat']['data'])

          # Unique grid
          S.za_grid = np.unique(S.za_grid)
          S.aa_grid = np.unique(S.aa_grid)

          nza = S.za_grid.size
          naa = S.aa_grid.size
          npha= 16
          S.ext_mat_data = np.empty((nf,nt,nza,1,3))*np.NAN
          S.abs_vec_data = np.empty((nf,nt,nza,1,2))*np.NAN
          S.pha_mat_data = np.empty((nf,nt,nza,naa,nza,1,npha))*np.NAN

          for f in np.arange(freq.size):
            for t in np.arange(temp.size):
              if (not nodata[f,t]):

                exti = spi.interp1d(SSD[f,t]['SingleScatteringData']['za_scat']['data'],
                                    SSD[f,t]['SingleScatteringData']['extMat_data']['data'],
                                    axis=2,kind=interpm,assume_sorted=True)
                S.ext_mat_data[f,t,...] = exti(S.za_grid).T

                absi = spi.interp1d(SSD[f,t]['SingleScatteringData']['za_scat']['data'],
                                    SSD[f,t]['SingleScatteringData']['absVec_data']['data'],
                                    axis=2,kind=interpm,assume_sorted=True)
                S.abs_vec_data[f,t,...] = absi(S.za_grid).T

                zasq,aasq,zaiq = np.meshgrid(S.za_grid,S.aa_grid,S.za_grid, indexing='ij')

                for ind in np.arange(npha):
                    f_npha = SSD[f,t]['SingleScatteringData']['phaMat_data']['data'][ind,0,:,:,:]
                    fi = spi.RegularGridInterpolator((SSD[f,t]['SingleScatteringData']['za_scat']['data'],
                                                      SSD[f,t]['SingleScatteringData']['aa_scat']['data'],
                                                      SSD[f,t]['SingleScatteringData']['za_inc']['data']),f_npha, method = interpm)

                    # Change the order of the points to match the input shape of the interpolation function.
                    extended_grid = np.rollaxis(np.array([zasq,aasq,zaiq]),0,4)
                    interpolated  = fi(extended_grid)
                    S.pha_mat_data[f,t,:,:,:,0,ind] = interpolated

    #####
    # unknown particle type
    #####
    else:
      raise Exception( "Ptype '%s' not (yet?) implemented." %S.ptype )

  else: # if (SSD.size>0)
    print('No SSD data available.\n'
          'Returning default-filled SingleScatteringData and ScatteringMetaData'
          ' (grids partly empty, data fields all empty).')

  return S, M


def assp_interp_t(S,new_t_grid=None,interpm='pchip',allow_extrap=False,min_tpoints=3):
  #2017-03-28 Jana Mendrok
  """Temperature interpolation of ARTS single scattering properties.
  
  The function allows to select interpolation method and if extrapolation is
  allowed or not. The input data ARE allowed to have missing data, as long
  as missing data are set to NAN in the ext_mat_data field. It is allowed that
  the temperature corresponding to missing data differs between the
  frequencies. The number of missing data points is considered when checking
  if *min_tpoints* is met.
  
  Parameters
  ----------
  S: object or list of objects (of SingleScatteringData)
    ARTS-type single scattering data of one or several scattering elements.
  new_t_grid: 1-D array
    New temperature grid. Default is to create *new_t_grid* as the union of all
    S.T_grid. An interpolation is always performed, ie time is wasted if the
    data are already on the desired grid.
  interpm: str
    Interpolation method. Allowed are: 'linear', 'nearest', 'zero', 'slinear',
    'quadratic','cubic' from scipy's interp1d (where the latter four use spline
    interpolation of zeroth to third order, respectively) as well as 'pchip'
    applying scipy's PchipInterpolator (supposed to provide similar, but not
    necessarily identical results to Matlab's pchip interpolation from interp1).
    Default is 'pchip'.
  allow_extrap: bool
    Flag to allow extrapolation. Default is False.
  min_tpoints: int
    Minimum number of temperatures with data. Default is 3.

  Returns
  -------
  CAUTION: Implicitly returns modified(!) S, ie modifies the original S (make a
           deepcopy before if the original is further needed).
  """
  if (new_t_grid is None):
    all_t = np.array([])
    all_t,count = flatappend_from_recursive_S(S,all_t,'T_grid',count=0)
    new_t_grid = np.unique( all_t )
    #print('New f_grid created from %i individual grids:' %count)
  #print(new_f_grid)
  nt = len(new_t_grid)

  if not( isinstance(interpm,str) ):
    raise Exception( 'Provided *interpm* needs to be a string' )
  if not( interpm in
          ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic', 'pchip'] ):
    raise Exception( \
      "Interpolation method '%s' unknown. See doc for all allowed methods." %interpm )

  flatS = []
  flatS = shallowappend_from_recursive_S(S,flatS,'ALL')
  #print('Flattened S has %i entries.' %(len(flatS)))

  for i,s in enumerate(flatS):
    if not allow_extrap:
      if not( new_t_grid[0]>=s.T_grid[0] ):
        raise Exception( \
          'Lowest temperature in element #%i in flattened S is %.1f K, i.e. ' \
          'extrapolation is needed for t_new=%.1f K, while *allow_extrap* is set to False.' \
          %(i,s.T_grid[0],new_t_grid[0]) )
      if not( new_t_grid[-1]<=s.T_grid[-1] ):
        raise Exception( \
          'Highest temperature in element #%i in flattened S is %.1f K, i.e. ' \
          'extrapolation is needed for t_new=%.1f K, while *allow_extrap* is set to False.' \
          %(i,s.T_grid[-1],new_t_grid[-1]) )

    if not( s.T_grid.size >= min_tpoints ):
      raise Exception( \
        'You demand at least %i temperature points to allow an interpolation, ' \
        'but the length of the T_grid of element #%i in flattened S is just %i.' \
        %(min_tpoints, i, s.T_grid.size) )

    # Copy existing data 
    old_t_grid = s.T_grid
    e = s.ext_mat_data
    a = s.abs_vec_data
    p = s.pha_mat_data

    # Reallocate 
    s.T_grid = new_t_grid
    sh = list(e.shape)
    sh[1] = nt
    sh = tuple(sh)
    s.ext_mat_data = np.zeros(sh)
    sh = list(a.shape)
    sh[1] = nt
    sh = tuple(sh)
    s.abs_vec_data = np.zeros(sh)
    sh = list(p.shape)
    sh[1] = nt
    sh = tuple(sh)
    s.pha_mat_data = np.zeros(sh)

    for f in np.arange(s.f_grid.size):
      # Non-NaN temperature index    
      if (s.ptype==a_ptype[20]): #totally random => ext is 2D
        iok = np.where(~np.isnan(e[f,:]))[0]
      elif (s.ptype==a_ptype[30]): #azimuthally random => ext is 5D
        iok = np.where(~np.isnan(e[f,:,0,0,0]))[0]
      else:
        raise Exception( \
          "Unknow *ptype* ('%s') found in element %i of flattened S." %(s.ptype,i) )

      if not(len(iok) >= min_tpoints):
        raise Exception( \
          'You demand at least %i temperature points to allow an interpolation, ' \
          'but element #%i in flattened S contains only %i valid temperatures for frequency %.1 GHZ.' \
          %(min_tpoints, i, len(iok), s.f_grid[f]*1e-9) )
      
      if (interpm=='pchip'):
        if not( old_t_grid[iok].size>2 ):
          raise Exception( \
            "Interpolation method '%s' requires at least 3 sample points for evaluating the function,\n" \
            "but element #%i in flattened S provides only %i valid temperature sample points for frequency %.1 GHZ.." \
            %(interpm, i, old_t_grid[iok].size, s.f_grid[f]*1e-9) )
        fi = spi.PchipInterpolator(old_t_grid[iok],e[f,iok,...],axis=0,
                                                 extrapolate=True)
        s.ext_mat_data[f,:] = fi(new_t_grid)
        fi = spi.PchipInterpolator(old_t_grid[iok],a[f,iok,...],axis=0,
                                                 extrapolate=True)
        s.abs_vec_data[f,:] = fi(new_t_grid)
        fi = spi.PchipInterpolator(old_t_grid[iok],p[f,iok,...],axis=0,
                                                 extrapolate=True)
        s.pha_mat_data[f,:] = fi(new_t_grid)
      else:
        fi = spi.interp1d(old_t_grid[iok],e[f,iok,...],axis=0,kind=interpm,
                                        bounds_error=False,fill_value='extrapolate',assume_sorted=True)
        s.ext_mat_data[f,:] = fi(new_t_grid)
        fi = spi.interp1d(old_t_grid[iok],a[f,iok,...],axis=0,kind=interpm,
                                        bounds_error=False,fill_value='extrapolate',assume_sorted=True)
        s.abs_vec_data[f,:] = fi(new_t_grid)
        fi = spi .interp1d(old_t_grid[iok],p[f,iok,...],axis=0,kind=interpm,
                                       bounds_error=False,fill_value='extrapolate',assume_sorted=True)
        s.pha_mat_data[f,:] = fi(new_t_grid)


def assp_interp_f(S,new_f_grid=None,interpm='pchip',allow_extrap=False):
  #2017-03-27 Jana Mendrok
  """Frequency interpolation of ARTS single scattering properties.
  
  The function allows to select interpolation method and if extrapolation
  is allowed or not. The input data are NOT allowed to have missing data.
  
  Parameters
  ----------
  S: object or list of objects (of SingleScatteringData)
    ARTS-type single scattering data of one or several scattering elements.
  new_f_grid: 1-D array
    New frequency grid. If None, *new_f_grid* is created as the union of all
    S.f_grid. Interpolation is skipped for S entries where the f_grid is
    already identical to the new grid.
  interpm: str
    Interpolation method. Allowed are: 'linear', 'nearest', 'zero', 'slinear',
    'quadratic','cubic' from scipy's interp1d (where the latter four use spline
    interpolation of zeroth to third order, respectively) as well as 'pchip'
    applying scipy's PchipInterpolator (supposed to provide similar, but not
    necessarily identical results to Matlab's pchip interpolation from interp1).
    Default is 'pchip'.
  allow_extrap: bool
    Flag to allow extrapolation. Default is False.

  Returns
  -------
  CAUTION: Implicitly returns modified(!) S, ie modifies the original S (make a
           deepcopy before if the original is further needed).
  """
  if (new_f_grid is None):
    all_f = np.array([])
    all_f,count = flatappend_from_recursive_S(S,all_f,'f_grid',count=0)
    new_f_grid = np.unique( all_f )
    #print('New f_grid created from %i individual grids:' %count)
  #print(new_f_grid)
  nf = len(new_f_grid)

  if not( isinstance(interpm,str) ):
    raise Exception( 'Provided *interpm* needs to be a string' )
  if not( interpm in
          ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic', 'pchip'] ):
    raise Exception( \
      "Interpolation method '%s' unknown. See doc for all allowed methods." %interpm )

  flatS = []
  flatS = shallowappend_from_recursive_S(S,flatS,'ALL')
  #print('Flattened S has %i entries.' %(len(flatS)))

  for i,s in enumerate(flatS):
    if not allow_extrap:
      if not( new_f_grid[0]>=s.f_grid[0] ):
        raise Exception( \
          'Lowest frequency in element #%i in flattened S is %.1f GHz, i.e. ' \
          'extrapolation is needed for f_new=%.1f GHz, while *allow_extrap* is set to False.' \
          %(i,s.f_grid[0]*1e-9,new_f_grid[0]*1e-9) )
      if not( new_f_grid[-1]<=s.f_grid[-1] ):
        raise Exception( \
          'Highest frequency in element #%i in flattened S is %.1f GHz, i.e. ' \
          'extrapolation is needed for f_new=%.1f GHz, while *allow_extrap* is set to False.' \
          %(i,s.f_grid[-1]*1e-9,new_f_grid[-1]*1e-9) )

    if ( (len(s.f_grid)==nf) and (s.f_grid==new_f_grid).all() ):
      print('For element %i, original f_grid identical to new one. No interpolation performed.' %i)
    else:
      if (interpm=='pchip'):
        if not( s.f_grid.size>2 ):
          raise Exception( \
            "Interpolation method '%s' requires at least 3 sample points for evaluating the function,\n" \
            "but element #%i in flattened S provides only %i sample points." \
            %(interpm,i,s.f_grid.size) )
        fi = spi.PchipInterpolator(s.f_grid,s.abs_vec_data,axis=0,
                                                extrapolate=True)
        s.abs_vec_data = fi(new_f_grid)
        fi = spi.PchipInterpolator(s.f_grid,s.ext_mat_data,axis=0,
                                                extrapolate=True)
        s.ext_mat_data = fi(new_f_grid)
        fi = spi.PchipInterpolator(s.f_grid,s.pha_mat_data,axis=0,
                                                extrapolate=True)
        s.pha_mat_data = fi(new_f_grid)
      else:
        fi = spi.interp1d(s.f_grid,s.abs_vec_data,axis=0,kind=interpm,
                                       bounds_error=False,fill_value='extrapolate',assume_sorted=True)
        s.abs_vec_data = fi(new_f_grid)
        fi = spi.interp1d(s.f_grid,s.ext_mat_data,axis=0,kind=interpm,
                                       bounds_error=False,fill_value='extrapolate',assume_sorted=True)
        s.ext_mat_data = fi(new_f_grid)
        fi = spi.interp1d(s.f_grid,s.pha_mat_data,axis=0,kind=interpm,
                                       bounds_error=False,fill_value='extrapolate',assume_sorted=True)
        s.pha_mat_data = fi(new_f_grid)
      s.f_grid = new_f_grid


def assp_extend_f(S,f_limit,zerofill=False):
  #2017-04-12 Jana Mendrok
  """Extension of frequency coverage of ARTS single scattering properties.
  
  The function provides a fast way to ensure that the frequency coverage of the
  ASSP data reaches the given limit. If there is a gap between highest
  frequency and *f_limit*, the data are extended in one of two simple manners.

  Default is to copy the scattering properties at highest existing frequency,
  and include with *f_limit* as the assigned frequency.

  If *zerofill* is set to true, the data are instead extended with zeros.
  Here two frequencies are added. The first point is the highest existing
  frequency + 1 GHz, and the second point is *f_limit*. The value 1 GHz is
  adjusted when necessary to always get an increasing frequency grid.
  
  Parameters
  ----------
  S: object or list of objects (of SingleScatteringData)
    ARTS-type single scattering data of one or several scattering elements.
  f_limit: float
    The frequency limit to meet.
  zerofill: bool
    Flag to instead fill with zeros, see above. Default is false.

  Returns
  -------
  CAUTION: Implicitly returns modified(!) S, ie modifies the original S (make a
           deepcopy before if the original is further needed).
  """
  flatS = []
  flatS = shallowappend_from_recursive_S(S,flatS,'ALL')

  if zerofill:
    for i,s in enumerate(flatS):
      if f_limit>s.f_grid[-1]:
        off_f = s.f_grid[-1]+1e9
        mean_f = np.mean([s.f_grid[-1],f_limit])
        s.f_grid = np.append(s.f_grid, [np.min([off_f,mean_f]),f_limit])
        s.ext_mat_data = np.append(s.ext_mat_data,np.zeros_like(s.ext_mat_data)[:2,...],axis=0)
        s.abs_vec_data = np.append(s.abs_vec_data,np.zeros_like(s.abs_vec_data)[:2,...],axis=0)
        s.pha_mat_data = np.append(s.pha_mat_data,np.zeros_like(s.pha_mat_data)[:2,...],axis=0)
  else:
    for i,s in enumerate(flatS):
      if f_limit>s.f_grid[-1]:
        s.f_grid = np.append(s.f_grid, f_limit)
        s.ext_mat_data = np.append(s.ext_mat_data,s.ext_mat_data[-1:,...],axis=0)
        s.abs_vec_data = np.append(s.abs_vec_data,s.abs_vec_data[-1:,...],axis=0)
        s.pha_mat_data = np.append(s.pha_mat_data,s.pha_mat_data[-1:,...],axis=0)


def assp_interp_za(S,new_za_grid=None,interpm='linear'):
  #2017-03-28 Jana Mendrok
  """Zenith angle interpolation of ARTS single scattering properties.
  
  Input data must cover zenith angles of 0 and 180 deg (ie extrapolation should
  never be required). The new grid can either be set by the user, or
  interpolation to the common grid is performed.
  
  Parameters
  ----------
  S: object or list of objects (of SingleScatteringData)
    ARTS-type single scattering data of one or several scattering elements.
  new_za_grid: 1-D array
    New zenith angle grid. If None, *new_za_grid* is created as the union of
    all S.za_grid. Interpolation is skipped for S entries where the za_grid is
    already identical to the new grid.
  interpm: str
    Interpolation method. Allowed are: 'linear', 'nearest', 'zero', 'slinear',
    'quadratic','cubic' from scipy's interp1d (where the latter four use spline
    interpolation of zeroth to third order, respectively) as well as 'pchip'
    applying scipy's PchipInterpolator (supposed to provide similar, but not
    necessarily identical results to Matlab's pchip interpolation from interp1).
    Default is 'linear'.

  Returns
  -------
  CAUTION: Implicitly returns modified(!) S, ie modifies the original S (make a
           deepcopy before if the original is further needed).
  """
  if (new_za_grid is None):
    all_za = np.array([])
    all_za,count = flatappend_from_recursive_S(S,all_f,'za_grid',count=0)
    new_za_grid = np.unique( all_za )
    #print('New f_grid created from %i individual grids:' %count)
  #print(new_f_grid)
  if not( new_za_grid[0]==0. ):
    raise Exception( \
      'First point in *new_za_grid* must be 0.0 deg, but is %.1e deg.' %new_za_grid[0] )
  if not( new_za_grid[-1]==180. ):
    raise Exception( \
      'Last point in *new_za_grid* must be 180.0 deg, but is %.1e deg.' %new_za_grid[-1] )
  nza = len(new_za_grid)

  if not( isinstance(interpm,str) ):
    raise Exception( 'Provided *interpm* needs to be a string' )
  if not( interpm in
          ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic', 'pchip'] ):
    raise Exception( \
      "Interpolation method '%s' unknown. See doc for all allowed methods." %interpm )

  flatS = []
  flatS = shallowappend_from_recursive_S(S,flatS,'ALL')
  #print('Flattened S has %i entries.' %(len(flatS)))

  for i,s in enumerate(flatS):
    assert( s.za_grid[0]==0. ), \
      'First point in *za_grid* of element #%i in flattened S must be 0.0 deg, but is %.1e deg.' \
      %(i,s.za_grid[0])
    assert( s.za_grid[-1]==180. ), \
      'Last point in *za_grid* of element #%i in flattened S must be 180.0 deg, but is %.1e deg.' \
      %(i,s.za_grid[-1])
    if ( (len(s.za_grid)==nza) and (s.za_grid==new_za_grid).all() ):
      print('For element %i, original za_grid identical to new one. No interpolation performed.' %i)
    else:
      if (interpm=='pchip'):
        if not( s.za_grid.size>2 ):
          raise Exception( \
            "Interpolation method '%s' requires at least 3 sample points for evaluating the function,\n" \
            "but element #%i in flattened S provides only %i zenith angle sample points." \
            %(interpm,i,s.za_grid.size) )
        if (s.ptype==a_type[30]):
          fi = spi.PchipInterpolator(s.za_grid,s.abs_vec_data,axis=2)
          s.abs_vec_data = fi(new_za_grid)
          fi = spi.PchipInterpolator(s.za_grid,s.ext_mat_data,axis=2)
          s.ext_mat_data = fi(new_za_grid)
          fi = spi.PchipInterpolator(s.za_grid,s.pha_mat_data,axis=2)
          s.pha_mat_data = fi(new_za_grid)
          fi = spi.PchipInterpolator(s.za_grid,s.pha_mat_data,axis=4)
          s.pha_mat_data = fi(new_za_grid)
        elif (s.ptype==a_type[20]):
          fi = spi.PchipInterpolator(s.za_grid,s.pha_mat_data,axis=2)
          s.pha_mat_data = fi(new_za_grid)
        else:
          raise Exception( \
            "Unknow *ptype* ('%s') found in element %i of flattened S." %(s.ptype,i) )
      else:
        if (s.ptype==a_type[30]):
          fi = spi.interp1d(s.za_grid,s.abs_vec_data,axis=2,kind=interpm,
                                          assume_sorted=True)
          s.abs_vec_data = fi(new_za_grid)
          fi = spi.interp1d(s.za_grid,s.ext_mat_data,axis=2,kind=interpm,
                                          assume_sorted=True)
          s.ext_mat_data = fi(new_za_grid)
          fi = spi.interp1d(s.za_grid,s.pha_mat_data,axis=2,kind=interpm,
                                          assume_sorted=True)
          s.pha_mat_data = fi(new_za_grid)
          fi = spi.interp1d(s.za_grid,s.pha_mat_data,axis=4,kind=interpm,
                                          assume_sorted=True)
          s.pha_mat_data = fi(new_za_grid)
        elif (s.ptype==a_type[20]):
          fi = spi.interp1d(s.za_grid,s.pha_mat_data,axis=2,kind=interpm,
                                          assume_sorted=True)
          s.pha_mat_data = fi(new_za_grid)
        else:
          raise Exception( \
            "Unknow *ptype* ('%s') found in element %i of flattened S." %(s.ptype,i) )
      s.za_grid = new_za_grid


def assp_interp_size(S,M,new_size_grid=None,size_type='dveq',interpm='pchip',
                     allow_extrap=False):
  #2017-03-28 Jana Mendrok
  """Size interpolation of ARTS single scattering properties.
  
  The function allows to perform an interpolation in size, with some
  constraints on the input data: All elements in *S* must have common
  frequency, temperature and angular grids, and have the same "ptype".

  The interpolation can be done for three size variables
    dveq : M.diameter_volume_equ is used as size parameter
    dmax : M.diameter_max is used as size parameter
    mass : M.mass is used as size parameter

  Input data can be unsorted in size.

  Data matching zero size are added automatically in the interpolation, with
  all optical properties, as well as diameters and mass, being set to 0 for
  size zero. This zero size is not considered when checking sizes with
  *allow_extrap* set to false.
  
  Parameters
  ----------
  S: object or list of objects (of SingleScatteringData)
    ARTS-type single scattering data of one or several scattering elements.
  M: object or list of objects (of ScatteringMetaData)
    ARTS-type scattering meta data of one or several scattering elements.
  new_size_grid: 1-D array
    New size grid. Values and unit shall match selection of *size_type*.
    Default is to create *new_size_grid* as the union of all incoming sizes.
    An interpolation is always applied and time is wasted if the data are
    already on the common size grid.
  size_type: str
    Size variable defining the interpolation grid, see above. Allowed options
    are 'dveq' (default), 'dmax' and 'mass'.
  interpm: str
    Interpolation method. Allowed are: 'linear', 'nearest', 'zero', 'slinear',
    'quadratic','cubic' from scipy's interp1d (where the latter four use spline
    interpolation of zeroth to third order, respectively) as well as 'pchip'
    applying scipy's PchipInterpolator (supposed to provide similar, but not
    necessarily identical results to Matlab's pchip interpolation from interp1).
    Default is 'pchip'.
  allow_extrap: bool
    Flag to allow extrapolation. Default is False.

  Returns
  -------
  S: list of objects (typhon.arts.scattering.SingleScatteringData
    As input, but in a flat list and with data interpolated to new_size_grid.
  M: list of objects (typhon.arts.scattering.ScatteringMetaData)
    As input, but in a flat list and with data interpolated to new_size_grid.

  CAUTION: Also implicitly modifies(!) the original (input) S and M (make a
           deepcopy before if the originals are further needed).
  """
  if not( isinstance(size_type,str) ):
    raise Exception( '*size_type* must be a string.' )
  size_type = size_type.lower()
  if not( size_type in ['dveq', 'dmax', 'mass'] ):
    raise Exception( \
      "Valid options for *size_type* are 'dveq', 'dmax', and 'mass'." + \
      "You selected %s." %size_type )
 
  if not( isinstance(interpm,str) ):
    raise Exception( 'Provided *interpm* needs to be a string' )
  if not( interpm in
          ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic', 'pchip'] ):
    raise Exception( \
      "Interpolation method '%s' unknown. See doc for all allowed methods." %interpm )

  if (size_type=='mass'):
    attname='mass'
  elif (size_type=='dmax'):
    attname='diameter_max'
  else:
    attname='diameter_volume_equ'

  flatS = []
  flatS = shallowappend_from_recursive_S(S,flatS,'ALL')
  flatM = []
  flatM = shallowappend_from_recursive_S(M,flatM,'ALL',instancetype='SMD')
  if not( len(S)==len(M) ):
    raise Exception( '*S* and *M* must have the same length.' )

  old_size_grid = np.array([])
  for m in flatM:
    old_size_grid = np.append(old_size_grid, getattr(m,attname))
  if not( old_size_grid.size==np.unique(old_size_grid).size ):
    raise Exception( \
      "Size duplicates found in input data. Interpolation can't handle this." )

  # Sort original data in size
  ind = old_size_grid.argsort()
  old_size_grid = old_size_grid[ind]
  flatS = list(np.array(flatS)[ind])
  flatM = list(np.array(flatM)[ind])

  if (new_size_grid is None):
    #new_size_grid = old_size_grid
    # in this case we do not need to do anything. apart from returning the sorted object lists
    return flatS,flatM
  
  # if new_size grid is given:
  if not( new_size_grid.size==np.unique(new_size_grid).size ):
    raise Exception( \
      "Size duplicates found in *new_size_grid*. Interpolation can't handle this." )
  new_size_grid.sort()

  # Check interpolation limits
  if not allow_extrap:
    if not( new_size_grid[0]>=old_size_grid[0] ):
      raise Exception( \
        'Smallest size in S is %.2e and an extraplotation is needed, while ' \
        '*allow_extrap* is set to be false.' %old_size_grid[0] )
    if not( new_size_grid[-1]<=old_size_grid[-1] ):
      raise Exception( \
        'Largest size in S is %.2e and an extraplotation is needed, while ' \
        '*allow_extrap* is set to be false.' %old_size_grid[-1] )

  # Put together all data into merged arrays
  # (this is only possible if/since requiring all data to be on identical grids)
  
  # First, prepare the empty (zeroed) arrays
  old_with0 = np.append(np.zeros(1),old_size_grid)
  nos1 = old_grid_size.size+1
  sh = flatS[0].ext_mat_data[np.newaxis].shape
  sh[0] = nos1
  e = np.zeros(sh)
  sh = flatS[0].abs_vec_data[np.newaxis].shape
  sh[0] = nos1
  a = np.zeros(sh)
  sh = flatS[0].pha_mat_data[np.newaxis].shape
  sh[0] = nos1
  p = np.zeros(sh)
  si = np.zeros((nos1,3))
  ar = np.zeros(nos1)
  
  # Now loop over all elements in flatS/M, check for identical grids and copy
  # the data if check successful
  for i,s in enumerate(S):
    if (i>0):
      if not( s.ptype==S[0].ptype ):
        raise Exception( \
          'The particle type of all elements of *S* must be the same.' )
      if not( (s.f_grid==S[0].f_grid).all() ):
        raise Exception( \
          'The frequency grid of all elements of *S* must be the same.' )
      if not( (s.T_grid==S[0].T_grid).all() ):
        raise Exception( \
          'The temperature grid of all elements of *S* must be the same.' )
      if not( (s.za_grid==S[0].za_grid).all() ):
        raise Exception( \
          'The zenith angle grid of all elements of *S* must be the same.' )
      if not( (s.aa_grid==S[0].aa_grid).all() ):
        raise Exception( \
          'The azimuth angle grid of all elements of *S* must be the same.' )
    
    e[i+1,...] = s.ext_mat_data
    a[i+1,...] = s.abs_vec_data
    p[i+1,...] = s.pha_mat_data
    si[i+1,:] = np.array([M[i].mass,M[i].diameter_max,M[i].diameter_volume_equ])
    ar[i+1] = M[i].diameter_area_equ_aerodynamical

  # Interpolate
  if (interpm=='pchip'):
    if not( old_size_grid.size>2 ):
      raise Exception( \
        "Interpolation method '%s' requires at least 3 sample points for evaluating the function,\n" \
        "but only %i size sample points given." \
        %(interpm, old_size_grid.size) )
    fi = spi.PchipInterpolator(old_size_grid,e,axis=0,
                                             extrapolate=True)
    e = fi(new_size_grid)
    fi = spi.PchipInterpolator(old_size_grid,a,axis=0,
                                             extrapolate=True)
    a = fi(new_size_grid)
    fi = spi.PchipInterpolator(old_size_grid,p,axis=0,
                                             extrapolate=True)
    p = fi(new_size_grid)
    fi = spi.PchipInterpolator(old_size_grid,si,axis=0,
                                             extrapolate=True)
    si = fi(new_size_grid)
    # Special treatment here as M[i].diameter_area_equ_aerodynamical often is
    # set to be NaN
    if (np.isnan(ar).any()):
      ar = ar*np.nan
    else:
      fi = spi.PchipInterpolator(old_size_grid,ar,axis=0,
                                               extrapolate=True)
      ar = fi(new_size_grid)
  else:
    fi = spi.interp1d(old_size_grid,e,axis=0,kind=interpm,
                                    bounds_error=False,fill_value='extrapolate',assume_sorted=True)
    e = fi(new_size_grid)
    fi = spi.interp1d(old_size_grid,a,axis=0,kind=interpm,
                                    bounds_error=False,fill_value='extrapolate',assume_sorted=True)
    a = fi(new_size_grid)
    fi = spi.interp1d(old_size_grid,p,axis=0,kind=interpm,
                                    bounds_error=False,fill_value='extrapolate',assume_sorted=True)
    p = fi(new_size_grid)
    fi = spi.interp1d(old_size_grid,si,axis=0,kind=interpm,
                                    bounds_error=False,fill_value='extrapolate',assume_sorted=True)
    si = fi(new_size_grid)
    if (np.isnan(ar).any()):
      ar = ar*np.nan
    else:
      fi = spi.interp1d(old_size_grid,ar,axis=0,kind=interpm,
                                      bounds_error=False,fill_value='extrapolate',assume_sorted=True)
      ar = fi(new_size_grid)

  # Create new S and M
  defs = {'ptype': flatS[0].ptype,
          'version': flatS[0].version,
          'description': flatS[0].version,
          'T_grid': flatS[0].T_grid,
          'f_grid': flatS[0].f_grid,
          #'aspect_ratio': np.nan,
          'za_grid': flatS[0].za_grid,
          'aa_grid': flatS[0].aa_grid,
          }
  if use_typhon:
    s = tas.SingleScatteringData.from_data(params=defs)
    m = tas.ScatteringMetaData()
    m.version = flatM[0].version
    m.description = flatM[0].description
    m.source = flatM[0].source
    m.refr_index = flatM[0].refr_index
    S = []
    M = []
    for i in np.arange(new_size_grid.size):
      S.append(s)
      S[-1].abs_vec_data = a[i,...]
      S[-1].ext_mat_data = e[i,...]
      S[-1].pha_mat_data = p[i,...]
      M.append(m)
      M[-1].mass = si[i,0]
      M[-1].diameter_max = si[i,1]
      M[-1].diameter_volume_equ = si[i,2]
      M[-1].diameter_area_equ_aerodynamical = ar[i]
      setattr(M[-1],attname,new_size_grid[i])
  else:
    raise Exception( 'Non-typhon *SingleScatteringData* not yet implemented' )
  
  return S,M



def assp_create_mix(S0,M0,W):
  #2017-04-12 Jana Mendrok
  """Sums up scattering data.
  
  This function can be used to create scattering data for a habit mix. The data
  are summed up according to provided weights. 

  The input arguments *S0*, *M0* and *W* must have the same size. Each row is
  assumed to correspond to a habit, and each column to a size. That is, the
  scattering data for habit i and size j are S0[i][j] and all *S0* in a column
  should have common size.

  Each S0[i][j] is given weight W[i][j]. Accordingly, each column of *W* should
  sum up to 1. Diameters and masses in *M0* are given the same weights, to
  create the matching data in final *M*.

  No interpolations are performed, and all data for each size must have
  common grids and have the same *ptype*. But note that grids and *ptype* are
  allowed to differ between sizes.

  Parameters
  ----------
  S0: 2-D list or array of objects (of SingleScatteringData)
    ARTS-type single scattering data for collections of habits.
  M0: 2-D list or array of objects (of ScatteringMetaData)
    ARTS-type scattering meta data corresponding to S0.
  W: 2-D list or array of float
    Weight matrix to apply on S0 and M0 characterising the habit mix.
  
  Returns
  -------
  S: object or list of objects (of SingleScatteringData)
    ARTS-type single scattering data of one or several scattering elements.
  M: object or list of objects (of ScatteringMetaData)
    ARTS-type scattering meta data of one or several scattering elements.
  """
  # Basic input checks
  S0 = np.array(S0)
  if not( len(S0.shape)==2 ):
    raise Exception( \
      'S0 does not seem to be a regular 2D object (2D-array or regular 2-level nested list).' )
  M0 = np.array(M0)
  if not( len(M0.shape)==2 ):
    raise Exception( \
      'M0 does not seem to be a regular 2D object (2D-array or regular 2-level nested list).' )
  W = np.array(W)
  if not( len(W.shape)==2 ):
    raise Exception( \
      'W does not seem to be a regular 2D object (2D-array or regular 2-level nested list).' )

  nhabits = S0.shape[0]
  if not( nhabits>1 ):
    raise Exception( \
      'Your data contain only one habit (or none at all). Why use this function?' )
  nsizes = S0.shape[1]
  if not( nsizes>1 ):
    raise Exception( \
      'Your data seems empty (at least no scattering elements found in habit #1).' )
  if not( M0.shape==S0.shape ):
    raise Exception( '*S0* and *M0* must have the same size.' )
  if not( W==S0.shape ):
    raise Exception( '*S0* and *W* must have the same size.' )
  if not( np.abs(W.sum(axis=0)-1.) < 1e-3 ):
    raise Exception( 'The sum of weights over each column in *W* shall be 1.' )

  # Check and harmonise content
  for j in np.arange(nsizes):
    ptype = np.zeros(nhabits)
    for i in np.arange(nhabits):
      if not( S0[0,j].f_grid==S0[i,j].f_grid ):
        raise Exception( \
          'The frequency grid of all habits in one column of *S0* must be the same.' )
      if not( S0[0,j].T_grid==S0[i,j].T_grid ):
        raise Exception( \
          'The temperature grid of all habits in one column of *S0* must be the same.' )
      if not( S0[0,j].za_grid==S0[i,j].za_grid ):
        raise Exception( \
          'The zenith angle grid of all habits in one column of *S0* must be the same.' )
      if not( S0[0,j].aa_grid==S0[i,j].aa_grid ):
        raise Exception( \
          'The azimuth angle grid of all habits in one column of *S0* must be the same.' )

      try:
        ptype[i] = ai_ptype[S0[i,j].ptype]
      except KeyError:
        raise Exception( \
          'Unknow *ptype* found in S0[%i,%i]: %s' %(i,j,S0[i,j].ptype) )
  
    # Do we have a mix of ptypes?
    if not( np.unique(ptype).size==1 ):
      raise Exception( 'Mixed ptypes are not yet handled.' )

  # Create data for weighted mx
  # Init output argument
  from copy import deepcopy
  S = deepcopy(S0[0,:])
  M = deepcopy(M0[0,:])

  for j in np.arange(nsizes):
    S[j].description     = 'A habit mix'
    S[j].abs_vec_data *= 0.
    S[j].ext_mat_data *= 0.
    S[j].pha_mat_data *= 0.

    M[j].description     = 'A habit mix'
    M[j].source          = ''
    M[j].refr_index      = ''
    M[j].mass                            = 0.
    M[j].diameter_max                    = 0.
    M[j].diameter_volume_equ             = 0.
    M[j].diameter_area_equ_aerodynamical = 0.
  
    for i in np.arange(nhabits):
      S[j].abs_vec_data += W[i,j] * S0[i,j].abs_vec_data
      S[j].ext_mat_data += W[i,j] * S0[i,j].ext_mat_data
      S[j].pha_mat_data += W[i,j] * S0[i,j].pha_mat_data

      M[j].mass                            += W[i,j] * M0[i,j].mass
      M[j].diameter_max                    += W[i,j] * M0[i,j].diameter_max
      M[j].diameter_volume_equ             += W[i,j] * M0[i,j].diameter_volume_equ
      M[j].diameter_area_equ_aerodynamical += W[i,j] * M0[i,j].diameter_area_equ_aerodynamical

  return S,M


def assp_bulkprops(S, M, size_unit, psd, *args):
  #2017-10-30 Jana Mendrok: ported from Matlab interface.
  """Calculates bulk properties for given PSD
  
  Particle properties are specified by lists/arrays of single scattering data
  and corresponding scattering meta data. These data can be provided unsorted in
  size.

  The quantity to consider as size is selected by *size_unit*. The options are:
    'dveq' : Size is taken from M.diameter_volume_equ
    'dmax' : Size is taken from M.diameter_max
    'area' : Size is taken from M.diameter_area_equ_aerodynamical
    'mass' : Size is taken from M.mass

  The particle size distribution (PSD) is specified by a function handle. The
  function selected shall take size as first input, and can have an arbitrary
  number of additional input arguments. The PSD arguments are given after *psd*.
  For example, if an exponential PSD exists having the format
     n = exp_psd(x,n0,la)
  this function could be called as
     B = assp_bulkprops(S,M,'dveq',exp_psd,n0,la)
  Another example, with PSD defined as an anonymous function 
     f = lambda x,n0,la: n0*np.exp(-la*x)
     B = assp_bulkprops(S,M,'dmax',f,1e6,2000);

  Of course, the PSD function and arguments selected must be consistent
  with *size_unit*. This consistency is up to the user, it can not be
  checked by the function.

  The derived bulk properties are returned using the SingleScatteringData
  format. The only difference to monodisperse single scattering data is that
  the extinction, absorption and scattering data now are not in terms of
  cross-sections anymore, but are integrated bulk properties in terms of
  ext/abs/sca coefficient (i.e. extinction per meter for B.ext_mat_data).

  No interpolations are performed, and all data for each individual size must
  have common grids and the same *ptype*.
  (Note: currently, only *ptype* 'totally_random' allowed).
  
  Parameters
  ----------
  S: 1-D list or array of objects (of SingleScatteringData)
    ARTS-type single scattering data for collections of habits.
  M: 1-D list or array of objects (of ScatteringMetaData)
    ARTS-type scattering meta data corresponding to S0.
  size_unit: str
    See above.
  psd: function handle
    See above.
  *args:
    PSD input arguments (arbitray number and type). Needs to be consistent with
    the function the psd handle points to.
  ...
  
  Returns
  -------
  B: object (of SingleScatteringData)
    Bulk scattering properties, in the format of SingleScatteringData.
  """
  # for simplicity, we convert S and M to numpy arrays.
  if not use_typhon:
    raise Exception( \
      'Non-typhon handling of *SingleScatteringData* and *ScatteringMetaData* not (yet?) implemented.' )

  try:
    s = np.array(S)
    m = np.array(M)
  except:
    raise Exception( \
      'Something is wrong with S or M. Failed to convert to arrays.' )

  if ( s[0].ptype!='totally_random' ):
    raise Exception( \
      'This function handles so far only totally random orientation data.' )
  if ( len(s.shape)>1 or
       not all([isinstance(i,tas.SingleScatteringData) for i in s]) ):
    raise Exception( 'S must be a 1D list or array of *SingleScatteringData*.' )
  if ( len(m.shape)>1 or
       not all([isinstance(i,tas.ScatteringMetaData) for i in m]) ):
    raise Exception( 'M must be a 1D list or array of *ScatteringMetaData*.' )
  if ( s.size != m.size ):
    raise Exception( 'S and M must have the same length.' )
  if ( s.size < 2 ):
    raise Exception( 'Length of S must be > 1.' )
  if ( not isinstance(size_unit, str) ):
    raise Exception( '*size_unit* must be a string.' )

  # Size grid
  if ( size_unit.lower()=='dveq' ):
    x = [i.diameter_volume_equ for i in M]
  elif ( size_unit.lower()=='dmax' ):
    x = [i.diameter_max for i in M]
  elif ( size_unit.lower()=='area' ):
    x = [i.diameter_area_equ_aerodynamical for i in M]
  elif ( size_unit.lower()=='mass' ):
    x = [i.mass for i in M]
  else:
    raise Exception( \
      "Valid options for *size_unit* are: 'dveq', 'dmax', 'area' and 'mass'." )
  if ( any(np.isnan(x)) ):
    raise Exception( \
      'There is at least one Nan in data matching size_unit = %s.' %size_unit )
  x = np.array(x)

  # Sort in size (otherwise trapz will not work properly)
  ind = x.argsort()
  x = x[ind] 
  s = s[ind]

  # Loop sizes and compile data
  ns  = s.size
  nf  = s[0].f_grid.size
  nt  = s[0].T_grid.size
  nza = s[0].za_grid.size

  abs_vec = np.zeros( (ns, nf, nt) )
  ext_mat = np.zeros( (ns, nf, nt) )
  pha_mat = np.zeros( (ns, nf, nt, nza, 1, 1, 1, 6) )

  for i in np.arange(ns):
    abs_vec[i,...] = s[i].abs_vec_data.squeeze()
    ext_mat[i,...] = s[i].ext_mat_data.squeeze()
    pha_mat[i,...] = s[i].pha_mat_data

  if ( np.isnan(abs_vec).any() ):
    raise Exception( 'There is at least one Nan in S.abs_vec_data.' )
  if ( np.isnan(ext_mat).any() ):
    raise Exception( 'There is at least one Nan in S.ext_mat_data.' )
  if ( np.isnan(pha_mat).any() ):
    raise Exception( 'There is at least one Nan in S.pha_mat_data.' )

  # Get PSD
  n = psd( x, *args );

  # Init output argument
  B = s[0]

  # Calculate bulk properties
  B.abs_vec_data[...,0,0,0] = np.trapz( np.multiply(abs_vec.T, n).T, x=x, axis=0 )
  B.ext_mat_data[...,0,0,0] = np.trapz( np.multiply(ext_mat.T, n).T, x=x, axis=0 )
  B.pha_mat_data = np.trapz( np.multiply(pha_mat.T, n).T, x=x, axis=0 )

  return B


def assp_save(S,M,filename):
  #2017-03-27 Jana Mendrok
  """Stores a combination of S and M from a (pickled) file.
  
  Parameters
  ----------
  S: 1-D array of SingleScatteringData
    Array of ARTS-type single scattering data.
  M: 1-D array of ScatteringMetaData
    Array of ARTS-type scattering meta data.
  filename: str
    Full path and file name of file to create.
  
  Returns
  -------
  None
  """
  try:
    import pickle
  except ImportError:
    raise Exception( 'Module *pickle* required but not found. Can not save data.' )

  with open(filename, 'wb') as f:
    pickle.dump([S,M], f)


def assp_load(filename):
  #2017-03-27 Jana Mendrok
  """Loads a combination of S and M from a (pickled) file.
  
  Parameters
  ----------
  filename: str
    Full path and file name of file to load.
  
  Returns
  -------
  S: 1-D array of SingleScatteringData
    Array of ARTS-type single scattering data.
  M: 1-D array of ScatteringMetaData
    Array of ARTS-type scattering meta data.
  """
  try:
    import pickle
  except ImportError:
    raise Exception( 'Module *pickle* required but not found. Can not load data.' )

  try:
    import os.path
  except ImportError:
    raise Exception( 'Module *os.path* required but not found. Can not load data.' )

  if os.path.isfile(filename):
    with open(filename, 'rb') as f:
      S,M = pickle.load(f)
  else:
    raise Exception( '%s does not exist. Can not load data.' %filename )

  return S,M



def flatappend_from_recursive_S(S,flat,attname,count):
  '''
  Note: *flat* can be initialised as numpy array or list.
  '''
  if use_typhon:
    if (isinstance(S,tas.SingleScatteringData)):
      flat = np.append(flat, getattr(S,attname))
      count += 1
      return flat, count
    else:
      try:
        for s in S:
          flat,count = flatappend_from_recursive_S(s,flat,attname,count)
      except TypeError:
        raise Exception( \
          'At least one final tree member is not of SingleScatteringData. Can not handle this.' )
      return flat, count
  else:
    raise Exception( 'Non-typhon *SingleScatteringData* not yet implemented' )

def shallowappend_from_recursive_S(S,shallow,attname,instancetype='SSD'):
  '''
  Note: *shallow* must be initialised as list.
  '''
  if not( isinstance(shallow,list) ):
    raise Exception( \
      'Output shallow list needs to be initialized as list (but is not).' )
  if not( instancetype in ['SSD','SMD'] ):
    raise Exception( "Wrong instancetype. Can only handle 'SSD' and 'SMD'." )
  if use_typhon:
    if ( (instancetype=='SSD' and isinstance(S,tas.SingleScatteringData)) or
         (instancetype=='SMD' and isinstance(S,tas.ScatteringMetaData)) ):
      if (attname=='ALL'):
        shallow.append(S)
      else:
        shallow.append([getattr(S,attname)])
      return shallow
    else:
      try:
        for s in S:
          shallow = shallowappend_from_recursive_S(s,shallow,attname,instancetype)
      except TypeError:
        if instancetype=='SSD':
          raise Exception( \
            'At least one final tree member is not of SingleScatteringData. Can not handle this.' )
        else:
          raise Exception( \
            'At least one final tree member is not of ScatteringMetaData. Can not handle this.' )
      return shallow
  else:
    raise Exception( 'Non-typhon *SingleScatteringData* not yet implemented' )



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
