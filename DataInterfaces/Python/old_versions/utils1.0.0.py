# -*- coding: utf-8 -*-

"""Collection of utility functions to access ARTS single scattering database."""


try:
  import os
except ImportError:
  raise Exception( 'Module *os* required but not found.' )

try:
  import netCDF4
except ImportError:
  raise Exception( 'Module *netCDF4* required but not found.' )

try:
  import numpy as np
except ImportError:
  raise Exception( 'Module *numpy* required but not found.' )

# for debugging
#import pdb



def ssdb_read_ncfile(nc_file, freq_range=None, temp_range=None, bool_onlyFolders=False):
  #2016-11-10 Robin Ekelund
  #2017-10-27 Robin Ekelund: Modified according to changes in data and folder format.
  #2017-10-30 Robin Ekelund: Fixed bug and changed output for "folders only"-option.
  """Reads input netCDF file for data.

  Parameters
  ----------
  nc_file: netCDF file name (incl. path).
  freq_range: 2-element list
    Frequency limits. Ignore data outside these limits. If not given, no limits
    are applied, i.e. all available data is returned.
  temp range: 2-element list
    Temperature limits. Ignore data outside these limits. If not given, no limits
    are applied, i.e. all available data is returned.
  bool_onlyFolders: Boolean. If true, will only extract frequency, temperatures
    and nodata of internal folders.
  Returns
  -------
  SSD: 2D numpy array
    Single scattering data. Dimensions are 'freq x temp".
  freq: list
    Frequencies corresponding to SSD rows.
  temp: list
    Temperatures corresponding to columns in SSD.
  nodata: numpy array
    Boolean vector/matrix, flagging empty SSD elements.
    """

  def get_ncData(dataSet): # Recurse data extraction
    data = {}
    # Extracts attributes
    for attribute in dataSet.ncattrs():
      data[attribute] = getattr(dataSet, attribute)
    # Extract variables
    for variable in dataSet.variables.keys():
      tmp_dict = {'data': dataSet.variables[variable][:], #.data,
                  'dim': dataSet.variables[variable].dimensions}
      if isinstance(tmp_dict['data'],np.ma.MaskedArray):
        tmp_dict['data'] = tmp_dict['data'].data
      data[variable] = tmp_dict
    # Extract groups
    for group in dataSet.groups.keys():
      data[group] = get_ncData(dataSet.groups[group])
    return data

  # Looping over the internal netCDF4 data structure
  def get_internalStruct(dataSet, freq_rang, temp_range, bool_onlyFolders): 
    folders = []
    for grp_key in dataSet.groups.keys():
      freq = 1e9*np.float(grp_key.partition('Freq')[-1].partition('GHz')[0])
      temp= np.float(grp_key.partition('T')[-1].partition('K')[0]) 
      if (freq<freq_range[0]) or (freq>freq_range[1]):
        continue
      if (temp<temp_range[0]) or (temp>temp_range[1]):
        continue
      group = dataSet.groups[grp_key]
      if bool_onlyFolders:
        data = {'freq':freq, 'temp':temp, 'name': grp_key}
      else:
        data = get_ncData(group)
      folders.append({'name': grp_key, 'data':data, 'freq':freq, 'temp':temp})
    return folders

  def order_data(folders):
    freq = [x['freq'] for x in folders]
    temp = [x['temp'] for x in folders]
    freq_ref = list(set(freq))
    freq_ref.sort()
    temp_ref = list(set(temp))
    temp_ref.sort()
    n_f = len(freq_ref)
    n_t = len(temp_ref)
    SSD = np.empty([n_f, n_t], dict)
    nodata = np.full([n_f, n_t], True)
    for i in range(0, n_f):
      for n in range(0, n_t):
        I = [x==freq_ref[i] and y==temp_ref[n] for (x,y) in zip(freq, temp)]
        if I.count(True) == 1:
          k = I.index(True)
          nodata[i, n] = False
          SSD[i, n] = folders[k]['data']
    return SSD, np.array(freq_ref), np.array(temp_ref), nodata
  
  # Pre-processing
  if (freq_range is not None):
    if not( (len(freq_range)==2) and (float(freq_range[0])==freq_range[0]) and
           (float(freq_range[1])==freq_range[1]) ):
      raise Exception( '*freq_range* must be a numeric vector of length 2.' )
  else:
    freq_range = [0,np.inf]

  if (temp_range is not None):
    if not( (len(temp_range)==2) and (float(temp_range[0])==temp_range[0]) and
           (float(temp_range[1])==temp_range[1]) ):
      raise Exception( '*temp_range* must be a numeric vector of length 2.' )
  else:
    temp_range = [0,np.inf]

  # Read data
  if not os.path.exists(nc_file):
    raise Exception( 'NetCDF data file does not exist:\n%s\n' %nc_file )
  data = netCDF4.Dataset(nc_file, 'r', 'NETCDF4')
  folders = get_internalStruct(data, freq_range, temp_range, bool_onlyFolders)
  
  # Organize data
  if not folders:
    SSD = np.array([], dict)
    freq = np.array([])
    temp = np.array([])
    nodata = np.array([], bool)
  else:
    SSD, freq, temp, nodata = order_data(folders)
  data.close()
  return SSD, freq, temp, nodata
  
def ssdb_data_files(data_folder):
  #2016-11-16 Jana Mendrok
  #2017-10-30 Jana Mendrok: Adapt for changes in data folder structure. Following Matlab interface changes.
  """Finds data files in given folder.
  
  Explores the actual data of a given data folder.
  
  Parameters
  ----------
  data_folder: str
    Full path to a data folder.

  Returns
  -------
  data_file: list of str
    Full path and file name for each data file found.
  dmax: numpy array
    Dmax, according to file name.
  dveq: numpy array
    Dveq, according to file name.
  mass: numpy array
    Mass, according to file name.
  """
  content = os.listdir( data_folder )
  ihit = [x for x in content if not os.path.isdir(os.path.join(data_folder,x))
          and x[:4]=='Dmax']
  if not (len(ihit)!=0):
    raise Exception( \
      'Not a single data file was found. Have you really provided a data folder?.' )
  
  n = len(ihit)
  data_file = []
  dmax = np.zeros(n)
  dveq = np.zeros(n)
  mass = np.zeros(n)
  for i,x in enumerate(ihit):
    data_file.append(os.path.join(data_folder,x))
    dmax[i] = round_sig(1e-6*np.float(x.partition('Dmax')[-1].partition('um')[0]),10)
    dveq[i] = round_sig(1e-6*np.float(x.partition('Dveq')[-1].partition('um')[0]),10)
    mass[i] = round_sig(np.float(x.partition('Mass')[-1].partition('kg')[0]),10)

  return data_file,dmax,dveq,mass


def ssdb_habits(habit_id, printinfo=False):
  #2016-11-17 Jana Mendrok
  #2017-10-30 Jana Mendrok: Adapted according to Matlab interface changes regarding new folder structure.
  """Keeps track on basic habit information.

  The function keeps track of defined habit id-s, and the folder for each
  habit. This information is stored by two internal (static) information
  variables. To create these information variables, call the function with the
  path to the database top folder (ie with a string parameter instead of an
  integer parameter):
  > ssdb_habits( topfolder )
  
  To obtain the folder for a habit use:
  > folder = ssdb_habits( habit_id )
  
  If the specified habit_id is not defined, there is an error.
  
  If printinfo==True, information on the habit_id's folder is printed on screen,
  otherwise the full path to the habit folder is returned.
  
  Setting habit id to -1 has a special meaning. In this case, the function
  returns all defined habit id-s as a vector.
  
  Parameters
  ----------
  habit_id: int or str
    if int: ID of the habit to extract basic information for.
    if str: Full path to top folder of the database.
  printinfo: bool
    Flag whether to screen print habit folder information (True) or to return
    the habit folder path as string (False).

  Returns
  -------
  folder: str
    Full path to folder for specified habit.
  """
  
  # some local helper functions
  def isint(x):
    try:
      ii = (x==int(x))
    except:
      ii = False
    return ii
  
  def scan_topfolder( topfolder ):
    main_folders = ['TotallyRandom', 'AzimuthallyRandom']
    print_folders = ''
    for s in main_folders:
      print_folders += (s+'\n')
    
    # Check that topfolder contains at least one of the "main folders"
    if not( any(os.path.isdir(os.path.join(topfolder,i)) for i in main_folders) ):
      raise Exception( \
        'There should at least one of the following folders be in *topfolder*:\n' + \
        print_folders )
    
    # Scan the main folders
    hs = []
    folders = []
    for i in main_folders:
      nextfolder = os.path.join(topfolder,i)
      if (os.path.isdir(nextfolder)):
        hs,folders = scan_folder( nextfolder, hs, folders )
    
    return hs,folders
  
  def scan_folder( datafolder, hs, folders ):
    content = os.listdir( datafolder )
    
    # Use report.pdf as indicator if we are in a habit folder, or higher up in
    # the folder tree
    if 'report.pdf' in content:
      # We have reached a habit folder
      summary_file = ssdb_summary_find( datafolder )
      S = ssdb_summary_read( summary_file )
      hs.append( S['HABIT_IDENT'] )
      folders.append( datafolder )
    else:
      # We have to continue downwards
      for i in content:
        # Only folders to be checked, any possible file and . and .. should be ignored
        nextfolder = os.path.join(datafolder,i)
        if os.path.isdir(nextfolder):
          hs,folders = scan_folder( nextfolder, hs, folders )
    
    return hs,folders


  #- Fill variables linking habits with a folder
  if isinstance(habit_id,str):
    # Scan the database folder structures (local sub-function)
    hs,folders = scan_topfolder( habit_id )
    if not(len(hs)>0):
      raise Exception( 'Not a single folder was found. Something must be wrong.' )
    
    # Save static variables
    ssdb_habits.habits = hs*1         # in py2.* .copy() not working on lists
    ssdb_habits.folders = folders*1   # hence, clone them by pseudo-modification
    
    # Check that no id numbers is duplicated
    n = len( hs )
    hs,counts = np.unique( ssdb_habits.habits, return_counts=True )
    if len(hs) < n:
      for i in np.arange(len(hs)):
        if counts[i] > 1:
         print( 'Habit id %d is used %d times.\n' %(hs[i], counts[i]) )
      assert(len(hs)==n), \
        'One or several habit numbers are duplicated. See above.'

  #- Extract folder for specified habit
  elif isint(habit_id):
    errmsg = 'The function has not yet been initialised. Do this by calling ' + \
             'ssdb_habits (or ssdb_init) with the top folder path as input.'
    if habit_id==-1:
      try:
        folder = ssdb_habits.habits
      except:
        raise Exception( errmsg )
    else:
      try:
        found = habit_id in ssdb_habits.habits
      except:
        raise Exception( errmsg )
      if not(found):
        raise Exception( 'There is no habit with id = %d.' %habit_id )
      folder = ssdb_habits.folders[ssdb_habits.habits.index(habit_id)]
    
    if printinfo:
      print( 'Habit %d is found in folder:\n%s\n' %(habit_id, folder) )
    else:
      return folder

  else:
    raise Exception( 'Invalid input argument.' )


def ssdb_display(habit_id=None, orientation=None, dveq=None):
  #2016-11-20 Jana Mendrok
  #2016-12-06 Jana Mendrok: adapted output formating to Matlab output changes.
  #2017-10-30 Jana Mendrok: Adapted according to changes in Matlab interface regarding changed folder and netcdf file structure.
  """Explore the database content.
  
  The function summarises the database content on the screen. Different levels
  of the database are displayed depending on which input arguments are
  specified. Three format versions ar at hand:
  
  To get an overview of habits and orientations at hand:  
  > ssdb_display()
  
  To list the particle sizes at hand for one habit and orientation combination:
  > ssdb_display( habit_id=integer, orientation=string )
  
  To list the frequencies and temperature combinations at hand for one habit,
  orientation and size combination:
  > ssdb_display(  habit_id=integer, orientation=string, dveq=float )
  
  Parameters
  ----------
  habit_id: int
    ID of the habit to explore.
  orientation: str
    descriptor of orientation to explore for the chosen habit.
  dveq: float (or int)
    size (in terms of volume equivalent sphere diameter) of the particle to
    explore [units: m].

  Returns
  -------
  None
  """
  # nothing specified: derive and print what habits and orientations are available
  if habit_id is None:
    if not(orientation is None):
      raise Exception( 'Specify either habit_id AND orientation or none of the two.' )
    
    main_folders = [ 'Ice', 'Mixed', 'Liquid' ]
    
    # Get all habit id-s
    habits = ssdb_habits( -1 )
    
    # Obtain folders names and extract the four levels for each habit
    fullpaths = ssdb_habits( habits[0] )
    topfolder = ''
    for i in main_folders:
      if i in fullpaths:
        topfolder = fullpaths.partition(i)[0]
    assert(len(topfolder)>0), 'Failed to identify topfolder.'
    fullpaths = [ fullpaths.partition(topfolder)[-1] ]
    
    for i in habits[1:]:
      fullpaths.append( ssdb_habits(i).partition(topfolder)[-1] )
    
    sorthab = np.argsort(fullpaths)
    
    #Print top lines
    print( '' )
    print( '---------------------- ARTS single scattering data ----------------------')
    print( '' )
    print( '  The table below gives an overview of habits at hand. The first three' )
    print( '  levels give a rough classification of each habit. The information found')
    print( '  in the two last levels is: ' )
    print( '' )
    print( '        habit id: habit name                           a         / b ' )
    print( '          orientations' )
    print( '' )
    print( '  where mass = a * Dmax^b.' )
    print( '' )
    print( '  Habits marked with "(x)" are still work in progress. Use with caution.' )
    print( '' )
    
    bs = ' '
    partpaths = fullpaths[sorthab[0]].split(os.sep)
    prevpaths = len(partpaths)*['']
    curhab = 0
    while curhab<len(sorthab):
      # derive changes between old and new habit folder
      partpaths = fullpaths[sorthab[curhab]].split(os.sep)
      identlev = 0
      ident = True
      while ident and (identlev<min(len(partpaths),len(prevpaths))):
        if partpaths[identlev]==prevpaths[identlev]:
          identlev += 1
        else:
          ident = False
      # print header
      for i in np.arange(identlev,len(partpaths)-1):
        if i==0:
          print( '%s\n---' %partpaths[0] )
        else:
          print( '%s%s' %(i*bs,partpaths[i]) )
      # derive and print habit summary
      hab_id = habits[sorthab[curhab]]
      S = ssdb_summary( hab_id, printinfo=False )
      try:
        if S['STATUS']=='WORKING':
          status_symbol = '(x)'
        elif S['STATUS']=='COMPLETED':
          status_symbol = ''
        else:
          status_symbol = '(?)'
      except:
        status_symbol='?'
      print( '   %3s %3d: %s%sa=%.1e / b=%.2f' \
        %( status_symbol, hab_id, S['HABIT_NAME'],
           (43-len(S['HABIT_NAME']))*bs,
           float(S['ALPHA']), float(S['BETA']) ) )
      # derive and print orientation info
      ofolders,orientations,tilt_angles = \
        ssdb_particle_folders( os.path.join(topfolder,fullpaths[sorthab[curhab]]) )
      for o in np.arange(len(orientations)):
        if (tilt_angles[o] is None):
          print( '          %s' %orientations[o] )
        else:
          print( '          %s, tilt=%0.3f deg' %(orientations[o],tilt_angles[o]) )

      # prepare next loop
      curhab += 1
      prevpaths = partpaths

    print( '' )
    print( '-------------------------------------------------------------------------')
    print( '' )

  # habit & orientation specified: derive and print info sizes available for this habit-orientation-combi
  else:
    if not(orientation is not None):
      raise Exception( 'Specify either habit_id AND orientation or none of the two.' )

    # Get the basic information by combining existing functions  
    S = ssdb_summary( habit_id, printinfo=False )
    ofolders,orientations,tilt_angles = \
      ssdb_particle_folders( ssdb_habits(habit_id) )

    # Filter out selected orientation
    ind = [i for i,x in enumerate(orientations) if orientation in x]
    assert(len(ind)<2), \
      'More than one matching orientation found. Should not happen.\n' + \
      'Complain to interface developers!'
    if not(len(ind)>0):
      raise Exception( \
        'No data could be located for specified orientation (%s).' %orientation )
    # Crop data according to orientation
    ofolder = ofolders[ind[0]]
    
    data_files,dmax_vec,dveq_vec,mass_vec = ssdb_data_files( ofolder )

    # Habit + orientation
    if dveq is None:
      print( '\n--- %d: %s - %s ---\n' %(habit_id, S['HABIT_NAME'], orientation) )
      print( '    a = %.1e / b = %.2f' %(float(S['ALPHA']), float(S['BETA'])) )
      print( '    N = %d\n' %len(dmax_vec) );
      print( '    Dveq [um]   Dmax [um]   Mass [g]' )
      print( '    ---------   ---------   --------' )
      
      sortd = np.argsort(dveq_vec)
      for i in np.arange(len(dmax_vec)):
        print( '%10.0f%12.0f%14.2e' %(dveq_vec[sortd[i]]*1e6, dmax_vec[sortd[i]]*1e6, mass_vec[sortd[i]]*1e3) )
      print()
    # Habit + orientation + size
    else:
      id = abs( dveq - dveq_vec ).argmin()
      if not(abs( dveq - dveq_vec )[id]<=0.5e-6):
        raise Exception( \
          'No particle found close to the given Dveq value. Allowed difference is +-0.5 um.' )

      infodata = ssdb_read_ncfile( data_files[id], bool_onlyFolders=True)
      freq_vec = np.array([i['freq'] for i in infodata[0].reshape(-1)])
      temp_vec = np.array([i['temp'] for i in infodata[0].reshape(-1)])

      # at least in python (as far as i understand), the freq_vec is already unique.
      freq = np.unique(freq_vec);
      if ( len(freq)<1 ):
        raise Exception( \
          'No data files could be found for given habit and orientation combination!' )
      
      print( '\n--- %d: %s - %s - %.0f um ---\n'
             %(habit_id, S['HABIT_NAME'], orientation, dveq_vec[id]*1e6) )
      print( ' Frequency [GHz]   x [-]   Temperatures [K] ' )
      print( ' ---------------   -----   ----------------' )
      for f in np.arange(len(freq)):
        ii = np.where(freq_vec==freq[f])[0]
        x = np.pi * dveq_vec[id] * freq[f] / 3e8
        temp = (temp_vec[ii])
        temp.sort()
        print( '%16.3f  %6.2f' %(freq[f]/1e9, x),  ''.join(['   %.2f' %i for i in temp]) )
      print()


def ssdb_summary( habit_id, printinfo=True ):
  #2016-11-18 Jana Mendrok
  """Obtains the habit summary data.

  The functions obtains the summary data for given habit id.
  
  If printinfo==True, the summary is printed on screen, else a dictionary holding
  the summary data is returned.
  
  Parameters
  ----------
  habit_id: int
    Habit id number.
  printinfo: bool
    Flag whether to screen print habit summary (True) or to return summary as
    dict (False).

  Returns
  -------
  S: dict
    Dictionary with summary data.
  """
  habit_folder = ssdb_habits( habit_id )
  
  summary_file = ssdb_summary_find( habit_folder )
  
  S = ssdb_summary_read( summary_file )
  
  if printinfo:
    fields = [i for i in S]
    fields.sort()
    try:
      fields.insert(0, fields.pop(fields.index('HABIT_NAME')))
    except:
      pass
    try:
      fields.insert(0, fields.pop(fields.index('HABIT_IDENT')))
    except:
      pass
    
    l2col = 15
    bs=' ' #single blank
    print()
    for i in fields:
      if i=='HABIT_IDENT':
        print('%s%s: %d' %((l2col-len(i))*bs,i,S[i]))
      else:
        print('%s%s: %s' %((l2col-len(i))*bs,i,S[i]))
    print()
  else:
    return S


def ssdb_summary_find( habit_folder ):
  #2016-11-17 Jana Mendrok
  """Locates the summary file inside a habit folder.

  The summary file must be placed in *habit_folder*. Sub-folders are not
  searched for any summary file.
  
  Parameters
  ----------
  habit_folder: str
    Full path to a habit folder.

  Returns
  -------
  summary_file: str
    Full path and file name of the summary file.
  """
  content = os.listdir( habit_folder )
  i = [x for x in content if 'DataSummary.' in x]
  if (len(i)==0):
    print('No data summary file of pattern "DataSummary.*" found in folder %s.' %habit_folder)
    print('Trying for obsolete pattern "DataSummary_*".')
    i = [x for x in content if 'DataSummary_' in x]
  if not(len(i)>0):
    raise Exception( \
      'No data summary file could be found in folder %s' %habit_folder )
  if not(len(i)<2):
    raise Exception( \
      'Multiple data summary files were found in folder %s' %habit_folder )
  return os.path.join( habit_folder, i[0])


def ssdb_summary_read( summary_file ):
  #2016-11-17 Jana Mendrok
  """Reads a habit summary file.

  The functions performs the actual reading of an identified summary file.
  
  The data are returned as a dictionary. The habit id is converted to a number,
  while all other fields are kept as the strings extracted from the file.
  
  Parameters
  ----------
  summary_file: str
    Full path and file name of a summary file.

  Returns
  -------
  S: dict
    Summary data dictionary.
  """
  if not(isinstance(summary_file,str)):
    raise Exception( \
      'Provided summary_file variable is not string. ' + \
      'Do not provide a list or array of string.' )
  if not(summary_file[-4:]=='.txt'):
    raise Exception( 'Data summary files are expected to have extension .txt.' )
  
  try:
    f = open(summary_file,'r')
  except:
    raise Exception( 'Could not open %s for reading.' %summary_file )
  
  S = {}
  
  line = f.readline().partition('=')
  while len(line[0])>0:
    if not(len(line[1])>0):
      raise Exception( \
        'Each line in a summary must contain at least one "=". ' + \
        'This is not the case for this line\n%s\nfound in file\n%s' \
        %(line, summary_file) )
    S[line[0]] = line[2].rstrip('\n')
    line = f.readline().partition('=')
  
  try:
    f.close()
  except:
    raise Exception( 'Could not close %s after reading.' %summary_file )
  
  try:
    id = S['HABIT_IDENT']
  except:
    raise Exception( 'No setting of HABIT_IDENT was found in %s' %summary_file )
  try:
    S['HABIT_IDENT'] = int(id)
  except:
    raise Exception( 'Invalid setting of HABIT_IDENT in %s' %summary_file )
  
  return S


def ssdb_particle_folders( habit_folder ):
  #2016-11-20 Jana Mendrok
  #2017-10-30 Jana Mendrok: Adapted to Matlab interface changes. Now looks for orientation folders.
  #2017-10-31 Robin Ekelund: Modified according to Matlab changes.
  """Finds folders for individual particles.

  Explores the particle data of a given habit folder.
  
  Parameters
  ----------
  habit_folder: str
    Full path to a habit folder.

  Returns
  -------
  orientation_folder: list of str
    Full path to each orientation folder found.
  orientations: list of str
    Orientation type (totally or azimuthally random) associated with each
    orientation folder.
  tilt_angles: list of float
    Tilt angle associated with each orientation folder.
  """
  content = os.listdir( habit_folder )
  ihit = [x for x in content if os.path.isdir(os.path.join(habit_folder,x))]
  if not (len(ihit)!=0):
    raise Exception( \
      'Not a single orientation folder found in %s.\n' + \
      'Have you really provided a habit folder?' %habit_folder )
  
  n = len(ihit)
  orientation_folder = []
  orientations = []
  tilt_angles = []
  
  for i,x in enumerate(ihit):
    orientation_folder.append(os.path.join(habit_folder,x))
    if ( 'TotallyRandom' in x ):
      orientations.append('totally_random')
      tilt_angles.append(None)
    elif ( 'AzimuthallyRandom' in x ):
      y = x.replace('AzimuthallyRandom', 'azimuthally_random') 
      orientations.append(y)
      tilt_angles.append(np.float(x.partition('beta')[-1].partition('deg')[0]))
    else:
      raise Exception( 'Folder %s not recognized as an orientation folder' %x )
  
  return orientation_folder, orientations, tilt_angles


def ssdb_init( topfolder ):
  #2016-11-18 Jana Mendrok
  """Initialises the SSDB interface.
  
  Assuming, the current folder is the interface folder, the function adds the
  interface folder to Python's search path and sets up the habit-folder table
  used by *ssdb_habits*.
  
  Parameters
  ----------
  topfolder: str
    Full path to top folder of the database.
  
  Returns
  -------
  None
  """
  #- Add current folder, assuming that's the interface folder, to search path
  # FIXME: do we actually need this? why/what for? and isn't this a bit of a
  #        far-fetched assumption anyways? wouldn't it make as much or more
  #        sense to add some relative path from the topfolder?
  os.path.sys.path.insert( 0, os.getcwd() )
  
  #- Create habit id table
  ssdb_habits( topfolder )


def ssdb_habit_logo( habit_id, show_image=True ):
  #2016-11-20 Jana Mendrok
  """Displays habit logo image.
  
  The function displays the image (optional) and returns the path to the logo file.
  
  Parameters
  ----------
  habit_id: int
    Habit id number.
  show_image: boolean
    Flag whether to displayed logo image or not.
  
  Returns
  -------
  logo_file: str
    Full path and filename of habit logo file.
  """
  habit_folder = ssdb_habits( habit_id )
  logo_file = os.path.join( habit_folder, 'shape_img.png' )
  if show_image:
    from matplotlib.pyplot import imread,imshow,show,axis,figure
    I = imread( logo_file )
    # make sure to not overwrite some older figure by opening a new figure
    figure()
    imshow( I )
    # switch off axes, frame, white padding
    axis('off')
    # show and return to script leaving open the figure window
    # (not tested whether this works in plain python shell or python from
    # OS shell command line)
    show(block=False)
  return logo_file


def ssdb_data_import(data_folder,
                     freq_range=None,temp_range=None,grid_sort=False):
  #2016-11-21 Jana Mendrok
  """Reads the scattering data of one particle instance folder.
  
  The reading can be restrictred in frequency and temperature by the optional
  arguments *freq_range* and *temp_range*.
  
  The default behaviour is to perform a plain reading and then returning
  *SSD* as an unsorted list of dictionary. *freq* and *temp* are in this case
  also vectors, having the same length as *SSD* (and *nodata* contains
  just False).
  
  By setting *grid_sort* to true, the scattering data are sorted into a matrix
  in terms of frequency (row dimension) and temperature (column dimension), and
  *SSD* is returned as 2-D array of dictionary. For example SSD[2,3] is the data
  for frequency 2 and temperature 3 (note: 0-indexing applies for all of them).
  *freq* and *temp* are in this case returned as the frequency and temperature
  grid of the data. If any element of SSD is unfilled, the corresponding element
  in *nodata* is set to True.
  
  Parameters
  ----------
  data_folder: str
    Full path to folder holding scattering netcdf data files.
  freq_range: 2-element list
    Frequency limits. Ignore data outside these limits. If not given, no limits
    are applied, i.e. all available data is returned.
  temp range: 2-element list
    Temperature limits. Ignore data outside these limits. If not given, no limits
    are applied, i.e. all available data is returned.
  grid_sort: boolean
    Flag whether to sort data in frequency and temperature and create matrix
    structure from them.
  
  Returns
  -------
  SSD: 1- or 2-D numpy array of dictionary
    Read single scattering data.
  freq: 1-D numpy array
    The frequency of each SSD element.
  temp: 1-D numpy array
    The temperature of each SSD element.
  nodata: numpy array
    Boolean vector/matrix, flagging empty SSD elements.
  """
  if (freq_range is not None):
    if not( (len(freq_range)==2) and (float(freq_range[0])==freq_range[0]) and
            (float(freq_range[1])==freq_range[1]) ):
      raise Exception( '*freq_range* must be a numeric vector of length 2.' )
  else:
    freq_range = [0,np.inf]

  if (temp_range is not None):
    if not( (len(temp_range)==2) and (float(temp_range[0])==temp_range[0]) and
            (float(temp_range[1])==temp_range[1]) ):
      raise Exception( '*temp_range* must be a numeric vector of length 2.' )
  else:
    temp_range = [0,np.inf]
  
  # Locate existing files and find index of ones inside given freq and temp
  # ranges
  data_files,freqs,temps = ssdb_data_files( data_folder )
  # we can do the temp and freq selection in one go, but lets separate for better error/warning info
  #ind = [i for i in np.arange(len(freqs))
  #       if (freqs[i]>=freq_range[0]) and (freqs[i]<=freq_range[1]) and
  #          (temps[i]>=temp_range[0]) and (temps[i]<=temp_range[1]) ]
  indf = [ i for i in np.arange(len(freqs))
           if (freqs[i]>=freq_range[0]) and (freqs[i]<=freq_range[1]) ]
  if ( len(indf)<1 ):
    print('  For %s,\n'
          '  no freqs in selected range (%.1f-%.1f GHz) found. Output data will be empty.'
          %(data_folder,freq_range[0]*1e-9,freq_range[1]*1e-9))
  indt = [ i for i in np.arange(len(freqs))
           if (temps[i]>=temp_range[0]) and (temps[i]<=temp_range[1]) ]
  if ( len(indt)<1 ):
    print('  For %s,\n'
          '  no temps in selected range (%.1f-%.1f K) found. Output data will be empty.'
          %(data_folder,temp_range[0],temp_range[1]))
  ind = list( set(indf).intersection(indt) )
  if ( len(ind)<1 and len(indf)*len(indt)>0 ):
    print('  For %s,\n'
          '  no data in intersection of selected freq (%.1f-%.1f GHz) and'
          ' temp (%.1f-%.1f K) ranges found. Output data will be empty.'
          %(data_folder,freq_range[0]*1e-9,freq_range[1]*1e-9,temp_range[0],temp_range[1]))
  
  # With grid sorting
  if grid_sort:
    freq = np.unique( freqs[ind] )
    temp = np.unique( temps[ind] )
    nodata = np.zeros( (len(freq),len(temp)), dtype=bool )
    
    SSD = []
    for fi,f in enumerate(freq):
      SSD.append([])
      for ti,t in enumerate(temp):
        i = np.where((freqs==f) & (temps==t) )[0]
        assert(len(i)<2), 'Duplicate of one frequency-temperature found!'
        if len(i)<1:
          nodata[fi,ti] = True
          SSD[fi].append({})
        else:
          SSD[fi].append( ssdb_read_ncfile( data_files[i[0]] ) )
    
    SSD = np.array( SSD ).reshape(freq.size,temp.size)
  
  # Without grid sorting
  else:
    freq = freqs[ind]
    temp = temps[ind]
    nodata = np.zeros( len(freq), dtype=bool )
    
    SSD = []
    for i in ind:
      SSD.append( ssdb_read_ncfile( data_files[i] ) )
    
    SSD = np.array( SSD )

  return SSD,freq,temp,nodata



#def ssdb_data_import_multiple(data_folders,
#                              freq_range=None,temp_range=None,
#                              verbosity=0):
  #2016-11-25 Jana Mendrok
  """Reads the scattering data of multiple folders.
  
  The function is an extended variant of ssdb_data_import(grid_sort=True) that
  considers more than one data folder, allowing to extract unique frequency and
  temperature grids over multiple datasets, e.g. several or all sizes of one
  habit.
  
  The reading can be restrictred in frequency and temperature by the optional
  arguments *freq_range* and *temp_range*.
  
  The scattering data are sorted into an array in terms of particle (page
  dimension), frequency (row dimension) and temperature (column dimension), and
  *SSD* is returned as 2-D array of dictionary. For example SSD[4,2,3] is the
  data for particle 4, frequency 2 and temperature 3 (note: 0-indexing applies
  for all of them).
  *freq* and *temp* are returned as the unique frequency and temperature grid
  of the combined data. If any element of SSD is unfilled, the corresponding
  element in *nodata* is set to True.
  
  Parameters
  ----------
  data_folders: list of str
    Full path to folders holding scattering netcdf data files.
  freq_range: 2-element list
    Frequency limits. Ignore data outside these limits. If not given, no limits
    are applied, i.e. all available data is returned.
  temp range: 2-element list
    Temperature limits. Ignore data outside these limits. If not given, no limits
    are applied, i.e. all available data is returned.
  
  Returns
  -------
  SSD: 3-D array of dictionary
    Read single scattering data.
  freq: 1-D numpy array
    The frequency of each SSD element.
  temp: 1-D numpy array
    The temperature of each SSD element.
  nodata: numpy array
    Boolean 3-D array, flagging empty SSD elements.
  """
"""
  #pdb.set_trace()
  
  if isinstance(data_folders,str):
    data_folder = data_folders
    data_folders = [data_folders]
  else:
    data_folder = None
  
  if (freq_range is not None):
    if not( (len(freq_range)==2) and (float(freq_range[0])==freq_range[0]) and
            (float(freq_range[1])==freq_range[1]) ):
      raise Exception( '*freq_range* must be a numeric vector of length 2.' )
  else:
    freq_range = [0,np.inf]

  if (temp_range is not None):
    if not( (len(temp_range)==2) and (float(temp_range[0])==temp_range[0]) and
            (float(temp_range[1])==temp_range[1]) ):
      raise Exception( '*temp_range* must be a numeric vector of length 2.' )
  else:
    temp_range = [0,np.inf]
  
  # Locate existing files and find index of ones inside given freq and temp
  # ranges
  data_files=[]
  freqs = []
  temps = []

  for x in data_folders:
    tfiles,tfreqs,ttemps = ssdb_data_files( x )
    data_files.append(tfiles)
    freqs.append(tfreqs)
    temps.append(ttemps)

  # flattened f and T arrays
  ffreqs = np.array([item for sublist in freqs for item in sublist])
  ftemps = np.array([item for sublist in temps for item in sublist])
  ind = [i for i in np.arange(len(ffreqs))
         if (ffreqs[i]>=freq_range[0]) and (ffreqs[i]<=freq_range[1]) and
            (ftemps[i]>=temp_range[0]) and (ftemps[i]<=temp_range[1]) ]
  
  freq = np.unique( ffreqs[ind] )
  temp = np.unique( ftemps[ind] )
  nodata = np.zeros( (len(data_folders),len(freq),len(temp)), dtype=bool )
  if verbosity>0:
    print( '%i unique frequencies and %i unique temperatures ' %(len(freq),len(temp)),
           'found in the given %i particle folders' %len(data_folders) )
    
  SSD = []
  for di,d in enumerate(data_folders):
    if verbosity>0:
      print( 'processing particle folder #%i out of %i' %(di,len(data_folders)))
    SSD.append([])
    for fi,f in enumerate(freq):
      SSD[di].append([])
      for ti,t in enumerate(temp):
        i = np.where((freqs[di]==f) & (temps[di]==t) )[0]
        assert(len(i)<2), 'Duplicate of one frequency-temperature found!'
        if len(i)<1:
          nodata[di,fi,ti] = True
          SSD[di][fi].append({})
        else:
          ssd = ssdb_read_ncfile( data_files[di][i[0]] )
          SSD[di][fi].append( ssd )
  
  if data_folder is not None:
    SSD = SSD[0]
    data_folders = data_folder
  
  SSD = np.array( SSD )
  return SSD,freq,temp,nodata
"""


def ssdb_import_habit(habit_id, orientation,
                      habit_folder=None,
                      size_range=None, size_type='dveq',
                      freq_range=None, temp_range=None,
                      allow_nodata=False):
  #2017-10-16 Jana Mendrok: Extracted from assp.assp_import_ssdb following Matlab example.
  #2017-10-16 Jana Mendrok: Added functionality for azimuthally random orientation following M.Brath's draft.
  #2017-10-30 Jana Mendrok: Adapted to Matlab interface changes related to changed folder and netcdf file structure.
  """Reads and compiles data from SSDB database for one habit.
  
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
  temp_range: 2-element list
    Temperature limits [unit: K]. Ignore data outside these limits. If not
    given, no limits are applied, i.e. all available data is considered.

  Returns
  -------
  data: list of dictionary
    Single scattering and auxiliary data per particle instances of the habit.
  """
  if (habit_folder is None):
    try:
      if not( len(habit_id)==1 ):
        raise Exception( \
          '*habit_id* only allowed to have length 1 (or be a scalar).' )
      habitID = habit_id[0]
    except: #habit_id is scalar
      habitID = habit_id
    if not( isinstance(habitID, int) ):
      raise Exception( '*habit_id* must be an integer.' )
  else:
    habitID = -1
    if not(os.path.isdir(habit_folder)):
      raise Exception( \
        'Provided habit_folder (%s) is not a valid folder.' %habit_folder )

  if not( isinstance(orientation, str) ):
    raise Exception( '*orientation* must be a string.' )

  if allow_nodata is not True:
    allow_nodata=False
  
  if (size_range is not None):
    if not( (len(size_range)==2) and (float(size_range[0])==size_range[0]) and
            (float(size_range[1])==size_range[1]) ):
      raise Exception( '*size_range* must be a numeric vector of length 2.' )
  else:
    size_range = [0,np.inf]

  if not( isinstance(size_type,str) ):
    raise Exception( '*size_type* must be a string.' )
  size_type = size_type.lower()
  if not( size_type in ['dveq', 'dmax', 'mass'] ):
    raise Exception( \
      "Valid options for *size_type* are 'dveq', 'dmax', and 'mass'." + \
      "You selected %s." %size_type )
  
  # freq_range and temp_range are checked in ssdb_read_ncfile
  
  # Make an inventory of existing data
  errmsg = 'No orientation folders could be located for specified habit (ID %i).' %habitID
  if (habit_folder is None):
    ofolders, orientations, tilt_angles = \
      ssdb_particle_folders( ssdb_habits(habitID) )
    if not (len(ofolders)!=0):
      raise Exception( errmsg )
  else:
    ofolders, orientations, tilt_angles = \
      ssdb_particle_folders( habit_folder )
    if not (len(ofolders)!=0):
      raise Exception( errmsg )
  
  # Filter out selected orientation
  ind = [i for i,x in enumerate(orientations) if orientation in x]
  assert(len(ind)<2), \
    'More than one matching orientation found. Should not happen.\n' + \
    'Complain to interface developers!'
  if (len(ind)<1):
    #orients = [item for sublist in orientations for item in sublist] #when we had nested lists...
    orients = [item for item in orientations]
    orients= np.unique( orients )
    orient_list=''
    for item in orients:
      orient_list += '   %s\n' %item
    raise Exception( \
      "No data could be located for specified orientation ('%s').\n" %orientation + \
      'The following orientations are at hand:\n' + orient_list )

  # Crop data according to orientation
  ofolder     = ofolders[ind[0]]
  tilt_angle  = tilt_angles[ind[0]]

  # Loop particle sizes and import data
  data = []
  print( 'Extracting data from SSDB.' )
  data_files,dmax,dveq,mass = ssdb_data_files( ofolder );

  # Crop data according to "size"
  if (size_type=='dveq'):
    size = dveq
  elif (size_type=='dmax'):
    size = dmax
  elif (size_type=='mass'):
    size = mass
  ind = [i for i,x in enumerate(size)
         if (x>size_range[0]) and (x<=size_range[1]) ]
  data_files = [data_files[i] for i in ind]
  size    = size[ind]
  dmax    = dmax[ind]
  dveq    = dveq[ind]
  mass    = mass[ind]
  # Sort
  #ind = size.argsort()
  #size = size [ind]
  #data_files = [data_files[i] for i in ind]
  #dmax    = dmax[ind]
  #dveq    = dveq[ind]
  #mass    = mass[ind]

  print( '%i scattering elements to process' %len(data_files) )
  for n in np.arange(len(data_files)):
    if ((n+1)%5==0):
      print( '  processing element %i' %(n+1) )
    data.append({}) # creates one entry per entry in size. remains empty in case no SSD avalable for it.
    SSD, freq, temp, nodata = \
      ssdb_read_ncfile( data_files[n], freq_range, temp_range )
    if (SSD.size>0):
      if (len(SSD.shape)>1):
        shapeinfo = SSD[0,0]['ShapeData']
      else:
        shapeinfo = SSD[0]['ShapeData']
      #data.append({}) # creates one entry per SSD-filled entry in size (ie when freq and temp available for this size)
      data[-1] = {'SSD': SSD,
                  'dveq': shapeinfo['diameter_vol_eq']['data'],
                  'dmax': shapeinfo['diameter_max']['data'],
                  'mass': shapeinfo['mass']['data'],
                  'freq': freq,
                  'temp': temp,
                  'nodata': nodata,
                  'orientation': orientations,
                  'tilt_angle': tilt_angle}

      if ( (not allow_nodata) and nodata.any() ):
        raise Exception( \
          'At least one frequency and temperature combination lacks data ' + \
          'for id %i, orientation %s and Dveq=%.2f um (Dmax = %.2f um).\n\n' \
          %(habitID, orientation, dveq[n]*1e6, dmax[n]*1e6) + \
          'Lacking data for a frequency-temperature combination, see above ' + \
          'for details. To allow lacking data, set allow_nodata=True.' )
    else:
      print('For this particle (Dveq=%.2f um, Dmax = %.2f um), no SSD remaining -'
            ' skipping this particle (ie. empty entry in data).\n'
            %(dveq[n]*1e6, dmax[n]*1e6))

  return data


def round_sig(x, sig=1):
  #2017-10-30 Jana Mendrok
  """Equivalent to Matlab's round(x, sig, 'significant') functionality.
  
  Source: https://stackoverflow.com/questions/3410976/how-to-round-a-number-to-significant-figures-in-python answered by indgar on Aug 5 '10 at 9:49
  
  Parameters
  ----------
  x: number (float, int, ...)
    The number to round.
  sig: int
    Number of significant digits to rounf to.
  
  Returns
  -------
  y: number (as x)
    x rounded to sig significant digits.
  """
  return np.round(x,sig-np.int(np.floor(np.log10(abs(x))))-1)


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
