# NOTE: 0) run this script from within (i)python with
#            >>> run demo_ssp4arts.py
#          or from terminal with
#            >>> python demo_ssp4arts.py
#       1) Requires the typhon package to be installed (and its location to be known).
#       2) To allow the script to locate your ARTS-SSDB Python interface, at least
#          one of the following needs to be fulfiled (used in this order of priority):
#          - You run this script from inside the Python interface folder.
#          - You run this script from inside the DataInterfaces folder.
#          - DataInterfaces folder (or, in general, the superfolder to the
#            interface Python folder) is in your PYTHONPATH.
#          - the DataInterfaces folder (in turn assumed to hold the Python interface
#            subfolder) is a subfolder to the specified SSDBpath.
#       3) Currently produces one scat and meta data file per scattering
#          element. In ARTS, you can use ScatSpeciesScatAndMetaRead to import
#          these data.


####################
# Input. Modify according to your needs.
####################

# location of SSDB in your system
SSDBpath = '/your/database/location/SSD'

# location (and name) of ARTS scat_data and scat_meta output
#   the filename is constructed as:
#     > outpath/sdBase_sdSize_size[.meta].sdExt
#   with
#     sdSize in ['dmax','dveq','mass']
#     size being the respective size of the particle in [um] or [kg]
#     sdExt in ['xml','xml.gz','xml.bin'] determining the file format
outpath = '.'
sdBase = 'EvansSnow_TotRand'
sdExt = 'xml.bin'
sdSize = 'dmax'

# properties of particles to extract SSP for
#   use utils.ssdb_display to see what data (habits, orientations, sizes,
#   frequencies, temperatures) are available.
#   see assp.assp_import_ssdb for further details on the parameters.

habID = 1 # habit ID
orient = 'totally_random' # orientation

minD = 0. # minimum size to extract
maxD = float("inf") # maximum size to extract
sizeparam = 'dmax' # size parameter to apply extraction limits on.
                    # available (param [unit]): 'dmax' [m], 'dveq' [m], 'mass' [kg]

fmin = 0. # minimum frequency to extract [Hz]
fmax = float("inf") # maximum frequency to extract

tmin = 0. # minimum temperature to extract [K]
tmax = float("inf") # maximum temperature to extract



####################
# Importing required modules. Do not modify (unless you know what you are doing).
####################
#try import from current location
#(ie see, whether we are inside the Python interface folder)
try:
  import utils
  import assp
except ImportError:
  #try import of Python interface package from current location (by inserting
  # current location at start of sys.path)
  # (ie see, whether we are inside a DataInterfaces folder)
  #implicitly also imports (with lower priority) if Python interface package is
  # a subfolder to any of the other PYTHONPATH entries.
  try:
    import os.path
    os.path.sys.path.insert(0,'')
    from Python import utils
    from Python import assp
  except ImportError:
    #finally, try DataInterfaces folder in SSDBpath
    db_interface_path =  os.path.join(os.path.dirname(SSDBpath),'DataInterfaces')
    #if os.path.isdir(db_interface_path):
    os.path.sys.path.append(db_interface_path)
    try:
      from Python import utils
      from Python import assp
    except ImportError:
      raise Exception(\
        'Script requires utils and assp from the SSDB python interface,\n' + \
        'but import failed.\n' + \
        'Maybe you are neither in the SSDB DataInterfaces folder\n' + \
        'nor in its Python subfolder nor have added the folder to your PYTHONPATH?' )
except Exception as e:
  print('Script requires utils and assp from the SSDB python interface,\n' + \
        'but import failed with following error message:\n' + \
        '  %s%\n' %str(e))

try:
  import typhon.arts.xml as tax
except ImportError:
  raise Exception(\
    'This script requires the typhon package. Retry after installing.')


####################
# Data extraction. (Likely) no need to modify.
####################

# Init database
utils.ssdb_init( SSDBpath )

# Import data, with some cropping in size and freq
S,M = assp.assp_import_ssdb( habID, orient, allow_nodata=False,
                             size_range=[minD, maxD], size_type=sizeparam,
                             freq_range=[fmin, fmax],
                             temp_range=[tmin, tmax] )

# Don't think, the following is needed. Hence, skip here.
# Some processing could be necessary. For example, to ensure that data are
# ordered in Dmax: 
#
#[dmax,ind] = unique( [ M.diameter_max ] )
#
#S = S(ind)
#M = M(ind)

# Convert S and M to ARTS internal format (assuming this is scat species 1)
#for i = 1 : length(S)
#  scat_data{1}{i} = S(i)
#  scat_meta{1}{i} = M(i)
#end

# Create files
assert(sdExt=='xml' or sdExt=='xml.gz' or sdExt=='xml.bin'), \
  'File extension needs to be either "xml", "xml.bin", or "xml.gz".'
extpart = sdExt.rpartition('.')
if extpart[-1]=='bin':
  sdExt = extpart[0]
  fmt = 'binary'
else:
  fmt = 'ascii'

filename = '%s/%s_%s' %(outpath,sdBase,sdSize)
if (sdSize=='dmax'):
  for i in range(len(S)):
    psize = M[i].diameter_max*1e6
    tax.save(S[i],'%s%04.0fum.%s' %(filename,psize,sdExt),format=fmt)
    tax.save(M[i],'%s%04.0fum.meta.%s' %(filename,psize,sdExt), format=fmt)
elif (sdSize=='dveq'):
  for i in range(len(S)):
    psize = M[i].diameter_volume_equ*1e6
    tax.save(S[i],'%s%04.0fum.%s' %(filename,psize,sdExt),format=fmt)
    tax.save(M[i],'%s%04.0fum.meta.%s' %(filename,psize,sdExt), format=fmt)
elif (sdSize=='mass'):
  for i in range(len(S)):
    psize = M[i].mass
    tax.save(S[i],'%s%.2ekg.%s' %(filename,psize,sdExt),format=fmt)
    tax.save(M[i],'%s%.2ekg.meta.%s' %(filename,psize,sdExt), format=fmt)
else:
  raise Exception(\
    "Size description parameter '%s' is unknown. Only 'dmax','dveq','mass' allowed." %sdSize)

