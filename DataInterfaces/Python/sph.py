# -*- coding: utf-8 -*-

"""Collection of functions for Spherical Harmonics transformation of the SSDB data."""


try:
  import numpy as np
except ImportError:
  raise Exception( 'Module *numpy* required but not found.' )

try:
  import scipy.interpolate as spi
except ImportError:
  raise Exception( 'Module *scipy.interpolate* required but not found.' )

try:
  import shtns
except ImportError:
  raise Exception( 'Module *shtns* required but not found.\n' +
                    'shtns is the python module of the SHTns library for' +
                    ' Spherical Harmonics Transformations.' )



#==============================================================================
#  Spherical Harmonics 
#==============================================================================


class Spharmt(object):
    #2017-10-13 Manfred Brath
    """
    wrapper class for commonly used spectral transform operations in
    atmospheric models.  Provides an interface to shtns compatible
    with pyspharm (pyspharm.googlecode.com). Taken from the shallow_water-example
    of the shtns and modified by Manfred Brath
    """
    def __init__(self, nlons, nlats, ntrunc, rsphere,
                 gridtype='gaussian', alpha=0):
        """initialize
        nlons:  number of longitudes
        nlats:  number of latitudes"""
        self._shtns = shtns.sht(ntrunc, ntrunc, 1, 
                                shtns.sht_orthonormal+shtns.SHT_NO_CS_PHASE)     
        if gridtype == 'gaussian':
            self._shtns.set_grid(nlats, nlons,
                                 shtns.sht_quick_init|shtns.SHT_PHI_CONTIGUOUS,1.e-10)  
        elif gridtype == 'regular':
            self._shtns.set_grid(nlats, nlons,
                                 shtns.sht_reg_poles|shtns.SHT_PHI_CONTIGUOUS,1.e-10)            
        self.lats = np.arccos(self._shtns.cos_theta)
        self.lons = (2.*np.pi/nlons)*np.arange(nlons)
        self.nlons = nlons
        self.nlats = nlats
        self.ntrunc = ntrunc
        self.nlm = self._shtns.nlm
        self.degree = self._shtns.l
        self.lap = -self.degree*(self.degree+1.0).astype(np.complex128)
        self.invlap = np.zeros(self.lap.shape, self.lap.dtype)
        self.invlap[1:] = 1./self.lap[1:]
        self.rsphere = rsphere
        self.lap = self.lap/rsphere**2
        self.invlap = self.invlap*rsphere**2
    def grdtospec(self, data):
        """compute spectral coefficients from gridded data"""
        return self._shtns.analys(data)
    def spectogrd(self, dataspec):
        """compute gridded data from spectral coefficients"""
        return self._shtns.synth(dataspec)
    def getuv(self, vrtspec, divspec):
        """compute wind vector from spectral coeffs of vorticity and divergence"""
        return self._shtns.synth((self.invlap/self.rsphere)*vrtspec,
                                 (self.invlap/self.rsphere)*divspec)
    def getvrtdivspec(self, u, v):
        """compute spectral coeffs of vorticity and divergence from wind vector"""
        vrtspec, divspec = self._shtns.analys(u, v)
        return self.lap*self.rsphere*vrtspec, self.lap*self.rsphere*divspec
    def getgrad(self, divspec):
        """compute gradient vector from spectral coeffs"""
        vrtspec = np.zeros(divspec.shape, dtype=np.complex128)
        u,v = self._shtns.synth(vrtspec,divspec)
        return u/self.rsphere, v/self.rsphere
    def zrotate(self, dataspec, alpha):
        """rotate spectral coefficients by angle alpha"""
        RLM=self._shtns.Zrotate(dataspec,alpha)
        return RLM
    def yrotate(self, dataspec, beta):
        """rotate spectral coefficients by angle beta (along y-axis)"""
        RLM=self._shtns.Yrotate(dataspec,beta)
        return RLM
    def gauss_wts(self):
        """computes gaussian weights for theta integration"""
        return self._shtns.gauss_wts()
    def idx(self, l, m):
        """get the index of the spectral array for coefficient l,m"""
        return self._shtns.idx(l,m)    
    def SH_to_point(self, dataspec, cos_theta, phi):
        """compute spatial point from spectral coefficients"""
        return self._shtns.SH_to_point(dataspec,cos_theta,phi)


            
#==============================================================================
#   get truncation (lmax) from sph coefficients      
#==============================================================================
    
def get_lmax(sph_coeffs):
    #2017-10-13 Manfred Brath
    ''' Function to get the truncation limit or in other words the maximum
        of the highest l coefficient'''
    
    #get truncation or lmax
    lmax=int(1/2+np.sqrt(1/4+2*len(sph_coeffs))-2)
    
    
    return lmax


#==============================================================================
#   map spherical expansion to another spherical expansion    
#==============================================================================

def map_sphexp(sph_coeffs, sph_object):
    #2017-10-13 Manfred Brath
    '''Function to map from one spherical harmonics expansion to another 
    Important the expansion, which is mapped to another has to be
    of smaller or same lmax as the expansion, to which it is mapped.'''
    
    #get truncation or lmax
    ntrunc=int(1/2+np.sqrt(1/4+2*len(sph_coeffs))-2)

    #create local spherical object
#    sph_local=shtns.sht(ntrunc, ntrunc, 1, shtns.sht_orthonormal+shtns.SHT_NO_CS_PHASE)
#    local_m=sph_local.m
#    local_l=sph_local.l

    #build array with the l coefficients   
    local_l=np.zeros(len(sph_coeffs), dtype=int)    
    start=0
    end=0
    for i in range(0,ntrunc+1):
        
        lnew=np.arange(i,ntrunc+1,1)
        
        end=end+len(lnew)
        
        local_l[start:end]=lnew
        
        start=end
        
    

    
    #build array with the m coefficients   
    local_m=np.zeros(len(sph_coeffs),dtype=int)    
    start=0
    end=0
    for i in range(0,ntrunc+1):
        
        mnew=np.ones(ntrunc+1-i)*i
        
        end=end+len(mnew)
        
        local_m[start:end]=mnew
        
        start=end


    #allocate
    sph_coeffs_new=np.zeros(sph_object.nlm,dtype=complex)
    
    
    if sph_object.nlm<len(sph_coeffs):
        
        RuntimeError('target expansion has less elements than source expansion!')
    
    #do the mapping
    for i in range(0,len(local_l)):
        
        idx=sph_object.idx(int(local_l[i]),int(local_m[i]))
        
        sph_coeffs_new[idx]=sph_coeffs[i]
        
        
        
        
    return sph_coeffs_new        

        
#==============================================================================
#     
#==============================================================================

def reshape_sphexp_2_common_exp(max_ntrunc_pha, max_len_coeffs_pha, phase_mat):
    #2017-10-13 Manfred Brath
    '''Function to reshape an array of sph_expansions to a common expansion'''

#    #create spherical harmonics object
#    if max_ntrunc_pha % 2 ==0:
#        min_grid_size=max_ntrunc_pha+2
#    else:
#        min_grid_size=max_ntrunc_pha+1
    
    sph_object=shtns.sht(int(max_ntrunc_pha), int(max_ntrunc_pha), 1,
                         shtns.sht_orthonormal+shtns.SHT_NO_CS_PHASE)
    
    # reshape phase_arrays to common sph_expansion
    #allocate
    phase_mat_array=np.zeros((len(phase_mat),len(phase_mat[0]),16,
                              int(np.max(max_len_coeffs_pha))),dtype=complex)
    
    for i in range(0,len(phase_mat)):
    
        phase_mat_i=phase_mat[i]
        
        for j in range(0,len(phase_mat_i)):
            
            phase_mat_ij=phase_mat_i[j]
            
            for k in range(0,len(phase_mat_ij)):
            
                phase_mat_array[i,j,k,:]=map_sphexp(phase_mat_ij[k],sph_object) 
    
    return phase_mat_array    


#==============================================================================
# create phase matrix in spatial representation
#==============================================================================

def create_regular_representation(phase_mat_sph, max_ntrunc,
                                  grid_size_theta=90, gridtype='gaussian'):
    #2017-10-13 Manfred Brath
    """create phase matrix in spatial grid representation"""
    
    if gridtype=='gaussian':
    
        if max_ntrunc % 2 ==0:
            min_grid_size=max_ntrunc+2
        else:
            min_grid_size=max_ntrunc+1
            
        if min_grid_size<grid_size_theta:
            min_grid_size=grid_size_theta 
        elif grid_size_theta<min_grid_size:
            grid_size_theta=min_grid_size
                
        if min_grid_size % 2 == 1:
            min_grid_size=min_grid_size+1
        
        
        #create spherical harmonics object
        x = Spharmt(min_grid_size*2,min_grid_size,max_ntrunc,1,gridtype='gaussian')
        
    elif gridtype=='regular': 
        
        grid_size_theta_old=grid_size_theta
        
        if max_ntrunc % 2 ==0:
            min_grid_size=2*max_ntrunc+2
        else:
            min_grid_size=2*(max_ntrunc+1)
            
        if min_grid_size<grid_size_theta:
            min_grid_size=grid_size_theta 
        elif grid_size_theta<min_grid_size:
            grid_size_theta=min_grid_size
                
        if min_grid_size % 2 == 1:
            min_grid_size=min_grid_size+1
        
        
        
        x = Spharmt(min_grid_size*2,min_grid_size,max_ntrunc,1,gridtype='regular')
        
    
    #create spatial grid
    #phi_s,theta_s = np.meshgrid(x.lons, x.lats)
    phi_s=x.lons*180/np.pi
    theta_s=x.lats*180/np.pi
    
    qm=np.shape(phase_mat_sph)
    #ql=np.shape(phi_s)
    
    #allocate
    phase_matrix_reg=np.ones((qm[0],qm[1],qm[2],len(theta_s),len(phi_s)))
    
    
    #transform to spatial grid
    for i in range(0,qm[0]):#loop over the beta angles
            
        for j in range(0,qm[1]):#loop over the incidence angles
            
            for k in range(0,qm[2]):# loop over the matrix elements
                    
                phase_matrix_reg[i,j,k,:,:]=x.spectogrd(np.complex128(phase_mat_sph[i,j,k,:]))
                
                
    
        
    if gridtype=='regular' and (
        grid_size_theta_old<grid_size_theta or grid_size_theta_old % 2==1):
        
        #to have uneven grid size or grid sizes that are smaller than the 
        #minimum grid size, we need to interpolate.
            
        f_int=spi.interp1d(theta_s,phase_matrix_reg,axis=3,kind='linear',
                                assume_sorted=True)
    
        theta_s=np.linspace(0,180,num=grid_size_theta_old)
        phase_matrix_reg=f_int(theta_s)
        
        
    #interpolate to assure that azimuth grid always include 0 and 180 and if
    #the grid size are smaller than the minimum grid size.    
    f_int=spi.interp1d(phi_s,phase_matrix_reg,axis=4,kind='linear',
                            assume_sorted=True)

    phi_s=np.linspace(0,180.,num=grid_size_theta_old)                
    phase_matrix_reg=f_int(phi_s)
    
    #filter out negative values of the (1,1) or in python (0,0)
    dummy=phase_matrix_reg[0,...]    
    logic0=dummy<0
    dummy[logic0]=0
    phase_matrix_reg[0,...]=dummy
    
    return phase_matrix_reg,theta_s,phi_s

