"""
  This program is free software: you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
  You should have received a copy of the GNU General Public License along with
  this program. If not, see <http://www.gnu.org/licenses/>.


  Author: Andreas Juettner (Andreas.Juttner@cern.ch, juettner@soton.ac.uk)

"""

import numpy as np
import h5py
def depth(d):
     if isinstance(d,dict):
      # https://stackoverflow.com/questions/23499017/know-the-depth-of-a-dictionary
      return 1 + (max(map(depth, d.values())) if d else 0)
#     elif isinstance(d,list):
#      return 1 + (max(map(depth, d))) if d else 0 
     else:
      return 0
def read_hdf5_file(file,d):
    for key,val in file.items():
        if type(val) == h5py._hl.dataset.Dataset:
            d[key] = np.array(val,d[key])
#             print(key,np.array(val))

        else:
            d[key] = read_hdf5_file(val)
    return d

def info_to_dict(h5in,out):
    for key,item in h5in.items():
     if isinstance(item,h5py.Dataset):
       #print(key,item[()],np.array(item))
       if h5in[key]:
        try:
         out[key] = item[()]
        except:
         
         print('cannot read ',key,item)
     elif isinstance(item,h5py.Group):
      tmp = {}
      out[key] = info_to_dict(item,{})
     else:
      print('info_to_dict: noideia')
      exit()
    return out

def dict_to_info(group,dictionary):
    d = depth(dictionary)
    if d == 0:
     print('dict_to_info: not a dictionary')
     exit()
    if d == 1:
     for key,value in dictionary.items():
        group[key] = value 
     return group
    else:
     for key, value in dictionary.items():
        if isinstance(value, dict):
            # Recursively decompose nested dictionaries
            subgroup = group.create_group(key)
            subgroup = dict_to_info(subgroup,value)
        else:
            # Handle non-dictionary values
            try:
              if isinstance(value,bytes):
               group[key] = value.decode()
              else:
               group[key] = value
            except:
               continue
    return group

def fmean(bs,bare_data):
    return np.mean(bare_data[bs,:],0)
try:
    import multiprocessing
    MULTIPROCESSING = True
except ImportError:
    print('Module multiprocessing not found -- proceed with out it')
    print('Consider installing the module for better performance')
    MULTIPROCESSING = False
    pass
def par_kernel(chunk):
    f           = chunk[1]
    args        = chunk[2]
    res         = np.apply_along_axis(f,1,chunk[0],*args)
    return res
# apply_along_axis for multprocessing -- not supported, under development
def my_apply_along_axis(f,samples,*args,Nproc=6,parallel=True):
    if (MULTIPROCESSING == False) or (parallel==False):
     # use single-core version
     return np.apply_along_axis(f,1,samples,*args)
    else:
     # parallelise over Nproc cores
     # Chunks for the mapping (only a few chunks):
     chunks     = [(sub_arr,f,args) for sub_arr in np.array_split(samples, Nproc)]
     pool       = multiprocessing.get_context("fork").Pool(Nproc)
     individual_results = pool.map(par_kernel, chunks)
     # Freeing the workers:
     pool.close()
     pool.join()
     return np.concatenate(individual_results)


class particle():
  """
   PARTICLES v0
  """
  def __new__(cls, *args, **kwargs):
    return super().__new__(cls)
  def __init__(self,seed=1,Nboot=1000,bare_data=[],parallel=False,Nproc=1):
   """
    initiatlisation
   """
   self.cov_type= 'sample_mean'
   self.central	= 0.
   self.bs_mean	= 0.
   self.bs_mean	= 0.
   self.info	= {}
   self.samples	= []
   self.bs_indices = []
   self.N_bs   	= Nboot
   self.bs_seed = seed
   self.bare_data= bare_data
   self.parallel = parallel
   self.Nproc	 = Nproc
   self.NEWINFO  = True # improved version for storing info
   # if class initalised with list or array, create samples
   if isinstance(bare_data, list):
    if bare_data:
     self.bare_data=np.array(self.bare_data)
     self.make_samples()
   elif isinstance(bare_data, np.ndarray):
     self.make_samples()
  def bs_val(self):
   """
    compute bootstrap average
   """
   return np.mean(np.array(self.samples),0)
  def bs_dval(self):
   """
    compute bootstrap error
   """
   return np.sqrt(np.diag(self.covariance()))
  def printnice(self):
   try: 
    import jn
    return jn.disperr(self.bs_val(),self.bs_dval())
   except:
    return('')
  def covariance(self):
    """ 
     compute covariance matrix
    """
    if self.cov_type == 'sample_mean':
     v   = self.bs_val()
    elif self.cov_type == 'ensemble_mean':
     v   = self.central
    else:
     print('particles.covariance(): unknown cov_type: %s'%self.cov_type)

    d   = np.array(self.samples)-v
    self.cov  = np.dot(d.T,d)/len(self.samples)
    return self.cov
  def get_all_data(self):
   return [self.bs_val(),self.bs_dval(),self.covariance(),self.samples]
  def make_samples(self):
   """
     Draw bootstrap samples. Assumes data to be NxT matrix,
     where N is the number of MC samples and T the number of observables
     (e.g. T for time in correlation functions).
   """
   if self.bs_indices == []: # Skip drawing BS indices if already allocated
    K = self.N_bs
    N = self.bare_data.shape[0]
    np.random.seed(self.bs_seed)
    # draw BS indices
    self.bs_indices = np.random.randint(0,high=N,size=(N,K)).T
   # I checked -- apply_along_axis takes about the same time as list comprehension 
   # (commented out below)
   #fmean2 = lambda bs: np.mean(self.bare_data[bs,:],0)
   #self.samples = np.apply_along_axis(fmean2,1,self.bs_indices)
   #print(self.samples)
   self.samples = my_apply_along_axis(fmean,self.bs_indices,self.bare_data,Nproc=self.Nproc,parallel=self.parallel)
   #print(self.samples)
   #exit()
#   tmp = np.array([np.mean(self.bare_data[bs,:],0) for bs in self.bs_indices])
   # clean up results after new data loaded
   self.cov = np.array([])
   self.bs_mean = []
   self.ens_mean = []
   self.bs_error = []
   self.computed_cov = 0
   self.computed_mean = 0
  def safe(self,f):
    """
     store content to disk
    """
    # store the samples to h5
    f=h5py.File(f,'w')
    f.create_dataset('samples',data=self.samples)
    f.create_dataset('central',data=self.central)
    f.create_dataset('bs_mean',data=self.bs_mean)
    if self.bare_data!=[]:
#     print "storing bare_data"
     f.create_dataset('bare_data',data=self.bare_data)
    #
    if self.NEWINFO:
     #
     # Create a group to store the dictionary
     group = f.create_group('info')
     group = dict_to_info(group,self.info)
     #
    else: # old version
     for item in self.info.keys():
      if isinstance(self.info[item],list):
       if any(isinstance(i, list) for i in self.info[item]): # check if nested list
        for ii,l in enumerate(self.info[item]):
         f.create_dataset('info_'+item+'_'+str(ii),data=l)
       else:
         #print(item,self.info[item])
         f.create_dataset('info_'+item+'_'+str(0),data=self.info[item])
      else:
       f.create_dataset('info_'+item,data=self.info[item])
     
    f.close
  def read(self,f):
    """
     read content from disk
    """
    f=h5py.File(f,'r')
    self.samples = np.array(f.get('samples'))
    self.central = np.array(f.get('central'))
    self.bs_mean = np.array(f.get('bs_mean'))
    if 'bare_data' in f:
     self.bare_data = np.array(f.get('bare_data'))
    info={}
    for item in f.keys():
     if 'info' in item:
      dat = np.array(f.get(item))
      info[item.replace('info_','')]=dat
    self.info = info
    f.close
  def read2(self,f):
    """
     read content from disk
    """
    f=h5py.File(f,'r')
    self.samples = np.array(f.get('samples'))
    self.central = np.array(f.get('central'))
    self.bs_mean = np.array(f.get('bs_mean'))
    if 'bare_data' in f:
     self.bare_data = np.array(f.get('bare_data'))
    info={}
    info = info_to_dict(f.get('info'),info)
    for item in f.keys():
     if 'info' in item:
      dat = np.array(f.get(item))
      info[item.replace('info_','')]=dat
    self.info = info
    f.close

