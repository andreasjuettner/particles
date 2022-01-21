import numpy as py
import h5py
class particles():
  """
   PARTICLES v0
  """
  def __init__(self):
   """
    initiatlisation
   """
   self.central	= 0.
   self.bs_mean	= 0.
   self.bs_mean	= 0.
   self.info	= {}
   self.samples	= []
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
  def covariance(self):
    """ 
     compute covariance matrix
    """
    v   = self.bs_val()
    d   = np.array(self.samples)-v
    self.cov  = np.dot(d.T,d)/len(self.samples)
    return self.cov
  def safe(self,f):
    """
     store content to disk
    """
    # store the samples to h5
    f=h5py.File(f,'w')
    f.create_dataset('samples',data=self.samples)
    f.create_dataset('central',data=self.central)
    f.create_dataset('bs_mean',data=self.bs_mean)
    for item in self.info.keys():
     if isinstance(self.info[item],list):
      if any(isinstance(i, list) for i in self.info[item]): # check if nested list
       for ii,l in enumerate(self.info[item]):
        f.create_dataset('info_'+item+'_'+str(ii),data=l[ii])
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
    info={}
    for item in f.keys():
     if 'info' in item:
      dat = np.array(f.get(item))
      info[item.replace('info_','')]=dat
    self.info = info
    f.close

