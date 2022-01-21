# PARTICLES -- data structure for resampling analysis

- License: GPL v3
- Author: Andreas Juettner, andreas.juettner@cern.ch, juettner@soton.ac.uk

- Requires: [numpy](https://numpy.org) and [h5py](https://www.h5py.org)

- Initiatlise: `dat = particle()`

- Assign samples: `dat.samples = <inputdata>`
  where input data is a $N_{\rm boot} \times N_{\rm observable}$ array

- Methods:
    - Output: 
        - `dat.bs_val()`: sample mean
        - `dat.bs_dval()`: sample standard deviation
        - `dat.covariance()`: samaple covarianace

    - Safe/Read data to h5:
        - `dat.safe(filename)`
        - `dat.read(filename)`
  
    - Other:
        - `dat.info` is a dictionary that can be used to store additional information like fit quality, fit range etc. 
     
