# CDFApprox
Gaussian CDF Approximation Functions <br />

The file "CDFNormalAproxPackNoC_0.1.5.tar.gz" is a simple package that contains the functions used to approximate the Gaussian CDF. This package contains multiple functions that allow these approximation to be done sequentially or in parallel.<br/>
The file "serial.R" is an example of approximating the Gaussian CDF sequentially. A 100x100 grid is simulated with an exponential covariance matrix.<br/>
The file "parallel.R" is an example of approximating the Gaussian CDF using parallel computing. This example requires the installation of the packages Rmpi and doMPI.
