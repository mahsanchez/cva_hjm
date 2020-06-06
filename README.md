# cva / interest rate swap / hjm interest rate model

Reduced memory usage and batched Monte Carlo simulation to be suitable for CUDA acceleration [02/06/2020]
See Method MonteCarloSimulation::generatePath in header simulation.h
       
Use boost function boost::math::erf_inv to generate gaussian variates based on a merssene twister uniform random generator [02/06/2020]

Use Intel MKL to generate gaussian variates based on a merssene twister uniform random generator [07/06/2020]

TODO - Modify YieldCurve calculation DF(t, T) and ForwardRates F(t;t, T). Data must be accessed using a linearized array [03/06/2020]

TODO - Test MonteCarlo Convergence for the CVA IRS value [03/06/2020]
