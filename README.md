# cva / interest rate swap / hjm interest rate model

Reduced memory usage and batched Monte Carlo simulation to be suitable for CUDA acceleration [02/06/2020]
See Method MonteCarloSimulation::generatePath in header simulation.h
       
Use boost function boost::math::erf_inv to generate gaussian variates based on a merssene twister uniform random generator [02/06/2020]

Use Intel MKL to generate gaussian variates based on a merssene twister uniform random generator [07/06/2020]

Fixed The Exposure Calculation profile for a given Interest Rate Swap Mark to Market [09/06/2020]

Added new dissertation notes with the Reference used during the study [09/06/2020]

Finished MonteCarlo Convergence for IRS Expected Exposure with a 0.01 accuracy and 4750 simulations on each simulation point.[14/06/2020]

TODO - Estimate the stddev and confidence interval for simulations

TODO - Apply a variance reduction method on HJM model for risk factor evolution to speed up the MC convergence

TODO - Transpose the Exposure matrix and check for memory bandwidth usage improvement

TODO - Modify YieldCurve calculation DF(t, T) and ForwardRates F(t;t, T). Data must be accessed using a linearized array 


