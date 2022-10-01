The main functions are provided by Nonlinear Analysis Core, Center for Human MovementVariability, University of Nebraska at Omaha. Here is the [link to their github](https://github.com/Nonlinear-Analysis-Core/NONANLibrary) but in some cases were modified by me (noted below).

This is a list of the included functions and the full name of the methods.

- **Ent_Samp** - Used to calculate the Sample Entropy of a time series.
- **Surr_findrho** - This should be used to find thee optimal noise level used in creating a Pseudo Period surrogate time series.
- **Surr_PseudoPeriodic** - This is used to find a Pseudo Period surrogate time series using the result from Surr_findrho.
- **Surr_Theiler** - This can be used to create different surrogates by the methods published by Theiler. *I modified the original function to include the iaaft shuffling for Algorithm 2. The original function simply did not have Algorithm 2 implemented.*
- **iaaft** - Performs iterative amplitude-adjusted fourier transform. The functions is from the Multifractal toolbox by E. Ihlen available [here](https://www.ntnu.edu/inb/geri/software)


COPYRIGHT

Copyright 2021 Nonlinear Analysis Core, Center for Human Movement Variability, University of Nebraska at Omaha
