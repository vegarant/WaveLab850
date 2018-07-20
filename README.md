# WAVELAB v 850

## About the project
This is a fork from the Stanford WaveLab library. The main purpose of this fork
is to fix some of the known bugs, resolve name conflicts, modernize the way the
WaveLab path is added and extend the support for orthogonal boundary wavelets
to a higher number of vanishing moments. The boundary wavelets does also use 
higher accuracy wavelet coefficients than what one usually finds in other 
implementations and they are computed for the minimum phase Daubechies wavelets. 
It does also support boundary wavelets for symmlets. 


## Installation

1. Download this repository.
2. Compile the MEX-files: Open Matlab and move to the WaveLab repository. Then 
make the function call `wl_install_mex` to compile all mex files.
3. Add WaveLab to your Matlab path: To do this in the current session, move to
the WaveLab repository and make a call to `wl_add_path.m`. To add WaveLab to
your Matlab path permanently do the following:  
    * Locate you Matlab `startup.m` file. You can do this by typing `userpath`
      in Matlab. If the userpath directory does not contain a `startup.m` file,
      create a file with this name i that folder.
    * Copy the content of the `wl_add_path.m` into your `startup.m` file. 
    * In the `startup.m` file modify the `wavelab_root` variable so that it
      points to the WaveLab root directory. To find the rigth value for this 
      variable you can move back to the WaveLab root directory and type
      `pwd`.
4. A few scripts are using the variable `WLVERBOSE` to determine how much
   information to print to the screen. By default verbose output is turned off.
   To turn this on type `global WLVERBOSE; WLVERBOSE = 'Yes';` in Matlab. To
   turn on verbosity permanently, add this in your `startup.m` file. 

