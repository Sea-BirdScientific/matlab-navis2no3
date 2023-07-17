# matlab-navis2no3
Melt of SBS and MBARI code for processing nitrate from ISUS/SUNA reduced binary spectra in NAVIS isus message files

It’s largely based on the MATLAB code from MBARI's Josh Plant. 

See the processing options I setup in NAVIS2NO3_mbari.m (lines 91-120).  With the ‘FW’ options, you should be able to reproduce the nitrate values produced by SBS SUNA firmware C code.  

There were some small tweaks to Josh’s code in calc_FLOAT_NO3.m to 
* include both Sakamoto 2009 and Josh / Nehir’s recent TCSS code (lines 196-230)
* address the CTD / SUNA offset (see line 185)

I’ve included the correct SUNA cal files for five floats (1114, 1115, 60, 61, 62).
