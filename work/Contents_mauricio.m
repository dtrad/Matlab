% Seismic Processing Toolbox.
% Version 1.0   10-July-1997
%
% Input/Output SEGY data
%   header    - Reads the complete header of the k-th segy trace
%   extract   - Extract a header word from a SEGY file
%   readsegy  - Reads a segy file and trace headers  (Headers are optional)
%   writesegy - Writes a segy file and trace headers (        "           )
%   
% Processing 
%   clip     - Clip the data before plotting
%   ssort    - Sort data according to a given header word (i.e. sort in cdp's)
%   down_constant_v - Downward extrapolation in a constant velocity medium
%   
%   spike    - Wavelet in the center of a grid to test migration codes
%   ricker   - Generates a Ricker wavelet 
%
% Plotting.
%   sgray    - set a gray colormap with b&w clustering
%   simage   - Plot seismic data as an image 
%   wigb     - Wigle plot of seismic data
%
%
%
% NOTE: The following structure are used for SEG I/O.
%   To see them use load xx.mat. Never touch these files
%   they are important!!!.
%   The header follows the SEGY format of Seismic Unix (CMS).
%
%     segy.mat  - A structure saved as a matlab file that
%                 contains the header names and precission.
%    count.mat  - A structure that contains the position in bytes
%                 within the size of each header word.
%
%
%  M.D.Sacchi, July 1997, Dept. of Physics, UofA.
%
%  sacchi@phys.ualberta.ca
%
%



