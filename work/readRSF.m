function vel=readRSF(fileName,NZ,NX)
% read .rsf file from madagscar to .mat
% last modified: 1/31/2018

fid=fopen(fileName,'r');
v=fread(fid,NZ*NX,'float32');
vel=reshape(v,[NZ,NX]);


end