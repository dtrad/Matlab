function [ A ] = readbinary( input,NZ,NX )
%Read binary file
fp=fopen(input,'r');
A=fread(fp,[NZ,NX],'float32');
fclose(fp);


end

