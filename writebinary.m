function [N]=writebinary(output,A)
[n m]=size(A);
fid=fopen(output,'w');
N=fwrite(fid,A,'float32');
if (n*m ~= N) display('WARNING: no all data was written')
fclose(fid);
end

