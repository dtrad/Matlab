function writesudata(name,d)
%
% Usage writesudata(name,d)
%

fid=fopen(name,'w')

fwrite(fid,d(:),'float32');

fclose(fid);

return