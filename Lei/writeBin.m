function n=writeBin(fileName,variable)
% write binary file for reading with other tools 
% last modified: 1/31/2018

fid=fopen(fileName,'wb');
n=fwrite(fid,variable,'float32')
return

end