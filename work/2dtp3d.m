fid=fopen('../data/ken2dto3d/section1','r')
d=fread(fid);
fclose(fid);
d=reshape(d,220,559);
wigb(d);