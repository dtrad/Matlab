fid=fopen('../data/ken2dto3d/section1','r')
T=fread(fid,inf,'float');
fclose(fid);
T=reshape(T,559,220);
%wigb(d2);
%figure;
%simage(T);

fid=fopen('../data/ken2dto3d/section2','r')
P=fread(fid,inf,'float');
fclose(fid);
P=reshape(P,559,220);
%wigb(d2);
%figure;
%simage(P);

