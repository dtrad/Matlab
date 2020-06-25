
load ../poland/spick1/fxstackvel.txt

cdp=fxstackvel(:,1);
tnmo=fxstackvel(:,2);
vnmo=fxstackvel(:,3);

fid = fopen('cdps.dat','w');
fprintf(fid,'cdp=');
%fprintf(fid,'%3.0f,',cdp(1));
ii=0;
for i=1:length(cdp)-1
  if (cdp(i)~=cdp(i+1)) 
    fprintf(fid,'%3.0f,',cdp(i));ii=ii+1;
  end
end
fprintf(fid,'%3.0f\n',cdp(i));ii=ii+1;

fclose(fid)


nmax=max(size(cdp))
i=1;
fid = fopen('tnmopar','w');

fprintf(fid,'tnmo=%d,',tnmo(1));
for i=2:1:nmax-1
  i,tnmo(i)
  if (tnmo(i)<tnmo(i+1))
    fprintf(fid,'%5.3f,',tnmo(i));
  else
    fprintf(fid,'%5.3f',tnmo(i));
    fprintf(fid,'\ntnmo=');    
  end
end


fclose(fid);


i=1;
fid = fopen('vnmopar','w');

fprintf(fid,'vnmo=%4.0f,',vnmo(1));
for i=2:1:nmax-1
  %i,vnmo(i)
  if (tnmo(i)<tnmo(i+1))
    fprintf(fid,'%4.0f,',vnmo(i));
  else
    fprintf(fid,'%4.0f',vnmo(i));
    fprintf(fid,'\nvnmo=');    
  end
end
fclose(fid);

% Merge the cdps, tnmo and vnmo files into a parfile
fid1 = fopen('tnmopar','r');
fid2 = fopen('vnmopar','r');
fid3 = fopen('parfile','w');
fid4 = fopen('cdps.dat','r');

S0=fgets(fid4);
fprintf(fid3,S0);
for i=1:ii;
 S1 = fgets(fid1);
 S2 = fgets(fid2);
 fprintf(fid3,'%s%s',S1,S2);
end 
fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);