load vel_line500.txt
x=vel_line500

nmax=max(size(x));
i=1;
fid = fopen('tnmopar','w');

fprintf(fid,'tnmo=%d,',x(1,1));
for i=2:1:nmax-1
  if (x(i,1)<x(i+1,1))
    fprintf(fid,'%d,',x(i,1));
  else
    fprintf(fid,'%d',x(i,1));
    fprintf(fid,'\ntnmo=');    
  end
end


fclose(fid);


i=1;
fid = fopen('vnmopar','w');

fprintf(fid,'vnmo=%d,',x(1,2));
for i=2:1:nmax-1
  if (x(i,1)<x(i+1,1))
    fprintf(fid,'%d,',x(i,2));
  else
    fprintf(fid,'%d',x(i,2));
    fprintf(fid,'\nvnmo=');    
  end
end


fclose(fid);

