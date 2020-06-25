   function sumwigbradar(filename,filenameoff,nt,nh,FLAG,NPLOT,mtitle,dt)
%  Plot seismic data. With SU use sustrip
%  function sumwigb(filename,filenameoff,nt,nh,FLAG,NPLOT,mtitle,dt)
%  
%  FLAG=0 ==> XT data (offset-t)
%  FLAG=1 ==> Radon data (q parameter-t)
%  FLAG=2 ==> Number of Trace - t
%  Example
%  sustrip < sudata15 > sudata15s outpar=/home/dtrad/radon/param
%  sugethw key=offset output=geom < sudata16 > sudata16off 
%  to produce offset file.
%  For Radon
%  sugethw key=f2 output=binary < sudata15rad > sudata15par
%  
%  Then plot in subplot(224) with sumwigb (better quality than supsiwgb)
%  x(512,64) title '(c)'
%  sumwigb('sudata15s','sudata15off',512,64,0,221,'(a)');
%  and for Radon
%  sumwigb('sudata15rads','sudata15par',512,100,1,222,'(b)');
%  Finally use print -dps name.ps to get the final psfile
%  Daniel Trad UBC- March 1999. 
 
   path='/home/dtrad/radar/';
   if (nargin<5) FLAG=0;end
   if (nargin<6) NPLOT=111;end
   if (nargin<7) mtitle='';end
   if (nargin<8) dt=0.004;end	
   subplot(NPLOT)
   filename=[path filename];
   
   if isempty(filenameoff) FLAG=2;end;
   filenameoff=[path filenameoff]; 	
   fid = fopen([filename],'rb');
   x=fread(fid,inf,'float32');
   x=reshape(x,nt,nh);   		
   fclose(fid)

   if (FLAG==0)
   fid = fopen([filenameoff],'rt'); 
   h=fscanf(fid,'%f');   
   fclose(fid);
   
   elseif (FLAG==1)
   fid = fopen([filenameoff],'rb') 
   h=fread(fid,inf,'float32');   
   fclose(fid)
   elseif (FLAG==2)
   h=1:nh;
   end

   if (length(h)~=nh) display('h and nh do not coincide');end  
   
   t=0:nt-1; t=t*dt;t=t-t(512);	
   wigb(x,1,h,t)
   xlabel('offset (m)');
   ylabel('depth (m)');
   title(mtitle);
   			












