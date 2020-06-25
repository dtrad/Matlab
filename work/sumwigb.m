   function sumwigb(filename,filenameoff,nt,nh,NPLOT,mtitle,dt)
   if (nargin<5) NPLOT=111;end
   if (nargin<6) mtitle='';end
   if (nargin<7) dt=0.004;end	
   subplot(NPLOT)	
   fid = fopen([filename],'rb');
   [filename 'recbs']
   x=fread(fid,inf,'float32');
   x=reshape(x,nt,nh);
   		
   fclose(fid)
   
   fid = fopen([filenameoff],'rt') 
   h=fscanf(fid,'%d');

   fclose(fid)

   t=0:nt-1; t=t*dt;	
   wigb(x,1,h,t)
   xlabel('offset (m)');
   ylabel('time (sec)');
   title(mtitle);
   			
