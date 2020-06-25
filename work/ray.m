 function[v]=ray(xl,yl,xr,yr,dx,dy,kx,ky)
% [v]=ray(xl,yl,xr,yr,dx,dy,kx,ky)
% gives a line in the tomography corresponds to transmitter in xl,yl, 
% and a receiver in xr,yr. The media is divided into kx times ky cells
% of size dx,dy
% Input: xl - x position of transmitter
%        yl - y position of transmitter
%        xr - x position of receiver
%        yr - y position of receiver
%        dx - size of a cell (x direction)
%        dy - size of a cell (y direction)
%        kx - How many cells in the x direction
%        ky - How many cells in the y direction
%
% Output: v - A vector size (kx*ky) of ray lengths.
% Remark The transmitters should always be on the left and the
%        receivers on the right as follows:
%
%        |                |
%       s*                |
%        |                |
%        |                >r
%        |                |
%


% Check that transmitters are on the left and receivers on the 
% right an if not change between them

%====== Preparations =====================================================
 if xr<xl, a=xr;xr=xl;xl=a;a=yr;yr=yl;yl=a; end;

% ==== First possibility transmitter and receiver have the same y ===========

  if yl==yr
       xpoint=[0:dx:kx*dx];
       ypoint=ones(1,kx+1)*yl;
  else

%=====================Step 1 ============================================= 
% Calculating the slope of each ray and its intersection

  a=(yr-yl)/(xr-xl);
  b=yl-a*xl;

%======================Step 2 ============================================ 
% Finding and Sorting all the intersections
     
  x=[ceil(xl/dx)*dx:dx:floor(xr/dx)*dx];
  if a>0,
        y=[ceil(yl/dy)*dy:dy:floor(yr/dy)*dy];
  elseif a<0,
        y=[floor(yl/dy)*dy:-dy:ceil(yr/dy)*dy];
  else, end;
  
  xy=(y-b)/a;
  xpoint=[x,xy];
  xpoint=sort(xpoint);

  
  ypoint=a*xpoint+b;
end 
%===================== Step 3 =============================================
% Assigning cell number for each intersection
 
  xmidpnt=xpoint+0.5*[diff(xpoint),0];
  ymidpnt=ypoint+0.5*[diff(ypoint),0];

  cellnum=floor(ymidpnt/dy)*kx+floor(xmidpnt/dx)+1;

%==================== Step 4 ==============================================  
% Putting all this in a vector
 
  v=zeros(1,kx*ky);
  
  for i=1:length(cellnum)-1, 
       
       if (cellnum(i)<=kx*ky & cellnum(i)>0)
             v(1,cellnum(i))=((ypoint(i+1)-ypoint(i))^2+...
               (xpoint(i+1)-xpoint(i))^2)^0.5;
       else
       end;
   end;


%================= Done Calculation ======================================   
% If there are any trouble take the remark of the next part  
% to see each ray.
             figure(1)
              pcolor(0:dx:dx*kx,0:dy:dy*ky,ones(ky+1,kx+1));
              colormap(1-gray);
              hold on 
              plot(xpoint,ypoint,'b') 
              plot(xmidpnt,ymidpnt,'o'); 
              axis([0,kx*dx,0,ky*dy]);
              axis('ij')
              hold off 

              
 
