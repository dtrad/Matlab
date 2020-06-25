	function [imp]=isitmp(x)

%  isitmp does what its name implies
%  x(nx) is the input
%  pac(nx) are the partial ac. coeficients
%  xr(nx) is the reconstituted trace  (x and xr should agree for MP)
%  ************************************************************
%	real*4 x(nx),pac(nx+1)
	nx=length(x);
	imp=0;
   pac=levfor(nx-2,x);

	for i=1:nx-2
	if(abs(pac(i+1))>=1)
		display('1+R(z) has |c|>or=1 at index '),i,
		imp=imp+1;
	end
   end

