function [vrms,t,vint]=vz2vrms(vel,z,dt,tmax)
%  VZ2VT: Compute V(x,t) and Vrms(x,t) from V(x,z) and a ZOS exploding
%         reflector model of V(x,z)
%
%  This function takes in a depth velocity model and a ZOS exploding
%  reflector model based on the velocity model and computes a velocity
%  matrix in time and an rms velocity matrix in time
%
%  This function is useful for converting depth models to time velocity models
%  to test the capabilities of the CREWES migration routines on ZOS images
%
%  [vrms,t]=vz2vt(vel,z,dt,tmax)
%
%  vel..........is the input velocity model in depth. Each row is a
%               constant depth.
%  z........depth coordinate for vel, length(z) = size(vel,1) 
%  dt.... desired time sample rate
%  tmax ... maximum two-way time desired
%
%
%  vrms....is the output RMS velocity matrix in time
%  t....... output two-way time coordinate
%  vint .... interval velocity versus time
%
%
%  Zoron Rodriguez, November 2006
%
% NOTE: It is illegal for you to use this software for a purpose other
% than non-profit education or research UNLESS you are employed by a CREWES
% Project sponsor. By using this software, you are agreeing to the terms
% detailed in this software's Matlab source file.

% BEGIN TERMS OF USE LICENSE
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its author (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) Use of this SOFTWARE by any for-profit commercial organization is
%    expressly forbidden unless said organization is a CREWES Project
%    Sponsor.
%
% 2) A CREWES Project sponsor may use this SOFTWARE under the terms of the
%    CREWES Project Sponsorship agreement.
%
% 3) A student or employee of a non-profit educational institution may
%    use this SOFTWARE subject to the following terms and conditions:
%    - this SOFTWARE is for teaching or research purposes only.
%    - this SOFTWARE may be distributed to other students or researchers
%      provided that these license terms are included.
%    - reselling the SOFTWARE, or including it or any portion of it, in any
%      software that will be resold is expressly forbidden.
%    - transfering the SOFTWARE in any form to a commercial firm or any
%      other for-profit organization is expressly forbidden.
%
% END TERMS OF USE LICENSE



if(nargin<4)
    error('not enough input arguments')
end
if(length(z)~=size(vel,1))
    error('vel and z sizes are not compatible')
end
nx=size(vel,2);
t=(0:dt:tmax);
vrms=zeros(length(t),nx);
vint=vrms;
yell=0;
for k=1:nx
   tv=vint2t(vel(:,k),z);
   %if(tv(end)<tmax)
   %    tv(end)=tmax;
   %end
   vint(:,k)=interpextrap(tv',vel(:,k),t/2);
   vrms(:,k)=vint2vrms(vel(:,k),tv,t/2);
   if(k>nx/2 && yell==0)
       disp('Hang on, almost done!')
       yell=1;
   end
end




