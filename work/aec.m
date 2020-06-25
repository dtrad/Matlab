function trout=aec(trin,t,op_length,trip)

%  trout=aec(trin,t,op_length,trip)
%  trout=aec(trin,t,op_length)
%
% AEC performs and automatic amplitude adjustment:
% Method:
%   1) Compute Hilbert envelope of trin
%   2) Convolve envelope with triangular smoother of half-length
%       op_length
%   3) Divide trin by smoothed envelope
%
% trin= input trace
% t= time coordinate vector for trin
% op_length= half-length of triangular smoother in seconds
% trip= front end time before which the smoothed envelope is
%        set to a constant ******** default= op_length/10 ******
% trout= output trace
%
% by G.F. Margrave, May 1991
%
% NOTE: It is illegal for you to use this software for a purpose other
% than non-profit education or research UNLESS you are employed by a CREWES
% Project sponsor. By using this software, you are agreeing to the terms
% detailed in this software's Matlab source file.
 
% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by 
% its author (identified above) and the CREWES Project.  The CREWES 
% project may be contacted via email at:  crewes@geo.ucalgary.ca
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
   
% set defaults
 if nargin<=3
   trip=op_length/10.;
 end
% double the operator length
 op2=op_length*2;
% form new trace padded to a power of 2
% trinew=padpow2(trin,0);
% compute the envelope
 env=abs(hilbm(trinew));
 env=env(1:length(trin));
% guard against accidental transpose
           if size(env)~=size(trin), env=env.';end
% compute the smoothed envelope
 nop=round(op2/(t(2)-t(1)))+1;
 envsm=conv(env,triang(nop));
% grab the central length(trin) samples
 envsm=envsm(round(nop/2):length(trin)+round(nop/2)-1);
% stabilize the envelope at the ends
 ntrip=round(trip/(t(2)-t(1)))+1;
 envsm=[envsm(ntrip)*ones(1,ntrip) envsm(ntrip+1:length(envsm))];
 envsm=[envsm(1:length(envsm)-ntrip) envsm(length(envsm)-ntrip)*ones(1,ntrip)];
% correct the trace
 trout=trin./envsm;
% balnce the output to have the same mean power as input
trout=balans(trout,trin); 








