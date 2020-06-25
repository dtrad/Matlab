function [sigmapfin]=findbeta(sigmap,JD,targetphi)
% Function to find the sigmap that produces the target misfit
% Input: 
% 		sigmap: values computed previously
% 		JD: data misfit computed previously.
% Output:
% 		sigmapfin: value of sigmap that produces JD close to the target misfit
% Geop 521B- Daniel Trad- 09/04/98
MM=max(size(sigmap));
% Interpolation of sigmap-JD using 
XID=logspace(log10(sigmap(1)),log10(sigmap(MM)),100);
YID=interp1(sigmap,JD,XID);

%=========================================================
% Find the beta for the misfit required
format short e;
display('Searching sigmap')
I=find(YID<=targetphi);
sigmapfin=min(XID(I));
%=========================================================
