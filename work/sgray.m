function sgray(alpha)
%SGRAY    SGTAY(alpha) set up a colormap for gray tones. 
%         alpha is the amount of clustering of white and black.
%         alpha=0.01 set negatives values to white and positive
%         to black. alpha=10 does a smooth transition.
%         This can be used before simage to enhance reflections
%         as an alternative to clipping.
%
%
%
%    M.D.Sacchi, July 1997, Dept. of Physics, UofA.
%        
%    sacchi@phys.ualberta.ca
%



i0=32;
i = 1:64;
t = (atan((i-i0)/alpha))';
s = t(64);
t = (t - min(t))*1./(max(t) -min(t));


m(1:64,2) = 1-t;
m(:,1) = 1-t;
m(:,3) =1-t;
colormap(m);
