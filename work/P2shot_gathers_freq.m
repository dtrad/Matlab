function [XX]=P2shot_gathers(P)
% [XX]=P2shot_gathers(P)
[ng,ns,nt]=size(P);
for ii=1:nt
   for ss=1:ns
      XX(ii,:,ss)=P(:,ss,ii);
   end
end

