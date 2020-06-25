function [P]=shot_gathers_freq2P(XX)
% [P]=shot_gathers_freq2P(XX)
[nt,ng,ns]=size(XX);

for ii=1:nt/2
   for ss=1:ns
      P(:,ss,ii)=XX(ii,:,ss);
   end
end

