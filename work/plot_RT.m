function [data,v,d,h0,h,t,q,tau]=plot_RT(file)
% Program to plot all Tau_p programs hr_ss, hr_vs, hr_parab.
% All the options imply 1=true, 0 = false
% Daniel Trad - UBC - Canada. 6/11/98
nmo=0; %0, 1 nmo.  
ss=0; % Slant Stack
mute=0; % Muting
residuals=0; % plots of residuals
% Original data set
pathfile='c:\daniel\radon\';

filename=[pathfile file];
nt=512;
h0=load ([filename '.off']);   %offset
nh=length(h0);

 
data=read_bin_double([filename '.bin']);
data=reshape(data,nt,nh);

[v,av,bv]=read_hb([pathfile 'vel_gather.out']);
if mute==1 
   [v2,av,bv]=read_hb([pathfile 'vel_gather2.out']);
   v2=reshape(v2,av(2),av(1));
end
[d,ad,bd]=read_hb([pathfile 'rec_data.out']);
%[dlu,ad,bd]=read_hb([pathfile 'rec_dataLU.out']);

if nmo==1
   [dnmo,adnmo,bdnmo]=read_hb([pathfile 'datanmo.out']);
   dnmo=reshape(dnmo,adnmo(2),adnmo(1));
	t=0:ad(2)-1;
	t=t*bd(2)+bd(4);	
	figure,wigb(dnmo,1,h0,t);title('NMO corrected data'),
	xlabel('offset (m)');ylabel('time (s)');
end

h=0:ad(1)-1;
h=h*bd(1)+bd(3);
t=0:ad(2)-1;
t=t*bd(2)-bd(4);

q=0:av(1)-1;
q=q*bv(1)+bv(3);
tau=0:av(2)-1;
tau=tau*bv(2)-bv(4);


figure,
data=normalize(data);
scaled=normalize_energy(data,d);
%subplot(221),wigb(data,1,h0,t)
title('Data');xlabel('offset');ylabel('t (sec)')

v=reshape(v,av(2),av(1));
%subplot(222),wigb(v,1,q,tau)
title('Velocity Gather');xlabel('q');ylabel('\tau')

max_data=max(max(abs(data)));

d=reshape(d,ad(2),ad(1));
%dlu=reshape(dlu,ad(2),ad(1));
%d=normalize(d);
%dlu=normalize(dlu);
%error=max(max(d-dlu))
max_d=max(max(abs(d)));

d=normalize(d);
if (ss==1)
	subplot(223),wigb(d,scaled,h0,t)
elseif ss~=1
%   subplot(223),wigb(d,scaled,h,t)
end
title('Data recovered');xlabel('offset');ylabel('t')
if (mute==1)
   subplot(224),wigb(v2,1,q,tau);
   title('Muted Radon domain');xlabel('p');ylabel('\tau')
end

if (residuals==1)
r=data-d;
subplot(224),wigb(r,max_data,h,t)
title('Residuals');xlabel('offset');ylabel('t')
end
scaled=1;
%figure,
wigb(d,scaled,h,t)
title('Data recovered');xlabel('offset');ylabel('t')

v=reshape(v,av(2),av(1));
figure,
wigb(v,1,q,tau)
title('Velocity Gather');xlabel('q');ylabel('\tau')

figure,
wigb(data,1,h0,t)
title('Data');xlabel('offset');ylabel('t (sec)')
