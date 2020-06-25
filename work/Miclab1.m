% Initialize stress tensor.
% Landers.
%t11=3294;
%t12=-45644;
%t22=81351;

% 1000 years.
t11=9.196e6;
t12=3.31e5;
t22=1.192e6;



tij=[t11,t12;t12,t22];

% Loop over azimuth. 
az=[0:10:170];
nz=length(az);
for iz=1:nz
  azm=az(iz)*pi/180;

% Unit vectors normal and parallel to plane.
  n=[cos(azm);-sin(azm)];
  s=[sin(azm);cos(azm)];

% Traction.
  ti=tij*n;

% Normal and shear stress.
  tn(iz)=ti'*n;
  ts(iz)=ti'*s;

% CFF
  cff(iz)=abs(ts(iz))+0.2*tn(iz);
  
end

%
figure

subplot(311);plot(az,tn,'o');title('Michael');
subplot(312);plot(az,ts,'o');,
subplot(313);plot(az,cff,'o');

