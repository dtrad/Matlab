clear all; close all; clc;
%setpath_Windows;
load linereg
%load data
dtrue = data;
dt = 0.004;
[nt,nx] = size(dtrue);

d = zeros(size(dtrue));
% define sampling operator: 50% missing traces
tmp = rand(1,nx);
T = zeros(size(tmp));

for n = 1 : nx
    if tmp(n) > 0.5;
        T(n) = 1;
        d(:,n) = dtrue(:,n);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for pocs_2d
f_low  = 0.1;
f_high = 80;
option = 3; %data-driven threshold schedule
perc_i = 99;
perc_f = 0.1;
N = 100;
a = 0.9;
tol = 1e-14;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
[dout] = pocs_2d(d,T,dt,f_low,f_high, option, perc_i, perc_f, N, a, tol);

figure(1);
subplot(141)
wigb(dtrue,1,[1:nx],([1:nt]-1)*dt,1);
xlabel('trace number'); ylabel('time (s)'); title('dtrue');
subplot(142)
wigb(d,1,[1:nx],([1:nt]-1)*dt,1);
xlabel('trace number'); ylabel('time (s)'); title('d');
subplot(143)
wigb(dout,1,[1:nx],([1:nt]-1)*dt,1);
xlabel('trace number'); ylabel('time (s)'); title('dout');
subplot(144)
wigb(dout-dtrue,1,[1:nx],([1:nt]-1)*dt,1);
xlabel('trace number'); ylabel('time (s)'); title('dout-dtrue');



