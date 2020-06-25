function [tau]=stress(lambda,mu,e)
% function [tau]=stress(lambda,mu,e)
% Function to compute the stress 
% Seismology course eosc 354
% Daniel Trad - UBC 

tau(1,1)=lambda*(e(1,1)+e(2,2))+2*mu*e(1,1);
tau(2,2)=lambda*(e(1,1)+e(2,2))+2*mu*e(2,2);
tau(1,2)=2*mu*e(1,2);
tau(2,1)=2*mu*e(2,1);
