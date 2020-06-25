function [AAR]=hessl2l1cost(x,AA,lambda)
Cmi=diag(1./(max(abs(x),1e-5)));
AAR=AA+lambda*Cmi;