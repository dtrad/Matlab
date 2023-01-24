function [ model , ft , eval ] = Gauss_Newton_optimization( model , eval_grad , eval_Hess , reg_func )

Wolfe1 = 1e-4;
Wolfe2 = 0.9;

[ H , g , ~ ] = eval_Hess(model);

[~,rg,rH]=reg_func(model);

H=H+(rH*mean(abs(diag(H))));
g = g + rg*mean(abs(diag(H)));

descent_d = -H\g;

%Starting step length
alpha0 = 1;
[alpha,ft,eval] = Linesearch_Nocedal_Full( model , descent_d , Wolfe1 , Wolfe2 , eval_grad , alpha0 );

model = model + alpha*descent_d;

end

