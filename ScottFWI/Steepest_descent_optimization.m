function [ model , ft , eval ] = Steepest_descent_optimization( model , eval_grad , reg_func )

Wolfe1 = 1e-4;
Wolfe2 = 0.1;

[ ~ , g ] = eval_grad(model,0);

[~,rg]=reg_func(model);

g=g+rg;

descent_d = -g;

eval_grad = @(model)eval_grad(model,1);

%Starting step length
alpha0 = max(model)/(100*max(g));
% save check
[alpha,ft,eval] = Linesearch_Nocedal_Full( model , descent_d , Wolfe1 , Wolfe2 , eval_grad , alpha0 );

model = model + alpha*descent_d;

end

