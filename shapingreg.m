function [ alpha, beta, p, m, r, res_2, norm_m ] = shapingreg( niter, r, m, p, H, L, lambda )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    res_2 = zeros(length(niter),1); % least-squares residual/objective function
    norm_m = zeros(length(niter),1);
for n = 1:niter(end)
    gm2 = transpose(L)*r-lambda*m; 
    gp2 = transpose(H)*gm2+lambda*p;
    gm2 = H*gp2;
    gr2 = L*gm2;
    p = transpose(gp2)*gp2;
    if n == 1
        sp2 = gp2;
        sm2 = gm2;
        sr2 = gr2;
        beta = 0;
    else
        beta = p/p_prep2;
        sp2 = gp2+beta*sp2;
        sm2 = gm2+beta*sm2;
        sr2 = gr2+beta*sr2;
    end
    p_prep2 = p;
    
    alpha = p/(transpose(sr2)*sr2+lambda*(transpose(sp2)*sp2-transpose(sm2)*sm2));
    p = p-alpha*sp2;
    m = m-alpha*sm2;
    r = r-alpha*sr2;
    
    norm_m(n) = transpose(m)*m;
    % Step 7: Update the least-squares norm of residuals and relative error
    res_2(n) = transpose(r)*r;
end
end

