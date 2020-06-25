function b=ndft(a,x,N,ft_opt,tflag)
%
% Computes the nonequispaced dft and its adjoint.
%
% b=ndft(a,x,N,options,tflag)
%
% a       Fourier coefficients or samples
% x       nodes \subset [-1/2,1/2)
% N       polynomial degree
% options .method   'direct'       computes the columns of the nonequispaced
%                                  Fourier matrix by direct calls of exp
%                   'horner'       avoids direct calls
%                   'matrix'       one matrix valued call of exp
%
% tflag   'notransp' evaluate trigonometric polynomial
%         'transp'   adjoint transform
%
% Author: Stefan Kunis

options=ft_opt.method;

M=length(x);
freq=(-N/2):(N/2-1);

if strcmp(tflag,'notransp')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nonequispaced dft
%
%         N/2-1  
% f(j) =  sum    f_hat(k+N/2+1)*exp(-2*pi*i*k*x(j)), 1 <= j <= M.
%        k=-N/2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  f_hat=a;
  switch options
    case 'direct'
      f=zeros(M,1);
      for k=1:N
        f=f+f_hat(k).*exp(-2*pi*i*x*freq(k));
      end;
    case 'horner'
      f=zeros(M,1);
      exp_x=exp(2*pi*i*x);
      for k=1:N
        f=(f+f_hat(k)).*exp_x;
      end;
      f=f.*exp(-N*pi*i*x);
    case 'matrix'
      f=exp(-2*pi*i*x*freq)*f_hat;
    otherwise
      disp('unknown option')
  end;
  b=f;

elseif strcmp(tflag,'transp')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjoint nonequispaced dft
%
%                    M  
% f_hat(k+N/2+1) =  sum    f(j)*exp(2*pi*i*k*x(j)), -N/2 <= k <= N/2-1.
%                   j=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  f=a;
  switch options
    case 'direct'
      f_hat=zeros(N,1);
      for k=1:N
        f_hat(k)=exp(-2*pi*i*x*freq(k))' *f;
        % exp(-2*pi*i*freq(k)*x)' *f differs by 10^-11!?
      end;
    case 'horner'
      f_hat=zeros(N,1);
      f=f.*exp(-N*pi*i*x);
      exp_x=exp(2*pi*i*x);
      for k=1:N
        f_hat(k)=sum(f);
        f=f.*exp_x;
      end;
    case 'matrix'
      f_hat=(exp(-2*pi*i*x*freq)') *f;
    otherwise
      disp('unknown option')
  end;
  b=f_hat;

end;
