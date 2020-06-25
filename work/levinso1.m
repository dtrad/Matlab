function a=levinso1(R,N,delay);
%LEVINSON  Levinson-Durbin Recursion.
%	A = LEVINSON(R,N) solves a symmetric toeplitz system of equations using
%	the Levinson-Durbin recursion.  R is a vector of autocorrelation
%	coefficients, starting with lag 0 as the first element. N is the
%	order of the recursion; A will be a length N+1 row, with A(1) = 1.
%
%	If you do not specify N, LEVINSON uses N = LENGTH(R)-1.
%
%	The equations solved are of the form:
%	    [  R(1)    R(2)  ...  R(N)  ] [  A(2)  ]  = [  -R(2)  ]
%	    [  R(2)    R(1)  ... R(N-1) ] [  A(3)  ]  = [  -R(3)  ]
%	    [  .       .           .    ] [   .    ]  = [    .    ]
%	    [ R(N-1)  R(N-2) ...  R(2)  ] [  A(N)  ]  = [  -R(N)  ]
%	    [  R(N)   R(N-1) ...  R(1)  ] [ A(N+1) ]  = [ -R(N+1) ]
%	If N is not large, LEVINSON will use the \ function to solve this
%	system, which is faster than the Levinson-Durbin recursion because
%	of its higher overhead.
%
%	See also LPC.

%	R is the auto correlation vector, R(1) = E(h(t)h*(t)), 
%		R(2) = E(h(t+1)h*(t)),...

%	Author(s): T. Krauss, 3-18-93
%	Copyright (c) 1984-94 by The MathWorks, Inc.
%	$Revision: 1.7 $  $Date: 1994/01/25 17:59:24 $

%	Reference(s):
% 	  [1] Lennart Ljung, "System Identification: Theory for the User",
%	      pp. 278-280

    error(nargchk(1,3,nargin))
    if nargin < 2, N = length(R)-1; end
    if length(R)<(N+1), error('Correlation vector too short.'), end

    if (N>195),     % use higher overhead, but asymptotically faster algorithm
        a = -R(2+delay-1)/R(1);
        V = R(1) - R(2+delay-1)^2/R(1);
        for n = 1:N-1,
            alfa = [1 a.']*R(n+2+delay-1:-1:2+delay-1);
            rho = -alfa / V;
            V = V + rho*alfa;
            a = [ a + rho*flipud(a); rho ];
        end
        a = [1; a].';
    else   % use the good old \ command
        b=-R(2+delay-1:N+1+delay-1);
        a = [ 1;  toeplitz(R(1:N))\(b(:)) ]';
    end

