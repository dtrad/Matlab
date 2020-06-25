function c=levfor(n,g)
%   c=levfor(n,g)
%   levfor determines the reflection coefficients associated
%   with a polynomial where g(0)=1 using the forward
%   Levinson recursion.
%
%*******************************************************************
%	PROPOSITO:	DADA O OPE E' FEITO USO DA RELACAO
%					_
%			Gp = [Gp-1 - gpp Gp-1]/[1-gpp]
%
%		PARA RECUPERAR TODOS OS COEFICIENTES DE REFLEXAO gpp
%
%	DESCRICAO DOS PARAMETROS:
%		ENTRADA:	N - COMPRIMENTO DO OPE,
%				    EX: (1,G21,G22), N=3
%				G - OPE
%
%		SAIDA:		c  - COEFICIENTES DE REFLEXAO
%				     {1,g11,g22,g33,....gnn}
%	OBS: O INPUT E' DESTRUIDO
%
%	AUTOR:	MILTON J. PORSANI  PPPG-UFBA  1/10/85
%
%************************************************************************

	g=g(1:n+1);
	c=g;
	
	for j=n-1:-1:2
		j1=j+1;
		eg=1-g(j1)*g(j1);
      for i=1:fix(j/2)
         i1=i+1;
			ji = j1-i;
         ga1  = g(ji);
         g(ji)= (ga1-g(j1)*g(i1))/eg;
			if(ji~=i1) g(i1)= (g(i1) - g(j1)*ga1)/eg;end
		end
	end

	temp=c;
	c=g;
	g=temp;





