	function [g1,g]=levbak(n,gpp)


%   levbak uses the levinson backward scheme to reconstitute the
%  polynomial using the reflection coefficients
%
%*******************************************************************
%	PROPOSITO:	DADA A SERIE DE COEFICIENTES DE REFLEXAO (GPP)
%			E' FEITO USO DA RELACAO DE LEVINSON
%					_
%			Gp = Gp-1 + gpp Gp-1	PARA GERAR TODOS OS
%
%			OPERADORES PREDITIVOS GP(s)
%	DESCRICAO DOS PARAMETROS:
%		ENTRADA:	N - N@ DE COEFICIENTES DE REFLEXAO
%				GPP - COEFICIENTES DE REFLEXAO
%				     {g11,g22,g33,....gnn}
%		SAIDA:		G1  - ARRANJO QUE ARMAZENA TODOS OS
%				     OPERADORES PREDITIVOS
%				     {g11;g21,g22;...;gn1,...,gnn}
%				G   - VETOR DE TRABALHO QUE AO FINAL
%				      ARMAZENARA' O ULTIMO GP
%
%	AUTOR:	MILTON J. PORSANI  PPPG-UFBA  1/10/85
%
%************************************************************************

	in=0;
	for j=1:n
	g(j)=gpp(j);
	for i=1:fix(j/2)
		ji = j-i;
		ga1 = g(ji);
		g(ji) = ga1	+ gpp(j)*g(i);
		if(ji~=i) g(i) = g(i)+ gpp(j)*ga1;end
	end
	for i=1:j
		ii = in+i;
      g1(ii) = g(i);
   end   
	in = in + j;
	end
	g(2:n)=g(1:n-1);g(1)=1;	