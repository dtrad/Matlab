	function [rxx]=acpeo(ncf,nrxx,g)
%****************************************************************
%	GERA (NRXX) COEFICIENTES DA FUNCAO DE AUTOCORRELACAO (RXX), DE
%	MODO QUE OS COEFICIENTES DO OPERADOR DE ERRO-PREDICAO (G) BURG
%	SEJAM SOLUCAO DAS EQUACOES NORMAIS ASSOCIADA A MATRIZ TOEPLITZ
%	DE AUTOCORRELACAO, QUANDO RESOLVIDA MEDIANTE RECURSAO DE LEVINSON.
%		NCF= N@ DE COEF. DO FILTRO AO ESTAGIO FINAL DA RECURSAO
%		NRXX = N@ DE COEFICIENTES DE RXX
%		RXX= VETOR P/ AUTOCORRELACAO
%		G  = VETOR P/ OPERADORES PREDITIVOS DE ERRO
%	OBS: O VETOR G ARMAZENA TODOS OS OPERADORES PREDITIVOS
%	     NA ORDEM CRESCENTE EM QUE FORAM GERADOS PELA RECURSAO
%	     G = {G11 ; G21,G22 ; G31,G32,G33 ; ...........}
%						MILTON J. PORSANI
%						PPPG-UFBA 24/02/85
%*****************************************************************
	rxx(1)=1;
	p = rxx(1);
	l = 0;
	for j=1:ncf
		sxx = 0;
		for k=1:j-1
			lk = l-k+1;
     		sxx = sxx + g(lk)*rxx(k+1);
   	end   
		l = l+j;
		rxx(j+1) = -(g(l)*p + sxx);
   	p = p*(1-g(l)^2);
   end
	in = l-ncf;
	for j=ncf+2:nrxx
		rxx(j) = 0;
		for k=1:ncf
			k1 = in+k;
			jk = j-k;
      	rxx(j) = rxx(j) -g(k1)*rxx(jk);
   	end
   end
	