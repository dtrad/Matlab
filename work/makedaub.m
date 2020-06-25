% toon0111 -- Wavelet Families
%
%   Wavelet analysis begins by choosing a specific family of wavelets
%   to work with.
%
%   The family is specified by a father and a mother wavelet, and
%   these generate a basis by translation and dilation.
% 
%   Here we illustrate four specific Mother wavelets
%
%       Haar          -- the first wavelet; a square-wave wavelet
%
%       Daubechies D4 -- the first continuous, compactly supported
%                        orthonormal wavelet family
%
%       Coiflet C3    -- orthonormal wavelets system where both father and
%                        mother have special vanishing moments properties
%
%       Symmlet S8    -- nearly-symmetric orthogonal wavelet of
%                        compact support with 8 vanishing moments.
% 
clear	
%
	t = (1:4096)*2;
%

%wave1 = MakeWavelet(4,2,'Daubechies',4,'Mother',4096);
%wave2 = MakeWavelet(6,24,'Daubechies',4,'Mother',4096);
%wave3 = MakeWavelet(8,140,'Daubechies',4,'Mother',4096);
%wave4 = MakeWavelet(10,750,'Daubechies',4,'Mother',4096);
%wave5 = MakeWavelet(6,30,'Daubechies',4,'Mother',4096);
%wave6 = MakeWavelet(6,50,'Daubechies',4,'Mother',4096);

wave1 = MakeWavelet(4,2,'Daubechies',4,'Mother',4096);
wave2 = MakeWavelet(5,12,'Daubechies',4,'Mother',4096);
wave3 = MakeWavelet(6,38,'Daubechies',4,'Mother',4096);
wave4 = MakeWavelet(7,100,'Daubechies',4,'Mother',4096);

wave=wave1+wave2+wave3+wave4;
subplot(111);
plot(t(:)./1000,wave(:));% title(' D4 Wavelet ');
axis([0 7 min(wave(:)) max(wave(:)) ]);
%axis([2 7 -0.12 0.12 ]);
xlabel('');ylabel('Amplitude');
figure(gcf)    
%wave=wave4+wave5+wave6;

    
%   
% Part of WaveLab version .701
% Built Tuesday, January 30, 1996 8:25:59 PM
% This is Copyrighted Material
% For copying permissions see copying.m
% Comments? e-mail wavelab@playfair.stanford.edu
%   
    
   
% Auto name mapping to DOS conventions  Wednesday, January 31, 1996 5:36:19 PM
