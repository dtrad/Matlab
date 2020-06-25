% in this script I test some ideas with respect to the use of beamlets to
% find the stratigraphy.

calc = 0;                   % do the beamlet transform
threshold = 0;              % do the thresholding
decor1  = 1;                % decorate attributes
decor2  = 1;                % decorate data

% load data
load testdata

data   = testdata(1:128,1:128);
atf    = testatf(1:128,1:128);
loc    = find(atf);
locatf = zeros(size(atf));
locatf(loc) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the beamlet transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if calc==1
  n      = 128;
  img    = locatf;
  img2   = data;
  [E,L]   = SlowBeamTrans_m(img);
  [E2,L2] = SlowBeamTrans_m(img2);

  save testshot
else
  load  testshot
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do thresholding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if threshold==1
   max_img = max(img(:));
   min_img = min(img(:));

   thr = 0;
   cnt = 2;
   
   Zs = zeros(size(L)); Zs(:,:,cnt+1)=1;
   EL = E.*Zs;
   
   edgelets = find(EL>=thr);
   c(:,:,1) = (img - min_img)./(max_img - min_img);
   c(:,:,2) = (img - min_img)./(max_img - min_img);
   c(:,:,3) = (img - min_img)./(max_img - min_img);
   image(c); 
   axis image; ax = axis; 
   hold on; PlotBeamlets4(EL,edgelets,n);
   axis(ax); axis off; drawnow
   title(sprintf('scale=%d',cnt));
    save testshot
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the decorating
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if decor1==1
  lambda = 1;
  mEC = E.*(1e-10 + sqrt(L)) - lambda;
  [btree,vtree,stree] = BeamletRDP(mEC);
  imgp=img;
end

if decor2==1
  lambda = 2;
  mEC = E2.*(1e-10 + sqrt(L2)) - lambda;
  [btree,vtree,stree] = BeamletRDP(mEC);
  imgp = img2;
end

% do the plotting

fprintf('\n');
disp('There are two figures.');
fprintf('\n');

figure(1); 
AutoImage(-imgp); ax = axis; title('Picasso');

figure(2); 
c(:,:,1) = ones(n);
c(:,:,2) = ones(n);
c(:,:,3) = ones(n);
image(c); 
axis image; ax = axis; hold on; 
colormap('default')
nedges=PlotBeamletRDP2(vtree,btree,stree,'g',ax,log2(n),0,'w',-lambda+1e-3);
title('Recovered by Beamlet-driven RDP'); hold off; 
xlabel(sprintf('Number of beamlets=%d',nedges));

figure(3)
%AutoImage(-imgp); 
imagesc(-imgp);
colormap('gray')
axis image; ax = axis; hold on; 
nedges=PlotBeamletRDP2(vtree,btree,stree,'g',ax,log2(n),0,'w',-lambda+1e-3);
hold off;
