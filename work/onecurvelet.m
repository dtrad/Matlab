% In this script we remove groundroll using curvelets. We also remove
% incoherent noise.

%
clear;

calc     =    0;
Ground   =    0;
Denoise  =    1;
printing =    0;

% load data

FileName = 'datagroll'

[x]  =  readsudata([pwd,'/',FileName,'.bin'],256,256);
x    =  x/(norm(x(:))/256);

dlevs     =   [1,1,2,2];
pfilt     =   'db8';
dfilt     =   'pkva';

% Compute the threshold

if calc,
  nvar      =   pdfb_nest(256, pfilt, dfilt, dlevs);  
  save CurveGrollNoise_nvar nvar
else
    load CurveGrollNoise_nvar
end

% Make the mask

rlevs           =   cell(length(dlevs),1);
rlevs{1}{1}     =   1;
rlevs{1}{2}     =   1;
rlevs{2}{1}     =   1;
rlevs{2}{2}     =   1;
rlevs{3}{1}     =   1;
rlevs{3}{2}     =   0; %=0 
rlevs{4}{3}     =   0;
rlevs{4}{4}     =   0;
rlevs{5}{1}     =   1;
rlevs{5}{2}     =   1;
rlevs{5}{3}     =   0;
rlevs{5}{4}     =   0;

if Ground,

  cy        =    pdfbdec(x,pfilt,dfilt,dlevs);

  cx        =    ThresBand(cy,dlevs,rlevs);
  
  xn        =    pdfbrec(cx,pfilt,dfilt);

  figure(1)
  imagesc(x);
  color(seiscolor)
  title('Original data')

  figure(2)
  h1 = ImageCurvelet(cy,dlevs,1);
  axes(h1(1));
  title('Curvelet coefficients')
  
  figure(3)
  h2 = ImageCurvelet(cx,dlevs,1);
  axes(h2(1))
  title('Curvelet coefficients groundroll and noise')
  
  figure(4)
  imagesc(xn);
  colormap(seiscolor);
  title('Predicted Groundroll');
  
  figure(5)
  imagesc(x-xn)
  colormap(seiscolor)
  title('Data after groundroll removal')
  
  if printing
    for ifig = 1 : 5
      figure(ifig)
        print('-depsc2',['CurveGroll' num2str(ifig) '.eps'])
    end
  end
  writesudata('/home/dtrad/work/datagroll.noise.bin',xn);
end

% Do the denoising

if Denoise,
  
  sigma     =   0.1;
  z         =   x + 0*sigma*randn(size(x));
% forward CT
  cz1       =    pdfbdec(z,pfilt,dfilt,dlevs); 
% eliminate bands that contain ground roll
  cz        =    ThresBand(cz1,dlevs,rlevs);
% convert both to vectors  
  [vcz,s]   =    pdfb2vec(cz);
  [vcz1,s]  =    pdfb2vec(cz1);

% difference is the signal only   
  vcz       =    vcz1-vcz;
  
  fac         =    1;
  pdfb_thr    =    fac * 3 * sigma * nvar;

  fssize      =    prod(s(end, :));	% finest scale size

  pdfb_thr(end-fssize+1:end) = 4 / 3 * pdfb_thr(end-fssize+1:end);

  sorh        =    'h';                   % hard thresholding
% threshold signal only
  tvcz        =    wthresh(vcz, sorh, pdfb_thr);
  tcz         =    vec2pdfb(tvcz, s);
  zes         =    pdfbrec(tcz,pfilt,dfilt);

% difference with original coeff is the removed part only
  vcn         =    vcz1 - tvcz;
  tcn         =    vec2pdfb(vcn, s);
  zn          =    pdfbrec(tcn,pfilt,dfilt);
  
  
  figure(1)
  imagesc(z);
  color(seiscolor)
  title('Original data')

  figure(2)
  h1 = ImageCurveletmod(cz1,dlevs,1);
  axes(h1(1));
  title_handle = title('Curvelet coefficients');
  set(title_handle,'Visible','on')

  figure(6)
  h1 = ImageCurvelet(cz,dlevs,1);
  axes(h1(1));
  title_handle = title('Curvelet coefficients after thresband');
  set(title_handle,'Visible','on')
 
  figure(7)
  h1 = ImageCurvelet(tcz,dlevs,1);
  axes(h1(1));
  title_handle = title('Curvelet coefficients after thresband and thresholding');
  set(title_handle,'Visible','on')
 
  
  
  figure(3)
  h2 = ImageCurvelet(tcn,dlevs,1);
  axes(h2(1))
  title('Curvelet coefficients coherent groundroll and noise')
  
  figure(4)
  imagesc(zn);
  colormap(seiscolor);
  title('Predicted Groundroll and noise');
  
  figure(5)
  imagesc(zes)
  colormap(seiscolor)
  title('Data after groundroll and noise removal')
  
  if printing
    for ifig = 1 : 5
      figure(ifig)
        print('-depsc2',['CurveGrollNoise' num2str(ifig) '.eps'])
    end
  end

  writesudata('/home/dtrad/work/datagroll.RNden.bin',zes);
end

if (1)
  figure(2)
  print -dps /home/dtrad/ps/datagroll.curvetransf.ps
end

return
 




