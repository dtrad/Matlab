% add some paths
addpath('C:\Documents and Settings\kdemeers\My Documents\Matlab\third_party\S-transform');
addpath('C:\Documents and Settings\kdemeers\My Documents\Matlab\third_party\nersc_toolbox\signal');
addpath('C:\Documents and Settings\kdemeers\My Documents\Matlab\seismic_utl');

% load the data
if (~exist('seismic_x') | ~exist('seismic_y') | ~exist('seismic_z'))
    load ../data
end

lmovel=1000

% Apply gain correction and filter the data
filterfreq=[1 2 20 25];
gain=(1:1:nsamples)'.^1.5;
seismic_x_filt=cosfilt(seismic_x(:,:).*(gain*sampint*ones(1,ntraces)),...
    sampint,filterfreq);
seismic_y_filt=cosfilt(seismic_y(:,:).*(gain*sampint*ones(1,ntraces)),...
    sampint,filterfreq);
seismic_z_filt=cosfilt(seismic_z(:,:).*(gain*sampint*ones(1,ntraces)),...
    sampint,filterfreq);

% Resample the data
resamp_fact=floor(1/(sampint*2)/(2*max(filterfreq)));
resamp_sampint=sampint*resamp_fact;
resamp_nsamples=ceil(nsamples/resamp_fact);
fprintf(['The frequency filtered raw data will be subsampled by a'...
    ' factor of %i. \n'],resamp_fact)

stinput_x=zeros(round(resamp_nsamples*1.5),ntraces);
stinput_y=zeros(round(resamp_nsamples*1.5),ntraces);
stinput_z=zeros(round(resamp_nsamples*1.5),ntraces);

stinput_x(1:resamp_nsamples,:)=seismic_x_filt(1:resamp_fact:nsamples,:);
stinput_y(1:resamp_nsamples,:)=seismic_y_filt(1:resamp_fact:nsamples,:);
stinput_z(1:resamp_nsamples,:)=seismic_z_filt(1:resamp_fact:nsamples,:);

resamp_nsamples=round(resamp_nsamples*1.5);

% apply divergence correction
fprintf('Applying lmo correction to the data... ')
stinput_x=lmo(stinput_x,offset,lmovel,resamp_sampint*1000);
stinput_y=lmo(stinput_y,offset,lmovel,resamp_sampint*1000);
stinput_z=lmo(stinput_z,offset,lmovel,resamp_sampint*1000); 
fprintf('done! \n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of data preparation part
%
% Start of filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize arrays to hold the noise estimate
noisefreq_x=zeros(1+ceil(max(filterfreq)*(resamp_sampint*resamp_nsamples)),ntraces);
noisefreq_y=zeros(1+ceil(max(filterfreq)*(resamp_sampint*resamp_nsamples)),ntraces);
noisefreq_z=zeros(1+ceil(max(filterfreq)*(resamp_sampint*resamp_nsamples)),ntraces); 

for freq=floor(min(filterfreq)*(resamp_sampint*resamp_nsamples)):1:ceil(max(filterfreq)*(resamp_sampint*resamp_nsamples))
    % Initialize the arrays for holding the S-transform frequency slice
    % Initialize the arrays for holdingv other parameters
    st_x=zeros(1,resamp_nsamples,ntraces);
    st_y=zeros(1,resamp_nsamples,ntraces);
    st_z=zeros(1,resamp_nsamples,ntraces);
    ellipticity=zeros(1+ceil(max(filterfreq)*(resamp_sampint*resamp_nsamples)),ntraces); 
    dip=zeros(1+ceil(max(filterfreq)*(resamp_sampint*resamp_nsamples)),ntraces); 
    azimuth=zeros(1+ceil(max(filterfreq)*(resamp_sampint*resamp_nsamples)),ntraces); 
    
    frequency=(freq)/(resamp_nsamples*resamp_sampint);
    
    % Compute S-transform frequecny slice
    for station=1:ntraces
        [st_x(:,:,station),t,f]=st(stinput_x(:,station),freq,freq,1,1);
        [st_y(:,:,station),t,f]=st(stinput_y(:,station),freq,freq,1,1);
        [st_z(:,:,station),t,f]=st(stinput_z(:,station),freq,freq,1,1);
    end
    st_x=squeeze(st_x);
    st_y=squeeze(st_y);
    st_z=squeeze(st_z);
    
    
    for trace=1:ntraces
        for sample=1:2/3*resamp_nsamples-round(offset(trace)/lmovel/resamp_sampint)
            % Do polarizatoin analysis, 
            vec=[st_x(sample,trace);st_y(sample,trace);st_z(sample,trace)];
            vec=vec/norm(vec);
            
            % Compute ellipticity, azimuth and dip of ellipse
            A=sum(real(vec).^2+imag(vec).^2);
            B=sum(real(vec).^2-imag(vec).^2);
            C=-2*sum(real(vec).*imag(vec));
            a=sqrt(A+sqrt(B^2+C^2))/sqrt(2);
            b=sqrt(A-sqrt(B^2+C^2))/sqrt(2);
            
            % Ellipticity
            ellipticity(sample,trace)=b/(a);
            dip(sample,trace)=atan(sqrt((real(vec(3))*imag(vec(2))-real(vec(2))*imag(vec(3)))^2+...
                                        (real(vec(3))*imag(vec(1))-real(vec(1))*imag(vec(3)))^2)/...
                                  (real(vec(2))*imag(vec(1))-real(vec(1))*imag(vec(2))));
            if dip(sample,trace) < 0
                dip(sample,trace)=pi+dip(sample,trace);
            end
            azimuth(sample,trace)=atan( (real(vec(3))*imag(vec(2))-real(vec(2))*imag(vec(3)))/...
                                        (real(vec(3))*imag(vec(1))-real(vec(1))*imag(vec(3))));
            
           
            
            % obtain noise estimate
            
        end
    end
    % Compute masking function
    %             mask=(ellipticity<0.5).*(abs(azimuth)>pi/2).*(dip > pi/4 & dip < 3*pi/4);
    
%     mask=(ellipticity>0.5).*(abs(azimuth)<pi/4).*(dip > pi/4 & dip < 3*pi/4);
    mask=(abs(azimuth)<pi/4).*(dip > pi/4 & dip < 3*pi/4);
    noisefreq_x(freq+1,:)=sum(st_x(1:size(mask,1),:).*mask);
    noisefreq_y(freq+1,:)=sum(st_y(1:size(mask,1),:).*mask);
    noisefreq_z(freq+1,:)=sum(st_z(1:size(mask,1),:).*mask);

    
    
    
    
    fprintf('Frequency: %3.2f Hz\n',frequency)
    if rem(frequency,2) <= 1/(resamp_nsamples*resamp_sampint)
        figure('windowstyle','docked')
        subplot(1,4,1)
        temp=cosfilt(seismic_x(:,:).*(gain*sampint*ones(1,ntraces)),sampint,frequency*[0.8 0.9 1.1 1.2]);
        imagesc(temp(1:resamp_fact:nsamples,:))
        title('Radial')
        caxis(.01*[-1 1])
        ylim([0 resamp_nsamples/1.5])
        colormap('gray')
        colorbar('location','southoutside')
        
        subplot(1,4,2)
        temp=cosfilt(seismic_y(:,:).*(gain*sampint*ones(1,ntraces)),sampint,frequency*[0.8 0.9 1.1 1.2]);
        imagesc(temp(1:resamp_fact:nsamples,:))
        title('Transverse')
        caxis(.01*[-1 1])
        ylim([0 resamp_nsamples/1.5])
        colorbar('location','southoutside')
        
        subplot(1,4,3)
        temp=cosfilt(seismic_z(:,:).*(gain*sampint*ones(1,ntraces)),sampint,frequency*[0.8 0.9 1.1 1.2]);
        imagesc(temp(1:resamp_fact:nsamples,:))
        title('Vertical')
        caxis(.01*[-1 1])
        ylim([0 resamp_nsamples/1.5])
        colorbar('location','southoutside')
        
        subplot(1,4,4)
        imagesc(ilmo(10*log(abs(st_x).^2+abs(st_y).^2+abs(st_z).^2),offset,lmovel,resamp_sampint*1000))
        str=sprintf('3C DB: Frequency: %3.2f Hz\n',frequency);
        title(str)
        colorbar('location','southoutside')
        ylim([0 resamp_nsamples/1.5])
        caxis([-150 0])
        
        figure('windowstyle','docked')
        subplot(1,4,1)
        imagesc(ilmo(real(ellipticity),offset,lmovel,resamp_sampint*1000))
        str=sprintf('ELLIPTICITY: Frequency: %3.2f Hz\n',frequency);
        title(str)
        caxis([0 1])
        colorbar('location','southoutside')
        ylim([0 resamp_nsamples/1.5])
        colormap('jet')
        
        subplot(1,4,2)
        imagesc(ilmo(real(dip),offset,lmovel,resamp_sampint*1000))
        title('DIP')
        caxis(pi*[0 1])
        colorbar('location','southoutside')
        ylim([0 resamp_nsamples/1.5])
        
        subplot(1,4,3)
        imagesc(ilmo(abs(azimuth),offset,lmovel,resamp_sampint*1000))
        title('AZIMUTH')
        caxis(pi*[0 1])
        colorbar('location','southoutside')
        ylim([0 resamp_nsamples/1.5])
       
        subplot(1,4,4)
        imagesc(ilmo(mask,offset,lmovel,resamp_sampint*1000))
        title('mask')
        caxis([0 1])
        colorbar('location','southoutside')
        ylim([0 resamp_nsamples/1.5])
        
    end
end
% Radial component
temp=noisefreq_x;
temp(resamp_nsamples*resamp_fact-size(noisefreq_x,1)+2:resamp_nsamples*resamp_fact,:)=flipud(conj(noisefreq_x(2:size(noisefreq_x,1),:)));
temp=ilmo(real(ifft(temp*resamp_fact)),offset,lmovel,resamp_sampint/resamp_fact*1000);
temp=temp(1:nsamples,:);

figure('windowstyle','docked')
subplot(1,3,1)
imagesc(seismic_x_filt)
title('INPUT: Radial')
caxis(.05*[-1 1])
colormap('gray')

subplot(1,3,2)
imagesc(seismic_x_filt-temp)
title('OUTPUT: Radial')
caxis(.05*[-1 1])
colormap('gray')

subplot(1,3,3)
imagesc(temp)
title('DIFFERENCE: Radial')
caxis(.05*[-1 1])
colormap('gray')

% Transverse component
temp=noisefreq_y;
temp(resamp_nsamples*resamp_fact-size(noisefreq_y,1)+2:resamp_nsamples*resamp_fact,:)=flipud(conj(noisefreq_y(2:size(noisefreq_y,1),:)));
temp=ilmo(real(ifft(temp*resamp_fact)),offset,lmovel,resamp_sampint/resamp_fact*1000);
temp=temp(1:nsamples,:);

figure('windowstyle','docked')
subplot(1,3,1)
imagesc(seismic_y_filt)
title('INPUT: Transverse')
caxis(.05*[-1 1])
colormap('gray')

subplot(1,3,2)
imagesc(seismic_y_filt-temp)
title('OUTPUT: Transverse')
caxis(.05*[-1 1])
colormap('gray')

subplot(1,3,3)
imagesc(temp)
title('DIFFERENCE: Transverse')
caxis(.05*[-1 1])
colormap('gray')


% Vertical component
temp=noisefreq_z;
temp(resamp_nsamples*resamp_fact-size(noisefreq_z,1)+2:resamp_nsamples*resamp_fact,:)=flipud(conj(noisefreq_z(2:size(noisefreq_z,1),:)));
temp=ilmo(real(ifft(temp*resamp_fact)),offset,lmovel,resamp_sampint/resamp_fact*1000);
temp=temp(1:nsamples,:);

figure('windowstyle','docked')
subplot(1,3,1)
imagesc(seismic_z_filt)
title('INPUT: Vertical')
caxis(.05*[-1 1])
colormap('gray')

subplot(1,3,2)
imagesc(seismic_z_filt-temp)
title('OUTPUT: Vertical')
caxis(.05*[-1 1])
colormap('gray')

subplot(1,3,3)
imagesc(temp)
title('DIFFERENCE: Vertical')
caxis(.05*[-1 1])
colormap('gray')
