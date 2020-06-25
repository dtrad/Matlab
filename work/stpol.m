% add some paths
addpath('C:\Documents and Settings\kdemeers\My Documents\Matlab\third_party\S-transform');
addpath('C:\Documents and Settings\kdemeers\My Documents\Matlab\third_party\nersc_toolbox\signal');
addpath('C:\Documents and Settings\kdemeers\My Documents\Matlab\seismic_utl');

% load the data
if (~exist('seismic_x') | ~exist('seismic_y') | ~exist('seismic_z'))
    load data
end

lmovel=1000

% filter and subsample the data
filterfreq=[1 2 20 25];
seismic_x_filt=cosfilt(seismic_x(:,:).*(  (1:1:nsamples).^1.5'*sampint*ones(1,ntraces)),sampint,filterfreq);
seismic_y_filt=cosfilt(seismic_y(:,:).*(  (1:1:nsamples).^1.5'*sampint*ones(1,ntraces)),sampint,filterfreq);
seismic_z_filt=cosfilt(seismic_z(:,:).*(  (1:1:nsamples).^1.5'*sampint*ones(1,ntraces)),sampint,filterfreq);
% seismic_x_filt=cosfilt(seismic_x(:,:),sampint,filterfreq);
% seismic_y_filt=cosfilt(seismic_y(:,:),sampint,filterfreq);
% seismic_z_filt=cosfilt(seismic_z(:,:),sampint,filterfreq);

resamp_fact=floor(1/(sampint*2)/(2*max(filterfreq)));
resamp_sampint=sampint*resamp_fact;
resamp_nsamples=ceil(nsamples/resamp_fact);
fprintf('The frequency filtered raw data will be subsampled by a factor of %i. \n',resamp_fact)

stinput_x=zeros(round(resamp_nsamples*1.5),ntraces);
stinput_y=zeros(round(resamp_nsamples*1.5),ntraces);
stinput_z=zeros(round(resamp_nsamples*1.5),ntraces);

stinput_x(1:resamp_nsamples,:)=seismic_x_filt(1:resamp_fact:nsamples,:);
stinput_y(1:resamp_nsamples,:)=seismic_y_filt(1:resamp_fact:nsamples,:);
stinput_z(1:resamp_nsamples,:)=seismic_z_filt(1:resamp_fact:nsamples,:);

resamp_nsamples=round(resamp_nsamples*1.5);

% apply divergence correction
fprintf('Applying divergence correction and lmo correction to the data. \n')
gain=((1:1:resamp_nsamples)'*resamp_sampint).^0;
for ii=1:ntraces
    stinput_x(1:resamp_nsamples,ii)=stinput_x(1:resamp_nsamples,ii).*gain;
    stinput_y(1:resamp_nsamples,ii)=stinput_y(1:resamp_nsamples,ii).*gain;
    stinput_z(1:resamp_nsamples,ii)=stinput_z(1:resamp_nsamples,ii).*gain;
end
% return
stinput_x=lmo(stinput_x,offset,lmovel,resamp_sampint*1000);
stinput_y=lmo(stinput_y,offset,lmovel,resamp_sampint*1000);
stinput_z=lmo(stinput_z,offset,lmovel,resamp_sampint*1000);
counter=0;
% compute the s-transforms one frequency slice at a time 
% (otherwise memory blows up)
noisefreq_x=zeros(1+ceil(max(filterfreq)*(resamp_sampint*resamp_nsamples)),ntraces);
noisefreq_y=zeros(1+ceil(max(filterfreq)*(resamp_sampint*resamp_nsamples)),ntraces);
noisefreq_z=zeros(1+ceil(max(filterfreq)*(resamp_sampint*resamp_nsamples)),ntraces); 
for freq=floor(min(filterfreq)*(resamp_sampint*resamp_nsamples)):1:ceil(max(filterfreq)*(resamp_sampint*resamp_nsamples))
    st_x=zeros(1,resamp_nsamples,ntraces);
    st_y=zeros(1,resamp_nsamples,ntraces);
    st_z=zeros(1,resamp_nsamples,ntraces);
    frequency=(freq)/(resamp_nsamples*resamp_sampint);
    ellipticity=zeros(1+ceil(max(filterfreq)*(resamp_sampint*resamp_nsamples)),ntraces); 
    for station=1:ntraces
        [st_x(:,:,station),t,f]=st(stinput_x(:,station),freq,freq,1,1);
        [st_y(:,:,station),t,f]=st(stinput_y(:,station),freq,freq,1,1);
        [st_z(:,:,station),t,f]=st(stinput_z(:,station),freq,freq,1,1);
    end
    st_x=squeeze(st_x);
    st_y=squeeze(st_y);
    st_z=squeeze(st_z);
    
    halvewidth=2;
    twohalvewidth=2*halvewidth;
    width=twohalvewidth+1;
    for ii=1:round(resamp_nsamples*2/3)
        for j=1+halvewidth:ntraces-halvewidth
            D=[st_x(ii,j-halvewidth:j+halvewidth);st_y(ii,j-halvewidth:j+halvewidth);st_z(ii,j-halvewidth:j+halvewidth)];
            [vec,val]=eig(D*D');
            [val,pos]=max(diag(val));
            vec=vec(:,pos);
            signal=(vec'*D);
%             phase=angle(signal(1:twohalvewidth)*signal(2:width)');
            
            dist=abs(soffset(j-halvewidth:j+halvewidth-1)-soffset(j-halvewidth+1:j+halvewidth));
            phasediffs=signal(1:twohalvewidth).*conj(signal(2:width));
            magnitude=abs(phasediffs);
            phasediffs=(phasediffs./magnitude).^(15./dist);
            
            phase=angle(sum(phasediffs));
%                phase=angle(sum(sum(D(:,1:twohalvewidth).*conj(D(:,2:width)))));

            vel=2*pi*frequency*15/(sign(soffset(j))*phase+frequency*2*pi/lmovel*15);
            if vel>200 & vel < 750
%                 phaserot=angle(sum(signal.*exp(i*(-halvewidth:1:halvewidth)*phase)));
%                 vec=vec*complex(cos(phaserot),sin(phaserot))
%                 signal=signal*complex(cos(phaserot),sin(phaserot));

                    
                    
%                 vec=mean(D.*([1;1;1]*exp(i*(-halvewidth:1:halvewidth)*phase)),2);
%                 scale=norm(vec);
%                 if (scale ~= 0)
%                     vec=vec/scale;
%                 end
                


%                  temp=vec'*D(:,halvewidth+1);
                %                  temp=signal*exp(-i*(-halvewidth:1:halvewidth)*phase)';
%                 temp=mean(signal.*exp(i*(-halvewidth:1:halvewidth)*phase));


                temp=signal*exp(-i*((soffset(j-halvewidth:j+halvewidth)-soffset(j))/15).*phase)';
                sc=(real(signal(halvewidth+1))*real(temp)+imag(signal(halvewidth+1))*imag(temp))/(temp*conj(temp));
                temp=sc*temp;
                noisefreq_x(freq+1,j)=noisefreq_x(freq+1,j)+vec(1)*(temp);
                noisefreq_y(freq+1,j)=noisefreq_y(freq+1,j)+vec(2)*(temp);
                noisefreq_z(freq+1,j)=noisefreq_z(freq+1,j)+vec(3)*(temp);

%                 noisefreq_x(freq+1,j)=noisefreq_x(freq+1,j)+st_x(ii,j);
%                 noisefreq_y(freq+1,j)=noisefreq_y(freq+1,j)+st_y(ii,j);
%                 noisefreq_z(freq+1,j)=noisefreq_z(freq+1,j)+st_z(ii,j);          
%             
            
            else
                counter=counter+1;
            end
%             rotangle=(atan((2*real(vec(:,pos)')*imag(vec(:,pos)))/(2*norm(imag(vec(:,pos)))^2-1)))/2;
%             vec(:,pos)=vec(:,pos)*complex(cos(rotangle),sin(rotangle));
%             if norm(real(vec(:,pos))) > norm(imag(vec(:,pos)))
%                 ellipticity(ii,j)=norm(imag(vec(:,pos)));
%             else
%                 ellipticity(ii,j)=norm(real(vec(:,pos)));
%             end
%                              test(ii,j)=phase;
            
            
             ellipticity(ii,j)=vel;
%               ellipticity(ii,j)=2*pi*frequency*15/(sign(soffset(j))*phase+frequency*2*pi/lmovel*15);
%             ellipticity(ii,j)=frequency*15/(sign(soffset(j))*phase);
        end
    end
    fprintf('Frequency: %3.2f Hz\n',frequency)
    if rem(frequency,2.5) <= 1/(resamp_nsamples*resamp_sampint)
        figure('windowstyle','docked')
        imagesc(ilmo(ellipticity,offset,lmovel,resamp_sampint*1000))
        str=sprintf('Frequency: %3.2f Hz\n',frequency);
        title(str)
        caxis(1000*[-1 1])
        colorbar
    end
    
    
    
%     imagesc(abs(squeeze(st_z)))
%     clim=caxis;
%     caxis(clim/10);
%     pause

    
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
