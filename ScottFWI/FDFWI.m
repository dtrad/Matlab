function [ vel ] = FDFWI( D, freq, step, fwave, FDFDfunc, nz, nx, vel0, R, optype, numits, PML_thick, rangevel ,vmin,vmax,vel_true,sx,sz,rx,rz)
% optype 1 = steepest descent
% optype 2 = Gauss Newton
%% Define regularization, starting model
regfac=1e-3;
regfac=1e-1;           %suggest by Scott
reg_func=@(model)Diff0_Reg(model,regfac,nz,nx); 

model = 1./(vel0.^2);

%% Plotting
vel_mean = mean(vel0);


%% Main loop
tic
for n=1:floor(length(freq)/step)
    freqind =  n*step;
    ifend=floor(freq(freqind));
    for m=1:numits 
        freqind = 1 + (n-1)*step : n*step;
        fprintf('n %d m %d freq %d %d %d %d %d  freq %f %f %f %f %f\n',n,m,freqind(1:step),freq(freqind(1):freqind(1)+step-1)); % ,freq(freqind+1),freq(freqind+2),freq(freqind+3),freq(freqind+4));
        % Calculate and apply update
        if optype == 1
            eval_grad = @(model,bound)Gradient( freq(freqind) , fwave(freqind) , 1./sqrt(model) , FDFDfunc , D(:,:,freqind) , R , nz , nx , PML_thick , bound );
            [ model , ft ] = Steepest_descent_optimization( model , eval_grad , reg_func );
        elseif optype == 2
            eval_grad = @(model)Gradient( freq(freqind) , fwave(freqind) , 1./sqrt(model) , FDFDfunc , D(:,:,freqind) , R , nz , nx , PML_thick);
            eval_Hess = @(model)Hessian( freq(freqind) , fwave(freqind) , 1./sqrt(model) , FDFDfunc , D(:,:,freqind) , R , nz , nx , PML_thick );
            [ model , ft ] = Gauss_Newton_optimization( model , eval_grad , eval_Hess , reg_func );
        end
        
        vel = 1./sqrt(model);
        
        for ii=1:length(model)
            if (model(ii)<0) fprintf('ii %d model %e\n',ii,model(ii));
            end
        end
        
        % Plot new model
        %{
        ifig = floor(1000+n);
        figure(ifig);
        imagesc(reshape(vel,nz,nx));
        caxis([vmin,vmax]);
        txt = sprintf('Updated model. Iteration %d freq %f to %f',n,freq(freqind(1)),freq(freqind(1)+step-1));
        title(txt);
        drawnow
        %}
        
        %Estimate time remaining
        fprintf('iteration %d of %d \n',(n - 1) * numits + m,floor(length(freq)/step)*numits)
        
        if n~=length(freq)
            t_cur = toc;
            t_est = t_cur*((numits*floor(length(freq)/step)) - ((n - 1) * numits + m)) /((n - 1) * numits + m);
            t_h = floor(t_est/3600);
            t_m = floor((t_est-3600*t_h)/60);
            t_s = round(t_est - 3600*t_h - 60*t_m);
            
            if t_h~=0
                fprintf('Estimated time remaining %i hours, %i minutes and %i seconds \n',t_h,t_m,t_s)
            elseif t_m~=0
                fprintf('Estimated time remaining %i minutes and %i seconds \n',t_m,t_s)
            else
                fprintf('Estimated time remaining %i seconds \n',t_s)
            end
        end
        
        %% Plot profiles
% close all;
        ifig = floor(1000+n);
        figure(ifig);
profileindex = floor(nx*0.5);  %25;
vel0mat = reshape(vel0,nz,nx);
vel0profile = vel0mat(:,profileindex)
velmat = reshape(vel,nz,nx);
velprofile = velmat(:,profileindex);
veltruemat = reshape(vel_true,nz,nx);
veltrueprofile = veltruemat(:,profileindex);
          
t_cur = toc;
t_h = floor(t_cur/3600);
t_m = floor((t_cur-3600*t_h)/60);
t_s = round(t_cur - 3600*t_h - 60*t_m);


if (optype==2)
 txt = sprintf('Gauss-Newton: Freq loop %d, steps 1 to %d',n,ifend);
else
 txt = sprintf('Steepest Descend: Freq loop %d,  steps 1 to %d',n,ifend);  
end    
txt_time=sprintf('Elapsed time: %i minutes and %i seconds',t_m,t_s);

close(ifig);
figure(ifig);
x0=10;
y0=10;
width=850;
height=400
set(gcf,'units','points','position',[x0,y0,width,height])

subplot(2,2,[1 1]), 
plot(1:nz,velprofile,'k-',1:nz,veltrueprofile,'b-',1:nz,vel0profile,'r+');
title({txt,txt_time});
ylim([1800 3400]);
legend('Fina; velocity','True velocity','Starting velocity');
grid; grid minor;
%txt=sprintf('Freq: %f to %f by %f',f1,endfreq(numbands),finc);

subplot(2,2,2),
imagesc(velmat); %colormap('gray'); 
caxis([vmin,vmax]);
colorbar;
title('True velocity');

subplot(2,2,3),
imagesc(vel0mat); %colormap('gray'); 
caxis([vmin,vmax]);
colorbar;
hold on;
plot(sx,sz,'ko');
plot(rx,rz,'k+');
caxis([vmin,vmax]);
title('Starting velocity');

subplot(2,2,4),
imagesc(veltruemat); %colormap('gray'); 
colorbar;
caxis([vmin,vmax]);
title('True velocity');

drawnow;

    
    end
   
end


end

