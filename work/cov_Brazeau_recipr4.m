clear;
close all;
%terminology
%    B   bin dimension
%    Br   in-line bin dimension (receiver direction)
%    Bs   cross-line bin dimension (source direction)
%    Lr   in-line dimension of the patch
%    Ls   cross-line dimension of the patch
%    NC   number of channels
%    Nr   number of receivers along inline in patch
%    NRL  number of receiver lines
%    NSL  number of source lines
%    RI   receiver interval
%    RLI  receiver line interval
%    SD   source density
%    SI   source interval
%    SLI  source line interval

%Assume receivers laid out in x-direction 
%Assume split spread with staggered shot, shot line is always in the middle of the patch
%R     R     R     R     R     R
%               S
%                 c  c  c
%Assume roll in /out
%fixed parameters
xline_offset=1000;  %multiplier used to calculate CDP
fig=1;              %1st plot
Vrms1=3000;
Vrms2=2700;
phi_eff=0;
t0=1.2;
t0_sq=t0*t0;

%reflectivity parameters
t0=1.5;
A=.1
Biso=-.2
Bani=.1
phi_iso=45/180*pi;
flg_refl=1;

flg_recipr=1
%%variables
RI=60         %Group interval
%RI=300
%SLI=300       %should be integer multiplier of RI
SLI=600       %should be integer multiplier of RI
SI=60         %Source interval
%SI=300
%RLI=300       %should be integer multiplier of SI
RLI=600       %should be integer multiplier of SI

%Parameters defining Source and Receiver Lines in cross-spread
NSL=8;           %number of source lines in the patch
NRL=8;            %number of receiver lines
%NSL=10;           %number of source lines in the patch
%NRL=10;            %number of receiver lines
%NSL=20;           %number of source lines in the patch
%NRL=20;            %number of receiver lines
%calculated
Nr=SLI/RI*NSL;  %number of receivers in cross-spread
Ns_per_xspread=RLI/SI*NRL

%Survey Parameters
NSLs=8         %number of source lines in survey, want at least Nr*RI/SLI to get full fold
%NSLs=10         %number of source lines in survey, want at least Nr*RI/SLI to get full fold
%NSLs=20         %number of source lines in survey, want at least Nr*RI/SLI to get full fold
NS_per_sline=(RLI/SI)*NSL  %note this assumes staggered shots
%NS_per_sline=100  %note this assumes staggered shots
%calculated
NC=Nr*NRL             %number of channels used in survey
NS=NS_per_sline*NSLs  %number of sources in the survey

%check using enough sources/(source line) and source lines to get full fold
if (NSLs < NSL)
    warning('for fold fold need at least Nr*RI/SLI source lines')
    [NSLs NSL]
    break
end

if (NS_per_sline < Ns_per_xspread)
    warning('for fold fold need at least RLI/SI*NRL shots in shot line')
    [NS_per_sline  Ns_per_xspread]
    break
end

%Ideal fold
Mx=Nr*RI/2/SLI;
My=NRL/2;  %My=NRL*RLI/2/RLI
Mfold=Mx*My

%create spiral mapping
[spiralx, spiraly, My_odd, Mx_odd]=spiral_ver1(Mx,My);
sline_bias=1; if (Mx_odd); sline_bias=0; end

%define BIN size
Bs=SI/2;
Br=RI/2;


% patch dimensions and roll-along
Lr=(Nr-1)*RI; %inline dimension of the patch
Ls=(NRL-1)*RLI;  %cross-line dimension of the patch
N_Rollx=SLI/RI;
N_Rolly=RLI/SI;  

%plot cross-spread
fig=plot_xspread(Nr,RI,RLI,NRL,SI,SLI,NSL,Ns_per_xspread,fig)
%fig=reciprocal_xspread3(Nr,RI,RLI,NRL,SI,SLI,NSL,Ns_per_xspread,fig);
[hx,hy,OVT_no,OVT,OVT_col,OVT_row,OVT_flg,fig]=reciprocal_xspread3(Nr,RI,RLI,NRL,Ns_per_xspread,SI,SLI,NSL,fig);

%Sx0  offset in x-dirrection of 1st receiver line
nSx0=floor(SLI/RI)-1/2;   %should be odd number
Sx0=RI*nSx0;
%Ry0  offset in x-dirrection of 1st receiver line
nRy0=floor(RLI/SI)-1/2;  %should be odd number
Ry0=SI*nRy0;
%TODO could write in terms of CDP then nSx0 and nRy0 will be integer
%(modify plot_xspread as well)

%reflectivity calculations
if (flg_refl)
    B=Biso+Bani/2;
    C=Bani*cos(2*phi_iso);  %where phi_sym=phi0=phi_iso+90
    D=Bani*sin(2*phi_iso);
    x=[A B C D]';
end

%size of plot area
YY=Ns_per_xspread*SI+NRL*RLI;
XX=(Nr*RI)+SLI*NSL;
YX=max(YY,XX);

%initialize array
Ntrace=NS*NC;
%trace.no=zeros(1,Ntrace);
%trace.chan=zeros(1,Ntrace);
trace.sn=zeros(1,Ntrace);
%trace.sx=zeros(1,Ntrace);
%trace.sy=zeros(1,Ntrace);
%trace.rx=zeros(1,Ntrace);
%trace.ry=zeros(1,Ntrace);
trace.cdpx=zeros(1,Ntrace);
trace.cdpy=zeros(1,Ntrace);
trace.hx=zeros(1,Ntrace);
trace.hy=zeros(1,Ntrace);
trace.h=zeros(1,Ntrace);
trace.azimuth=zeros(1,Ntrace);
trace.H=zeros(1,Ntrace);
trace.AZIMUTH=zeros(1,Ntrace);
trace.inline=zeros(1,Ntrace);
trace.xline=zeros(1,Ntrace);
trace.cdp=zeros(1,Ntrace);
trace.cov=zeros(1,Ntrace);
trace.cov_row=zeros(1,Ntrace);
trace.cov_col=zeros(1,Ntrace);

nshot=0;  %shot counter
itrace=0; %trace counter
for ishot=1:NS_per_sline;  %  move source keeping x-position constant    
    %  move source line keeping source y-position constant, this will
    %  simulate a 2D geometry,  need to move patch in conjuction with
    %  source
    for iSL=1:NSLs;
        nshot=nshot+1;
        Sx=(iSL-1)*SLI+(floor(NSL/2)-sline_bias)*SLI+Sx0;  %same as in plot_xspread
        
        iRolly=floor((ishot-1)/N_Rolly);
        iRolly_remainder=ishot-iRolly*N_Rolly;
       
        Sy=2*floor(NRL/4)*RLI+(ishot-1)*SI;   %shot y-position                                
        %[ishot Sx Sy]
        
        %TODO try to use same subroutine to calculate patch as plot_xspread
        %define receiver patch (shot centered)
        ichan=0;
        for iRL=1:NRL;  %define patch in y-direction
           Ry=RLI*(iRolly+(iRL-1))+Ry0;  
           for n=1:Nr  %define patch in x-direction
              ichan=ichan+1;
              itrace=itrace+1;
              Rx=(iSL-1)*SLI+(n-1)*RI;
                     
              %trace.no(itrace)=itrace;
              %trace.chan(itrace)=ichan;
              trace.sn(itrace)=nshot;
              %trace.sx(itrace)=Sx;
              %trace.sy(itrace)=Sy;
              %trace.rx(itrace)=Rx;
              %trace.ry(itrace)=Ry;
              trace.cdpx(itrace)=0.5*(Rx+Sx);
              trace.cdpy(itrace)=0.5*(Ry+Sy);
              trace.hx(itrace)=0.5*(Rx-Sx);
              trace.hy(itrace)=0.5*(Ry-Sy); 
              trace.h(itrace)=sqrt(trace.hx(itrace)*trace.hx(itrace)+trace.hy(itrace)*trace.hy(itrace));
              trace.azimuth(itrace)=atan2(trace.hy(itrace),trace.hx(itrace));
              
              inline=floor(trace.cdpx(itrace)/Br)+1;
              xline=floor(trace.cdpy(itrace)/Bs)+1;
              trace.inline(itrace)=inline;
              trace.xline(itrace)=xline;
              trace.cdp(itrace)=xline*xline_offset+inline;
          
              if (0)  %check to see my CDP logic is correct
                cdpx=(trace.inline(itrace)-0.5)*Br;
                if (abs(cdpx-trace.cdpx(itrace))>Br/2); 
                    warning('problem with calculation of inline'); 
                end
                cdpy=(trace.xline(itrace)-0.5)*Bs;
                if (abs(cdpy-trace.cdpy(itrace))>Bs/2); 
                    warning('problem with calculation of xline'); 
                end
              end
          
              %cov_col=floor(trace.hx(itrace)/SLI);
              %cov_row=floor(trace.hy(itrace)/RLI);
              J=find((trace.hx(itrace) == hx) & (trace.hy(itrace) == hy));
              if not(length(J) == 1)
                error('did not find unique mapping between hx,hy and COVID')
              end
              iOVT=OVT_no(J);
              trace.cov(itrace)=iOVT;
              trace.cov_col(itrace)=OVT.col(iOVT);
              trace.cov_row(itrace)=OVT.row(iOVT);
              trace.H(itrace)=OVT.H(iOVT);
              trace.AZIMUTH(itrace)=OVT.AZM(iOVT);              
           end
        end
        if (0)
            I=find(trace.sn==nshot);
            Sx=trace.cdpx(I(1))-trace.hx(I(1));
            Sy=trace.cdpy(I(1))-trace.hy(I(1));
            Rx=trace.cdpx(I)+trace.hx(I);
            Ry=trace.cdpy(I)+trace.hy(I);

            figure(201); bold_title('source&receivers locations')
            plot(Rx,Ry,'*'); axis([0 YX 0 YX])
            hold; plot(Sx,Sy,'r*')
            plot(trace.cdpx(I),trace.cdpy(I),'k*');
            hold off;
            [Sx Sy]
        end
    end
end

if (0)
    movie2avi(F1,'COV_shots.avi','compression','None','fps',2)
    break
end

if (0)
    figure(2); bold_title('offset')
    plot(trace.hx,trace.hy,'k*'); 
end

%sort by CDP and create fold map
inline1=min(trace.inline); inlineN=max(trace.inline);
xline1=min(trace.xline); xlineN=max(trace.xline);
%COV_fold=zeros(xline,inline,Mfold);
%COV_h=zeros(xline,inline,Mfold);
%COV_dh=zeros(xline,inline,Mfold);
%COV_fh=zeros(xline,inline,Mfold);
%COV_az=zeros(xline,inline,Mfold);
%COV_daz=zeros(xline,inline,Mfold);
%COV_dt=zeros(xline,inline,Mfold);
COV_Rh=zeros(xline,inline,Mfold);
COV_RH=zeros(xline,inline,Mfold);
COV_dR=zeros(xline,inline,Mfold);
%COV_theta=zeros(xline,inline,Mfold);
COV_A=zeros(xline,inline);
COV_Biso=zeros(xline,inline);
COV_Bani=zeros(xline,inline);
COV_Phi_iso=zeros(xline,inline);
for xline=xline1:xlineN;
    for inline=inline1:inlineN;
        iCDP=xline*xline_offset+inline;
        I=find(trace.cdp==iCDP);
        Nfold=length(I);
        fold(xline,inline)=Nfold;
        GG=zeros(Nfold,4);
        data=zeros(Nfold,1);
        for k=1:Nfold        
            cov_number=trace.cov(I(k));
            %COV_fold(xline,inline,cov_number)=COV_fold(xline,inline,cov_number)+1;
            
            %h=sqrt((spiralx(cov_number)*SLI)^2+(spiraly(cov_number)*RLI)^2);
            %[h trace.h(I(k))]
            %COV_h(xline,inline,cov_number)=trace.h(I(k)); %only makes sense since single fold
            %COV_dh(xline,inline,cov_number)=trace.h(I(k))-trace.H(I(k)); %only makes sense since single fold
            %COV_fh(xline,inline,cov_number)=(trace.h(I(k))-trace.H(I(k)))/trace.H(I(k)); %only makes sense since single fold

            if (0)
                COV_az(xline,inline,cov_number)=trace.azimuth(I(k)); %only makes sense since single fold
                COV_daz(xline,inline,cov_number)=trace.azimuth(I(k))-trace.AZIMUTH(I(k)); %only makes sense since single fold
            
                Vnmo_sq=azimuthal_Vnmo(Vrms1,Vrms2,phi_eff,trace.azimuth(I(k)));
                th=sqrt(t0_sq+4*trace.H(I(k))*trace.H(I(k))/Vnmo_sq);
            
                Vnmo_sq=azimuthal_Vnmo(Vrms1,Vrms1,phi_eff,trace.AZIMUTH(I(k)));
                tH=sqrt(t0_sq+4*trace.H(I(k))*trace.H(I(k))/Vnmo_sq);
                dt=th-tH;
                COV_dt(xline,inline,cov_number)=dt;
            end
            
            if (flg_refl)
                %ASSUMING VRMS1=V0!!!!!!!!!!
                %theta=atan2(t0*Vrms1/2,trace.h(I(k)));   
                theta=atan(2*trace.h(I(k))/(t0*Vrms1));   
                G=avo_HTI_kernal(theta,trace.azimuth(I(k)));
                %GG(k,:)=G;
                COV_Rh(xline,inline,cov_number)=G*x;
                data(k)=COV_Rh(xline,inline,cov_number);
                
                theta=atan(2*trace.H(I(k))/(t0*Vrms1));   
                %COV_theta(xline,inline,cov_number)=theta;
                G=avo_HTI_kernal(theta,trace.AZIMUTH(I(k)));
                GG(k,:)=G;
                COV_RH(xline,inline,cov_number)=G*x;
                COV_dR(xline,inline,cov_number)=COV_RH(xline,inline,cov_number)-COV_Rh(xline,inline,cov_number);
            end
        end
        if ((Nfold > 3) & cond(GG) <= 10^3)
            y=inv(GG'*GG)*GG'*data;
        else
            y=zeros(4,1);
        end
        C=y(3); D=y(4);
        Bani_est=sqrt(C^2+D^2);
        phi_iso_est=atan2(D,C)/2;
        COV_A(xline,inline)=y(1);
        COV_B(xline,inline)=y(2)-Bani_est/2;
        COV_Bani(xline,inline)=Bani_est;
        COV_Phi_iso(xline,inline)=phi_iso_est*180/pi;
    end
end
figure(100); imagesc(fold); colorbar; flipy
if (0)
    figure(110); plot(squeeze(COV_dt(200,200,:))*1000); bold_xlabel('OVT number (snail plot)'); bold_ylabel('time (msec)');  bold_title('time error due to OVT azimuth dispersion (msec)')
    figure(111); imagesc(squeeze(max(COV_dt))*1000); colorbar; bold_xlabel('OVT number (snail plot)'); bold_ylabel('in-line');  bold_title('Max. time error due to OVT azimuth dispersion (msec)')
end
if (1)
    figure(112); imagesc(COV_A,[0 .118]); colorbar; bold_xlabel('in-line'); bold_ylabel('cross-line');  bold_title('Intercept')
    figure(113); imagesc(COV_B,[-.29 0]); colorbar; bold_xlabel('in-line'); bold_ylabel('cross-line');  bold_title('Isotropic Gradient')
    figure(114); imagesc(COV_Bani,[0 .25]); colorbar; bold_xlabel('in-line'); bold_ylabel('cross-line');  bold_title('Anisotropic Gradient')
    figure(115); imagesc(COV_Phi_iso,[-90 90]); colorbar; bold_xlabel('in-line'); bold_ylabel('cross-line');  bold_title('Isotropy plane')
end
%figure(102);
%scrsz = get(0,'ScreenSize');
%figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])

for k=1:Mfold
    %fractional error in offset
    %figure(102); imagesc(COV_dt(:,:,k),[-.005,.005]); flipy; axis square
    %figure(102); imagesc(COV_fh(:,:,k),[-1,1]); flipy; axis square
    %figure(103); imagesc(COV_h(:,:,k)*2,[-1,1]); colorbar; axis
    
    %figure(103); imagesc(COV_h(:,:,k)*180/pi,[0,230]); colorbar; axis square

    %figure(102); imagesc(COV_az(:,:,k)*180/pi,[-180,180]); flipy; axis square
    %figure(102); imagesc(COV_daz(:,:,k)*180/pi,[-45,45]); flipy; axis square
    %figure(102); imagesc(COV_Rh(:,:,k),[0,.15]); flipy; axis square
    %figure(102); imagesc(COV_RH(:,:,k),[0,.15]); flipy; axis square
    figure(102); imagesc(COV_dR(:,:,k),[-.03,.03]); flipy; axis square
    %figure(102); imagesc(COV_theta(:,:,k)*180/pi,[0,40]); flipy; axis square
    F(k) = getframe;
end

if (1)
    % Play the movie ten times
    fps=2
    movie(F,1,fps)
end

%movie2avi(F,'COV_RH.avi','compression','None','fps',2)
%movie2avi(F,'COV_theta_H.avi','compression','None','fps',2)
movie2avi(F,'COV_dR.avi','compression','None','fps',2)

%save COV_SLI300
%save COV_SLI600
break