%plot OVT's from a cross-spread assuming reciprocity
function [hx,hy,OVT_no,OVT,OVT_col,OVT_row,OVT_flg,fig]=reciprocal_xspread3(Nr,RI,RLI,NRL,Ns,SI,SLI,NSL,fig);
% creates mapping between hx,hy and OVT_col and OVT_row assuming 
%    1) reciprocity  (generalize)
%    2) receiver and source lines are staggered
%    3) a specific cross-spread defined by input parameters
%          assume working with cross-spread centered on 0,0
%          !!! assume symmetric cross-spread  (i.e. fold is evenly divisible by 4)%
%optionally plot OVT's from a cross-spread assuming reciprocity

%INPUT
%    Nr   number of receivers along receiver line
%    RI   receiver interval
%    RLI  receiver line interval
%    NRL  number of receiver lines in a cross-spread
%    Ns  number of sources along source line in cross-spread
%    SI   source interval
%    SLI  source line interval
%    NSL  number of source lines
%
%OUTPUT
%   hx
%   hy
%   OVT_no
%   OVT        structure containing
%       OVT.Hx    Average offset x
%       OVT.Hy    Average offset y
%       OVT.H     average offset = sqrt(Hx^2+Hy^2)
%       OVT.AZM   Azimuth in degrees
%       OVT.col   column number of OVT
%       OVT.row   row number of OVT
%   OVT_col
%   OVT_row
%   OVT_flg:   set if longer of two reciprocal offsets
%   OVT_H
%   OVT_AZM

%arrays used to define OVT_no
col= [ 1 -1  1 -1  2  2  1 -1 -2 -2  2  2  1 -1 -2 -2  3  3  3  2  1 -1 -2 -3 -3 -3 3 3 3 2 1 -1 -2 -3 -3 -3 ];     
row= [ 1  1  1  1  1  2  2  2  2  1  1  2  2  2  2  1  1  2  3  3  3  3  3  3  2  1 1 2 3 3 3  3  3  3  2  1 ];
larg=[ 0  0  1  1  0  0  0  0  0  0  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0 1 1 1 1 1  1  1  1  1  1 ];

%Predicted fold
Mx=Nr*RI/2/SLI;
My=NRL/2;  
M=Mx*My;

%initialize arrays
hx=zeros(M,1);
hy=zeros(M,1);
hhx=zeros(M,1);
hhy=zeros(M,1);
OVT_row=zeros(M,1);
OVT_col=zeros(M,1);
OVT_flg=zeros(M,1);

%  !!! assume symmetric cross-spread  (i.e. fold is evenly divisible by 4)
Ncol=ceil(Mx/2); sline_bias=1;
Nrow=ceil(My/2); xline_bias=1;
if (abs(2*Ncol-Mx)>0)
    error('currently fold must be divisible by 4')
    sline_bias=0;
end
if (abs(2*Nrow-My)>0)  %odd
    error('currently fold must be divisible by 4')
    xline_bias=2;    %for odd fold need asymmetric spread
end

%values use for plotting source line and receiver line used in cross-spread
nRy0=floor(RLI/SI)-1/2;    %used in calculating stagger
Ry0=SI*nRy0;               %stagger 
%y-position of the center of the cross-spread (same as in plot_xspread)
y0=(floor(NRL/2)-xline_bias)*RLI+Ry0;  
nSx0=floor(SLI/RI)-1/2;   %should be odd number
Sx0=RI*nSx0;
%x-position of the center of the cross-spread (same as in plot_xspread)
x0=(floor(NSL/2)-sline_bias)*SLI+Sx0;  %same as in plot_xspread

%setup plotting
flg_plot=0
if (flg_plot)
    figure(fig); hold;
    axis square;
    
    %PLOT receiver line of cross-spread
    n=[1:Nr];
    Rx_xspread=(n-1)*RI;
    
    %PLOT source line of cross-spread
    m=[1:Ns];
    Sy_xspread=(m-1)*SI;   
    
    %center x-spread on origin   y0=x0=0
    %rline_x=[Rx_xspread(1) Rx_xspread(end)]; rline_y=[y0 y0];
    rline_x=[Rx_xspread(1)-x0 Rx_xspread(end)-x0]; rline_y=[0 0];
    plot(rline_x,rline_y,'b','LineWidth',2);  %plot receiver line
    sline_x=[0 0]; sline_y=[Sy_xspread(1)-y0 Sy_xspread(end)-y0];
    plot(sline_x,sline_y,'r','LineWidth',2);  %plot source line
    bold_xlabel('receiver_x(red), CDP_x(green), offset_x(cyan)')
    bold_ylabel('receiver_y(red), CDP_y(green), offset_y(cyan)')
    bold_title('Cross-spread')
    grid minor
end

itrace=1;
for ishot=1:Ns/2;  %  move source keeping x-position constant    
    Sx=0;
    Sy=(ishot-1)*SI-y0;   %shot y-position  
    if (flg_plot) plot(Sx,Sy,'r*','MarkerSize',10);  end;%plot source positions along source line
    %loop over receivers in cross-spread
    for n=1:Nr  %define patch in x-direction
        Ry=0;
        Rx=(n-1)*RI-x0;
        if (flg_plot); plot(Rx,Ry,'bs');  end; %plot receiver positions along receiver line
        cdpx=0.5*(Rx+Sx);
        cdpy=0.5*(Ry+Sy);
        if (flg_plot); plot(cdpx,cdpy,'g.','MarkerSize',15');  end; %plot receiver positions along receiver line
        hx_1=0.5*(Rx-Sx);
        hy_1=0.5*(Ry-Sy); 
        if (flg_plot); plot(hx_1,hy_1,'c.','MarkerSize',15');  end; %plot receiver positions along receiver line
        h_1=sqrt(hx_1*hx_1+hy_1*hy_1);
        OVTc_1=sign(hx_1)*ceil(abs(hx_1/SLI));
        OVTr_1=sign(hy_1)*ceil(abs(hy_1/RLI));
        
        %calculate information for CDP with opposite azimuth (will have
        %different offset).  This will lie on a cross-spread defined by
        if (OVTc_1 > 0)
            OVTc_2=2*OVTc_1-1;
        else
            OVTc_2=2*OVTc_1+1;
        end
        OVTr_2=-2*OVTr_1+1;
        Sx_2=OVTc_2*SLI+Sx;
        Ry_2=OVTr_2*RLI+Ry;
            
        %know cdp(cdpx,cdpy), center coordinates of xspread (Sx_2,Ry_2)
        Rx_2=2*cdpx-Sx_2;
        Sy_2=2*cdpy-Ry_2;
        hx_2=0.5*(Rx_2-Sx_2);
        hy_2=0.5*(Ry_2-Sy_2); 
        h_2=sqrt(hx_2*hx_2+hy_2*hy_2);
                
        %if (0); 
        if (flg_plot); 
            plot(Rx_2,Ry_2,'cs');   %plot receiver positions along receiver line
            plot(Sx_2,Sy_2,'k*','MarkerSize',10);  %plot source positions along source line
            cdpx_2=0.5*(Rx_2+Sx_2);
            cdpy_2=0.5*(Ry_2+Sy_2);
            plot(cdpx_2,cdpy_2,'m.','MarkerSize',15');
            plot(hx_2,hy_2,'k.','MarkerSize',15');  
        end; %plot receiver positions along receiver line
                
        hx(itrace)=hx_1;
        hy(itrace)=hy_1;
        OVT_col(itrace)=OVTc_1;
        OVT_row(itrace)=OVTr_1;
                
        hx(itrace+1)=hx_2;
        hy(itrace+1)=hy_2;
        OVT_row(itrace+1)=OVTr_1;
        OVT_col(itrace+1)=OVTc_1;
        %identify larger of the two offsets
        if (h_2 > h_1); 
            OVT_flg(itrace)=0;
            OVT_flg(itrace+1)=1;       
            
            J=find((OVTc_1 == col) & (OVTr_1 == row) & (0 == larg));
            OVT_no(itrace)=J;
        
            J=find((OVTc_1 == col) & (OVTr_1 == row) & (1 == larg));
            OVT_no(itrace+1)=J;
        else
            OVT_flg(itrace)  =1; 
            OVT_flg(itrace+1)=0;
            
            J=find((OVTc_1 == col) & (OVTr_1 == row) & (1 == larg));
            OVT_no(itrace)=J;
        
            J=find((OVTc_1 == col) & (OVTr_1 == row) & (0 == larg));
            OVT_no(itrace+1)=J;
        end                

        %store reciprocal offset for mean offset calculation
        hhx(itrace)=hx_1;
        hhy(itrace)=hy_1;
        hhx(itrace+1)=-hx_2;  %reciprocal offset
        hhy(itrace+1)=-hy_2;  %reciprocal offset 
        itrace=itrace+2;



        if (0)
            [hx_1 hy_1 OVTc_1 OVTr_1 h_1 atan2(hy_1,hx_1)*180/pi]
            [-hx_2 -hy_2 OVTc_2 OVTr_2 h_2 atan2(-hy_2,-hx_2)*180/pi]
            
            OVTc_2=sign(-hx_2)*ceil(abs(hx_2/SLI));
            OVTr_2=sign(-hy_2)*ceil(abs(hy_2/RLI));        
            
            if not((OVTr_2 == OVTr_1) & (OVTc_2 == OVTc_1))
                error('reciprocal OVT is not correct')
            end
        end
        
    end
end


OVT.Hx=zeros(M,1);
OVT.Hy=zeros(M,1);
OVT.H=zeros(M,1);
OVT.AZM=zeros(M,1);
OVT.col=zeros(M,1);
OVT.row=zeros(M,1);
OVT.flg=zeros(M,1);
for iOVT=1:M
    %mapping from COVID to col,row numbers
    OVT.col(iOVT)=col(iOVT);
    OVT.row(iOVT)=row(iOVT);
    OVT.flg(iOVT)=larg(iOVT);
    
    J=find(OVT_no==iOVT);
    OVT.Hx(iOVT)=mean(hhx(J));
    OVT.Hy(iOVT)=mean(hhy(J));
    OVT.H(iOVT)=sqrt(OVT.Hx(iOVT)*OVT.Hx(iOVT)+OVT.Hy(iOVT)*OVT.Hy(iOVT));
    OVT.AZM(iOVT)=atan2(OVT.Hy(iOVT),OVT.Hx(iOVT)); %*180/pi;    
    if (flg_plot);
        plot(hhx(J),hhy(J),'r.');  %plot receiver positions along receiver line
        plot(hhx(J),hhy(J),'c.');  %plot receiver positions along receiver line
        plot(OVT.Hx(iOVT),OVT.Hy(iOVT),'k.');  
    end; %plot receiver positions along receiver line
end

%[Hx/1000 Hy/1000 H/1000 AZM]

hold off
fig=fig+1   

return