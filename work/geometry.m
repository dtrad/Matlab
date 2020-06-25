% Program geometry
% testing geometry building for interpolation given full cmp/azimuth tables.
% Given requirements for azimuth/offset tables at each cmp
% the shot and receiver locations are calculated.
% Basic design parameters are:
%                   nx: number of inlines
%                   ny: number of crosslines 
%                   dx: sampling interval for inlines
%                   dy: sampling interval for crosslines
%                   dh: sampling interval for offset
%                   da: sampling interval for azimuth
%                   maxh: maximum offset
%
% Daniel Trad- Veritas. Jan 2005
%
% needs functions:
%                countfold
%                getshift
%                getUniqueID
%                getXYfromID
%                getShotPatch
%************************************************************************** 

clear;
close all;
nfig =1; % comment out this line to keep previous plots

% survey parameters
% cmp grid
nx = 100; % number of inlines and crosslines (cmpx and y)
ny = 100;
dx = 25; % bin size
dy = 25;

% azimuth and offset axis
da = 2*pi/12; % assume reciprocity, varies between 0-180 degrees only
dh = 120;  % offset sampling
maxh = 2400; % maximum offset

% define staggering for offset as a function of cmp and azimuth
offsetShift=0;%60; 
% define staggering for azimuth as a function of cmp and offset
azimuthShift=0;

%binning size for shots and rcvrs
binSx = dx/2;
binSy = dy/2;
binRx = dx/2;
binRy = dy/2;

% start defining cmp axes.
cmpx = 1:nx;
cmpx = (cmpx-1)*dx;

cmpy = 1:nx;
cmpy = (cmpy-1)*dy;

%create grid of cmps
cmpgridx = ones(ny,1)*cmpx(1:nx);cmpgridx = cmpgridx(:);
cmpgridy = cmpy(1:ny)'*ones(1,nx);cmpgridy = cmpgridy(:);

% azimuth
azimuth = -pi+da:da:pi;
na = length(azimuth);

% offsets
offset = dh/2:dh:maxh;
nh = length(offset);

%offset azimuth table.
% The number of cmp families depends on whether getShift depends on cmp or
% not. If it doesn't there is only one family.
for icmp=1:2 % here use as many families of cmps as required (normally 1 or 2)
    for ia=1:na
        for ih=1:nh
                imap =(ia-1)*nh+ih;
                [shifth,shiftA] = getShift(icmp,ia,ih,offsetShift,azimuthShift);
                offsetList(icmp,imap)=offset(ih)+shifth;
                azimuthList(icmp,imap)=(azimuth(ia)+shiftA)*180/pi;
        end
    end
end
cmpfold = na*nh;


    
% Calculate shot and rcvr positions given the previous 
% cmpx,cmpy,offset and azimuth values.
% Create traces for every cmpxy, azimuth and offset;
ntraces = nx*ny*na*nh;
shotx = zeros(1,ntraces);
shoty = zeros(1,ntraces);
rcvrx = zeros(1,ntraces);
rcvry = zeros(1,ntraces);

for ix=1:nx    
    for iy=1:ny     
        for ia=1:na
            for ih=1:nh
                imap = (ix-1)*ny*na*nh + (iy-1)*na*nh+(ia-1)*nh+ih;
                icmp = (ix-1)*ny +iy;
                [shifth,shiftA] = getShift(icmp,ia,ih,offsetShift,azimuthShift);

                offsetx = (shifth + offset(ih))*cos(azimuth(ia)+shiftA);
                offsety = (shifth + offset(ih))*sin(azimuth(ia)+shiftA);
                
                shotx(imap) = cmpx(ix) + offsetx/2;
                shoty(imap) = cmpy(iy) + offsety/2;
                rcvrx(imap) = cmpx(ix) - offsetx/2;
                rcvry(imap) = cmpy(iy) - offsety/2;
                
                
            end
        end
    end
end

% bin shot and rcvr coordinates
shotBinX = round(shotx/binSx)*binSx;
shotBinY = round(shoty/binSy)*binSy;
rcvrBinX = round(rcvrx/binRx)*binRx;
rcvrBinY = round(rcvry/binRy)*binRy;

%count number of shots and receivers and fold
%create a single number with both shotx and y
[shotUniqueId,rangey,minx,miny] = getUniqueId(shotBinX,shotBinY);

%shotxy= shotBinX*1e9 + shotBinY;
[shotFold,shotIDNew,ix] = countfold(shotUniqueId);
% unfold shotnew in shot x andy
[shotNewX,shotNewY] = getXYfromID(shotIDNew,rangey,minx,miny);

% plot a few shot patches.
for N=1:200:1000
    [rpatchX,rpatchY,thisShotX,thisShotY,shotFold]=getShotPatch(shotUniqueId,shotBinX,shotBinY,rcvrBinX,rcvrBinY,shotUniqueId(N));
    shotFold; 
    figure(nfig)
    plot(rpatchX,rpatchY,'o',thisShotX,thisShotY,'+');text= sprintf('shot patch ID = %d with fold = %d',shotUniqueId(N),shotFold);title(text); 
    nfig= nfig+1
end


% to plot fold, need to create a sparse matrix
[foldMap]=createSparseMatrix(shotNewX,shotNewY,shotFold);

%plots
figure(nfig);plot(cmpgridx,cmpgridy,'+');title('original cmp grid');
nfig=nfig+1;

figure(nfig);plot(offsetList(1,:),azimuthList(1,:),'o');title('azimuth-offset table for odd cmp');
xlabel('offset');ylabel('azimuth');
nfig=nfig+1;

figure(nfig);plot(offsetList(2,:),azimuthList(2,:),'o');title('azimuth-offset table for even cmp');
xlabel('offset');ylabel('azimuth');
nfig=nfig+1;

figure(nfig);plot(shotx,shoty,'.');title('actual shots');
xlabel('shotx');ylabel('shoty');
nfig=nfig+1;

figure(nfig);plot(shotBinX,shotBinY,'.');title('shots after binning on a regular grid');
xlabel('shotx');ylabel('shoty');
nfig=nfig+1;

% this figure seems to require too much memory
% save data for later plotting
save foldmap.mat foldMap
%figure(nfig);imagesc(foldMap);colorbar;title('fold per shot');
%nfig=nfig+1;

%figure(nfig);plot3(shotNewX,shotNewY,shotFold,'o');title('fold per shot');
%xlabel('shotx');ylabel('shoty');zlabel('fold')
%nfig=nfig+1;

% define axes to plot for zoom effect
ax = [560 640 560 640];

figure(nfig);plot(shotx,shoty,'o');title('actual shots (zoomed)');axis(ax);
nfig=nfig+1;

figure(nfig);plot(shotBinX,shotBinY,'o');title('shots after binning on a 25x25 grid (zoomed)');axis(ax);
nfig=nfig+1;

% Actual midpoints after binning
offsetFinalX = zeros(1,ntraces);
offsetFinalY = zeros(1,ntraces);
cmpFinalX = zeros(1,ntraces);
cmpFinalY = zeros(1,ntraces);
offsetFinal = zeros(1,ntraces);
azimuthFinal = zeros(1,ntraces);

for imap = 1:ntraces
    offsetFinalX(imap)=shotBinX(imap)-rcvrBinX(imap);
    offsetFinalY(imap)=shotBinY(imap)-rcvrBinY(imap);
    cmpFinalX(imap)=(shotBinX(imap)+rcvrBinX(imap))/2;
    cmpFinalY(imap)=(shotBinY(imap)+rcvrBinY(imap))/2;
    offsetFinal(imap)=sqrt(offsetFinalX(imap)^2+offsetFinalY(imap)^2);
end
%display('start azimuth calculation');

for imap = 1:ntraces
    if (offsetFinalX(imap)~=0)
        azimuthFinal(imap)=(atan(offsetFinalY(imap)/offsetFinalX(imap)))*180/pi;
    elseif (offsetFinalY(imap) > 0)
        azimuthFinal(imap)=90;
    else
        azimuthFinal(imap)=-90;    
    end  
end
    
figure(nfig);plot(cmpFinalX,cmpFinalY,'+');title('midpoints after binning');
nfig = nfig + 1;

figure(nfig);plot(offsetFinal,azimuthFinal,'+');title('azimuth-offset table after binning');
nfig = nfig + 1;

