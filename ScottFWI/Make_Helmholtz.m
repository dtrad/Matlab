function [ A , M , N , S , R , L ] = Make_Helmholtz( vel_in , vel0 , omega , PML_thick , nz , nx , dz , dx , R_theor )
%% Rotate to correct dimension, define important terms
nzPML = nz + 2*PML_thick;
nxPML = nx + 2*PML_thick;
NPML = nzPML * nxPML;

%% Add velpad part

vel_in = reshape(vel_in, nz, nx);
vel = zeros(nzPML, nxPML);
vel(PML_thick+1:PML_thick+nz, PML_thick+1:PML_thick+nx) = vel_in;

% Edges
vel(1:PML_thick, PML_thick+1 : PML_thick+nx) = ones(PML_thick,1)*vel_in(1, :); 
vel(1 + PML_thick + nz : end, PML_thick+1 : PML_thick+nx) = ones(PML_thick,1)*vel_in(nz, :); 

vel(:, 1:PML_thick) = vel(:, PML_thick + 1) * ones(1,PML_thick);
vel(:, 1+PML_thick+nx : end) = vel(:, PML_thick+nx) * ones(1,PML_thick);

vel = vel(:);

gammax = zeros(size(vel));
gammaz = zeros(size(vel));

%% PML implementation
LFT = 1:nzPML;
RIT = NPML-nzPML+1:NPML;
TOP = 1:nzPML:NPML;
BOT = nzPML:nzPML:NPML;

LFTTOP = 1;
RITTOP = NPML-nzPML+1;
LFTBOT = nzPML;
RITBOT = NPML;

PMLx = (PML_thick-1)*dx;
PMLz = (PML_thick-1)*dz;
% R_theor = 0;
if nargin < 9
    R_theor = 1e-3;
end

vmax = max(vel0);
xiL0 = log(1/R_theor)*(3*vmax)./(2*PMLx);
xiR0 = log(1/R_theor)*(3*vmax)./(2*PMLx);
xiT0 = log(1/R_theor)*(3*vmax)./(2*PMLz);
xiB0 = log(1/R_theor)*(3*vmax)./(2*PMLz);
% xiL0 = log(1/R_theor)*(3*vel(LFT))./(2*PMLx);
% xiR0 = log(1/R_theor)*(3*vel(RIT))./(2*PMLx);
% xiT0 = log(1/R_theor)*(3*vel(TOP))./(2*PMLz);
% xiB0 = log(1/R_theor)*(3*vel(BOT))./(2*PMLz);

% if nargin<9
    for n=1:PML_thick
        gammax(LFT+(n-1)*nzPML) = gammax(LFT+(n-1)*nzPML) + xiL0.*((PML_thick-(n-1))/(PML_thick-1)).^2;
        gammax(RIT-(n-1)*nzPML) = gammax(RIT-(n-1)*nzPML) + xiR0.*((PML_thick-(n-1))/(PML_thick-1)).^2;
        gammaz(TOP+(n-1)) = gammaz(TOP+(n-1)) + xiT0.*((PML_thick-(n-1))/(PML_thick-1)).^2;
        gammaz(BOT-(n-1)) = gammaz(BOT-(n-1)) + xiB0.*((PML_thick-(n-1))/(PML_thick-1)).^2;
    end
% end

xix = 1 - (1i*gammax)./omega;
xiz = 1 - (1i*gammaz)./omega;

%% Make Helmholtz matrix
deltax = 1./(xix);
deltaz = 1./(xiz);

delta_pad_x = zeros(nzPML+2,nxPML+2);
delta_pad_x(2:end-1,2:end-1) = reshape(deltax,nzPML,nxPML);
delta_pad_x(1,:) = delta_pad_x(2,:);
delta_pad_x(end,:) = delta_pad_x(end-1,:);
delta_pad_x(:,1) = delta_pad_x(:,2);
delta_pad_x(:,end) = delta_pad_x(:,end-1);

deltaxplus = 0.5*(delta_pad_x(2:end-1,2:end-1)+delta_pad_x(2:end-1,3:end));
deltaxplus = deltaxplus(:);
deltaxminus = 0.5*(delta_pad_x(2:end-1,2:end-1)+delta_pad_x(2:end-1,1:end-2));
deltaxminus = deltaxminus(:);

delta_pad_z = zeros(nzPML+2,nxPML+2);
delta_pad_z(2:end-1,2:end-1) = reshape(deltaz,nzPML,nxPML);
delta_pad_z(1,:) = delta_pad_z(2,:);
delta_pad_z(end,:) = delta_pad_z(end-1,:);
delta_pad_z(:,1) = delta_pad_z(:,2);
delta_pad_z(:,end) = delta_pad_z(:,end-1);

deltazplus = 0.5*(delta_pad_z(2:end-1,2:end-1)+delta_pad_z(3:end,2:end-1));
deltazplus = deltazplus(:);
deltazminus = 0.5*(delta_pad_z(2:end-1,2:end-1)+delta_pad_z(1:end-2,2:end-1));
deltazminus = deltazminus(:);

%%

M = (omega.^2)./vel.^2 - (1./((dx^2)*xix)).*(deltaxplus + deltaxminus) - (1./((dz^2)*xiz)).*(deltazplus + deltazminus) ;

R = (deltaxplus)./((dx^2)*xix);

L = (deltaxminus)./((dx^2)*xix);

S = (deltazplus)./((dz^2)*xiz);

N = (deltazminus)./((dz^2)*xiz);

NL = zeros(size(M));

NR = zeros(size(M));

SL = zeros(size(M));

SR = zeros(size(M));

% %% Clayton 3
% for n=1:length(BOT)
%     ind = BOT(n);
%     M(ind) = 1i*omega^3/vel(ind) + omega^2/(dz*xiz(ind)) - 3i*omega*vel(ind)/(2*(dx*xix(ind))^2) - vel(ind)^2/(2*(dz*xiz(ind))*(dx*xix(ind))^2);
%     N(ind) = -(omega^2)/(dz*xiz(ind)) + vel(ind)^2/(2*(dz*xiz(ind))*(dx*xix(ind))^2);
%     R(ind) = 3i*omega*vel(ind)/(4*(dx*xix(ind))^2) + vel(ind)^2/(4*(dz*xiz(ind))*(dx*xix(ind))^2);
%     L(ind) = 3i*omega*vel(ind)/(4*(dx*xix(ind))^2) + vel(ind)^2/(4*(dz*xiz(ind))*(dx*xix(ind))^2);
%     NR(ind) = -vel(ind)^2/(4*(dz*xiz(ind))*(dx*xix(ind))^2);
%     NL(ind) = -vel(ind)^2/(4*(dz*xiz(ind))*(dx*xix(ind))^2);
% end
% 
% for n=1:length(TOP)
%     ind = TOP(n);
%     M(ind) = 1i*omega^3/vel(ind) + omega^2/(dz*xiz(ind)) - 3i*omega*vel(ind)/(2*(dx*xix(ind))^2) - vel(ind)^2/(2*(dz*xiz(ind))*(dx*xix(ind))^2);
%     S(ind) = -(omega^2)/(dz*xiz(ind)) + vel(ind)^2/(2*(dz*xiz(ind))*(dx*xix(ind))^2);
%     R(ind) = 3i*omega*vel(ind)/(4*(dx*xix(ind))^2) + vel(ind)^2/(4*(dz*xiz(ind))*(dx*xix(ind))^2);
%     L(ind) = 3i*omega*vel(ind)/(4*(dx*xix(ind))^2) + vel(ind)^2/(4*(dz*xiz(ind))*(dx*xix(ind))^2);
%     SR(ind) = -vel(ind)^2/(4*(dz*xiz(ind))*(dx*xix(ind))^2);
%     SL(ind) = -vel(ind)^2/(4*(dz*xiz(ind))*(dx*xix(ind))^2);
% end
% 
% for n=1:length(RIT)
%     ind = RIT(n);
%     M(ind) = 1i*omega^3/vel(ind) + omega^2/(dx*xix(ind)) - 3i*omega*vel(ind)/(2*(dz*xiz(ind))^2) - vel(ind)^2/(2*(dx*xix(ind))*(dz*xiz(ind))^2);
%     L(ind) = -(omega^2)/(dx*xix(ind)) + vel(ind)^2/(2*(dx*xix(ind))*(dz*xiz(ind))^2);
%     N(ind) = 3i*omega*vel(ind)/(4*(dz*xiz(ind))^2) + vel(ind)^2/(4*(dx*xix(ind))*(dz*xiz(ind))^2);
%     S(ind) = 3i*omega*vel(ind)/(4*(dz*xiz(ind))^2) + vel(ind)^2/(4*(dx*xix(ind))*(dz*xiz(ind))^2);
%     NL(ind) = -vel(ind)^2/(4*(dx*xix(ind))*(dz*xiz(ind))^2);
%     SL(ind) = -vel(ind)^2/(4*(dx*xix(ind))*(dz*xiz(ind))^2);
% end
% 
% for n=1:length(LFT)
%     ind = LFT(n);
%     M(ind) = 1i*omega^3/vel(ind) + omega^2/(dx*xix(ind)) - 3i*omega*vel(ind)/(2*(dz*xiz(ind))^2) - vel(ind)^2/(2*(dx*xix(ind))*(dz*xiz(ind))^2);
%     R(ind) = -(omega^2)/(dx*xix(ind)) + vel(ind)^2/(2*(dx*xix(ind))*(dz*xiz(ind))^2);
%     N(ind) = 3i*omega*vel(ind)/(4*(dz*xiz(ind))^2) + vel(ind)^2/(4*(dx*xix(ind))*(dz*xiz(ind))^2);
%     S(ind) = 3i*omega*vel(ind)/(4*(dz*xiz(ind))^2) + vel(ind)^2/(4*(dx*xix(ind))*(dz*xiz(ind))^2);
%     NR(ind) = -vel(ind)^2/(4*(dx*xix(ind))*(dz*xiz(ind))^2);
%     SR(ind) = -vel(ind)^2/(4*(dx*xix(ind))*(dz*xiz(ind))^2);
% end
% 
% %% Corners
% NL([LFTTOP,RITTOP,LFTBOT,RITBOT]) = 0;
% L([LFTTOP,RITTOP,LFTBOT,RITBOT]) = 0;
% SL([LFTTOP,RITTOP,LFTBOT,RITBOT]) = 0;
% N([LFTTOP,RITTOP,LFTBOT,RITBOT]) = 0;
% S([LFTTOP,RITTOP,LFTBOT,RITBOT]) = 0;
% NR([LFTTOP,RITTOP,LFTBOT,RITBOT]) = 0;
% R([LFTTOP,RITTOP,LFTBOT,RITBOT]) = 0;
% SR([LFTTOP,RITTOP,LFTBOT,RITBOT]) = 0;
% 
% M(LFTTOP) = 1i*omega/vel(ind) + 1/sqrt((dz*xiz(ind))*(dx*xix(ind)));
% SR(LFTTOP) = -1/sqrt((dz*xiz(ind))*(dx*xix(ind)));
% 
% M(RITTOP) = 1i*omega/vel(ind) + 1/sqrt((dz*xiz(ind))*(dx*xix(ind)));
% SL(RITTOP) = -1/sqrt((dz*xiz(ind))*(dx*xix(ind)));
% 
% M(LFTBOT) = 1i*omega/vel(ind) + 1/sqrt((dz*xiz(ind))*(dx*xix(ind)));
% NR(LFTBOT) = -1/sqrt((dz*xiz(ind))*(dx*xix(ind)));
% 
% M(RITBOT) = 1i*omega/vel(ind) + 1/sqrt((dz*xiz(ind))*(dx*xix(ind)));
% NL(RITBOT) = -1/sqrt((dz*xiz(ind))*(dx*xix(ind)));

SR = [zeros(nzPML+1,1);SR(1:end-nzPML-1)];
R = [zeros(nzPML,1);R(1:end-nzPML)];
NR = [zeros(nzPML-1,1);NR(1:end-nzPML+1)];
S = [0;S(1:end-1)];

N = [N(2:end);0];
SL = [SL(nzPML:end);zeros(nzPML-1,1)];
L = [L(nzPML+1:end);zeros(nzPML,1)];
NL = [NL(nzPML+2:end);zeros(nzPML+1,1)];



% A = spdiags([NL,L,SL,N,M,S,NR,R,SR],[-nx-1:-nx+1,-1:1,nx-1:nx+1],NN,NN);
A = spdiags([NL,L,SL,N,M,S,NR,R,SR],[-nzPML-1:-nzPML+1,-1:1,nzPML-1:nzPML+1],NPML,NPML); %% Correct for non-square



