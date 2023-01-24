function [ S, R ] = Define_Acquisition( sz, sx, rz, rx, nz, nx, PML_thick )

nzPML = nz+2*PML_thick;
nxPML = nx+2*PML_thick;
NPML = nzPML*nxPML;

sind = ((sx+PML_thick-1)*(nzPML)+PML_thick+sz);
S = sparse(1:length(sind),sind,ones(size(sind)),length(sind),NPML);
S=S';

PlotS = reshape(S*S'*ones(size(S,1),1),nzPML,nxPML);
imagescf(PlotS(PML_thick + 1: PML_thick + nz, PML_thick + 1: PML_thick +nx))
colorbar;
title('Sources')
drawnow

rind = ((rx+PML_thick-1)*(nzPML)+PML_thick+rz);
R = sparse(1:length(rind),rind,ones(size(rind)),length(rind),NPML);

PlotR = reshape(R'*R*ones(size(R,2),1),nzPML,nxPML);
imagescf(PlotR(PML_thick + 1: PML_thick + nz, PML_thick + 1: PML_thick +nx))
colorbar;
title('Receivers')
drawnow

end

