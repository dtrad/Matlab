function [cmap] = seiscolor;
% produces seismic color map
% Input:
%        ncolor - # of colors if empty ncolor = 32
%
% Output:
%        cmap   - seismic color map of length ncolor
%
% Usage:
%        [cmap] = seiscolor;

  if nargin <1,
    ncolor = 32;
  end
  tmp1 = [ones(1,ncolor) ; [0:(ncolor-1)]/(ncolor-1) ; [0:(ncolor-1)]/(ncolor-1)]';
  tmp2 = [[(ncolor-1):-1:0]/(ncolor-1) ; [(ncolor-1):-1:0]/31 ; ones(1,ncolor)]';
  cmap = [tmp1 ; tmp2];

% $$$ tmp1= [ ones(1,16); [0:15]/15; [0:15]/15 ];
% $$$ tmp2= [ [15:-1:0]/15; [15:-1:0]/15; ones(1,16) ];
% $$$ tmp3= [ zeros(1,16); [0:15]/15; [15:-1:0]/15 ];
% $$$ tmp4= [ [0:15]/15; [15:-1:0]/15; zeros(1,16) ];
% $$$ cmap = [tmp1 tmp2 tmp3 tmp4]';
