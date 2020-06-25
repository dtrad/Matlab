% In this script we remove the groundroll by taking the curvelet transform
% followed by a thresholding of the most important coefficients. The
% inverse curvelet transform of these most important coefficients is the
% predictor for the groundroll

% set the parameters for the curvelet transform

pfilt     =   'db8';
dfilt     =   'pkva';
dlevs     =   [1,1,2,2];

% make the mask

rlevs           =   cell(length(dlevs),1);
rlevs{1}{1}     =   1; % =1
rlevs{1}{2}     =   1; % =1 
rlevs{2}{1}     =   1; % =1
rlevs{2}{2}     =   1; % =1 
rlevs{3}{1}     =   1; % =1 
rlevs{3}{2}     =   0;
rlevs{4}{1}     =   1; % =1 
rlevs{4}{2}     =   1; % =1 
rlevs{4}{3}     =   0;
rlevs{4}{4}     =   0;
rlevs{5}{1}     =   1;
rlevs{5}{2}     =   1;
rlevs{5}{3}     =   0;
rlevs{5}{4}     =   0;


% test delete later
rlevs{3}{2}     =   0;
rlevs{4}{3}     =   0;
rlevs{4}{4}     =   0;
rlevs{5}{3}     =   0;
rlevs{5}{4}     =   0;


% run the script

[yd,y]    =   curvegroll('ozdata0win',0.011,pfilt,dfilt,dlevs,rlevs,0.001);
[yd,y]    =   curvegroll('ozdata1win',0.011,pfilt,dfilt,dlevs,rlevs,0.002);
[yd,y]    =   curvegroll('ozdata2win',0.011,pfilt,dfilt,dlevs,rlevs,0.001);
[yd,y]    =   curvegroll('ozdata3win',0.011,pfilt,dfilt,dlevs,rlevs,0.001);
[yd,y]    =   curvegroll('ozdata4win',0.011,pfilt,dfilt,dlevs,rlevs);
[yd,y]    =   curvegroll('ozdata5win',0.011,pfilt,dfilt,dlevs,rlevs);
[yd,y]    =   curvegroll('ozdata6win',0.011,pfilt,dfilt,dlevs,rlevs);
[yd,y]    =   curvegroll('ozdata7win',0.011,pfilt,dfilt,dlevs,rlevs);


% Save the curvelet transform image in a ps file for the paper
if (1)
  figure(2)
  print -dps /home/dtrad/ps/oz25.curvetransf.ps
end

return





