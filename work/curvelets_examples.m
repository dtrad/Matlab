 fid=fopen('sss.raw','wb');
 fwrite(fid,sss,'char');
 fclose(fid)

 % copy sss.raw to the BeamLab200\.... directory in place of BigMac
PrepareThreshCurvelet('c:\daniel\datasets','BigMac','ReadRaw1(''BigMac'');')
DoUnifThreshCurvelet('c:\daniel\datasets','Thres','BigMac', 2000,40,0,0);
ViewThreshCurvelet('c:\daniel\datasets','Thres','BigMac',0);
ViewThreshCurvelet2('c:\daniel\datasets','Thres','BigMac',0);
