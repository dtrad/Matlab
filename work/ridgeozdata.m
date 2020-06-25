% needs to read binary data  obtained by sustrip in matlab
% These scripts produce files called ozdataNwin.noise.bin in
% /home/dtrad/work
% Look at the routine ridgegroll.m for details.


[yd,y,wy,twy]=ridgegroll('ozdata0win',0.011);

[yd,y,wy,twy]=ridgegroll('ozdata1win',0.021);
[yd,y,wy,twy]=ridgegroll('ozdata2win',0.021);
[yd,y,wy,twy]=ridgegroll('ozdata3win',0.041);
[yd,y,wy,twy]=ridgegroll('ozdata4win',0.111);
[yd,y,wy,twy]=ridgegroll('ozdata5win',0.111);
[yd,y,wy,twy]=ridgegroll('ozdata6win',0.111);
[yd,y,wy,twy]=ridgegroll('ozdata7win',0.111);



return