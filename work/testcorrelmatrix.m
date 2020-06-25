% Test correlmatrix

if (1) 
  load data2m.dat;
  data2m=data2m.';
  load wav2d.dat
  wav2d=wav2d.';
end,


for sx=600:600:3000
  for st=0.005:0.005:0.005
    dataout=correlmatrix(data2m,wav2d);

    datascale=datascale+correlmatrix(dataout,wav2d);
  end,
end,

    
figure(1)
wigb(dataout);
