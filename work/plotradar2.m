orient landscape
dt=0.03;
NP=170;NH=100;
sumwigb('suradar2s','suradar2boff',501,NH,2,221,'(a)',dt);
sumwigb('suradar2radhrs','radar2radoff',501,NP,1,222,'(b)',dt);
sumwigb('suradar2recbhrs','suradar2boff',501,NH,2,223,'(c)',dt);
sumwigb('suradar2radbhrs','radar2radoff',501,NP,1,224,'(d)',dt);
subplot(222),xlabel('slope')
subplot(224),xlabel('slope')
