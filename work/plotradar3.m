orient landscape
dt=0.03;
NP=170;NH=100;
NT=1024;
sumwigbradar('suradar2cs','suradar2boff',NT,NH,2,221,'(a)',dt);
sumwigbradar2('suradar2crads','radar2radoff',NT,NP,1,222,'(b)',dt);
sumwigbradar('suradar2crecbs','suradar2boff',NT,NH,2,223,'(c)',dt);
sumwigbradar2('suradar2cradbs','radar2radoff',NT,NP,1,224,'(d)',dt);
subplot(222),xlabel('slope')
subplot(224),xlabel('slope')
