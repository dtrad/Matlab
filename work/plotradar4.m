orient landscape
dt=.9;
NP=170;NH=26;NT=640;
sumwigb('suradar3s','suradar3off',NT,NH,0,221,'(a)',dt);
sumwigb('suradar3rads','suradar3radoff',NT,NP,1,222,'(b)',dt);
sumwigb('suradar3recbs','suradar3off',NT,NH,0,223,'(c)',dt);
sumwigb('suradar3radbs','suradar3radoff',NT,NP,1,224,'(d)',dt);
subplot(222),xlabel('p (ns/m)')
subplot(224),xlabel('p (ns/m)')
subplot(221),ylabel('time (nsec)')
subplot(222),ylabel('time (nsec)')
subplot(223),ylabel('time (nsec)')
subplot(224),ylabel('time (nsec)')
