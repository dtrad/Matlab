orient landscape
dt=.9;
NP=170;
NH=26;
NT=640;
sumwigb('suradar3recbs','suradar3off',NT,NH,0,121,'(a)',dt);
sumwigbradar3('suradar3prads','suradar3pradoff',NT,NP,1,122,'(b)',dt);

subplot(121),xlabel('offset (m)')
subplot(122),xlabel('q  (ns/m)')
subplot(121),ylabel('time (nsec)')
subplot(122),ylabel('time (nsec)')

