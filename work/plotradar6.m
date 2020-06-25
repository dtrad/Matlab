figure
orient landscape
dt=.9;
NP=200;
NH=26;
NT=640;
sumwigb('suradar3phrrecbs','suradar3off',NT,NH,0,122,'(b)',dt);
sumwigbradar3('suradar3phrrads','suradar3ptradoff',NT,NP,1,121,'(a)',dt);

subplot(122),xlabel('offset (m)')
subplot(121),xlabel('q  (ns/m^2)')
subplot(122),ylabel('time (nsec)')
subplot(121),ylabel('time (nsec)')

