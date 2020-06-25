figure
orient landscape
dt=.9;
NP=100;
NH=26;
NT=640;
sumwigb('suradar3recs','suradar3off',NT,NH,0,121,'(a)',dt);
sumwigb('suradar3phrrecbs','suradar3off',NT,NH,0,122,'(b)',dt);

subplot(121),ylabel('time (nsec)')
subplot(122),ylabel('time (nsec)')

