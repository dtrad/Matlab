!cd c:\dtrad\Matlab7\work
close all
clear all
A = load('zone_headersb.txt');
i = A(:,1);
cmp = A(:,2);
az = A(:,3);
cmpx = A(:,4);
cmpy = A(:,5);
offset = A(:,6);
shotx  = A(:,7);
shoty  = A(:,8);
rcvrx  = A(:,9);
rcvry  = A(:,10);
offsetx = A(:,11);
offsety = A(:,12);
midx = A(:,13);
midy = A(:,14);
alive = A(:,15);
figure;
!plot offset/azimuth for original
subplot(211);plot(offset(find(alive == 1)),az(find(alive == 1)),'.')

!plot offset/azimuth for interpolated
subplot(212);plot(offset(find(alive == 0)),az(find(alive == 0)),'.')

figure;
! plot midpoint for regular output
subplot(211);plot(midx(find(alive == 0)),midy(find(alive == 0)),'.')
subplot(212);plot(midx(find(alive == 1)),midy(find(alive == 1)),'.')

figure;hold on;
for i=100:100:1000
   
! plot offset and azimuth for one cmp (whatever cmp is in possition i)
    plot(offset(find(cmp == cmp(i))),az(find(cmp == cmp(i))),'.');
    pause
end
hold off;

clear A;
