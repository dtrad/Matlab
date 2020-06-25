load zone_int_brazeau
figure;
!plot offset/azimuth for original
subplot(211);plot(offset(find(alive == 1)),az(find(alive == 1)),'.')

!plot offset/azimuth for interpolated
subplot(212);plot(offset(find(alive == 0)),az(find(alive == 0)),'.')

figure;
! plot midpoint for regular output
subplot(211);plot(midx(find(alive == 0)),midy(find(alive == 0)),'.')
subplot(212);plot(midx(find(alive == 1)),midy(find(alive == 1)),'.')


for i=100:100:100
    figure(i);
    index = find(cmp == cmp(i));
    [x,ind]=binning(midy(index),60);plot(ind,x,'.')
    [x,ind]=binning(midx(index),30);plot(ind,x,'.')
    [x,ind]=binning(offset(index),10);plot(ind,x,'.')
    [x,ind]=binning(az(index),30);plot(ind,x,'.')
    ! plot offset and azimuth for one cmp (whatever cmp is in possition i)
    figure(i+1),
    subplot(111);plot(offset(index),az(index),'.');
    
    pause
end


clear A;
