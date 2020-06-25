cd V:\prg\dtrad\projects\narrasyn
load ampafter.dat
cdp=356153;
ampint=ampafter(find(ampafter(:,1)==cdp),:);
clear ampafter;

printflag=1

[n m]=size(ampint);
orig=zeros(size(ampint));
new=[];
ori=[];
for i=1:n
    if (ampint(i,5)==0)
        new=[new,i];
    else
        ori=[ori,i];
    end
end

ampnew=ampint(new,:);
ampori=ampint(ori,:);

avo_ori = polyfit(ampori(:,2),ampori(:,4),1);
avo_int = polyfit(ampint(:,2),ampint(:,4),1);
avo_new = polyfit(ampnew(:,2),ampnew(:,4),1);

line_ori = avo_ori(1)*ampori(:,2)+avo_ori(2);
line_int = avo_int(1)*ampint(:,2)+avo_int(2);
line_new = avo_new(1)*ampnew(:,2)+avo_new(2);

size(new)
size(ori)

v=[0 6000 3000 9000];
figure(1)
subplot(111);plot(ampori(:,2),ampori(:,4),'.',ampori(:,2),line_ori);title('original','FontSize',16);
ylabel('amp','FontSize',14);xlabel('offset','FontSize',14);axis(v);
text(1000,8000,['slope=' num2str(avo_ori(1))])
text(1000,7500,['intercept=' num2str(round(avo_ori(2)))])

figure(2)
subplot(111);plot(ampint(:,2),ampint(:,4),'.',ampint(:,2),line_int);title('interpolated','FontSize',16);
ylabel('amp','FontSize',14);xlabel('offset','FontSize',14);axis(v);
text(1000,8000,['slope=' num2str(avo_int(1))])
text(1000,7500,['intercept=' num2str(round(avo_int(2)))])

figure(3),
subplot(111);plot(ampnew(:,2),ampnew(:,4),'.',ampnew(:,2),line_new);title('only interpolated','FontSize',16);
xlabel('offset','FontSize',14);ylabel('amp','FontSize',14);axis(v);
text(1000,8000,['slope=' num2str(avo_new(1))])
text(1000,7500,['intercept=' num2str(round(avo_new(2)))])

if (printflag)
    figure(1);print -dpng avoorig.png
    figure(2);print -dpng avointerp.png
    figure(3);print -dpng avointeronly.png
end

%if azimuth variation use this
if (0==1)
    subplot(223);plot(ampori(:,3),ampori(:,4),'.');title('original');
    ylabel('amp');
    subplot(224);plot(ampint(:,3),ampint(:,4),'.');title('interpolated');
    xlabel('azimuth');ylabel('amp');
end
