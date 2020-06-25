close all
clear net
clear all
% P = [0 1 2 3 4 5 6 7 8 9 10];
% T = [0 1 2 3 4 3 2 1 2 3 4];

nn2dto3d %load data

n1=50;
n2=30;
n10=261;

P=P(n10:n10+n1-1,1:n2); % 2D prediction as input
T=T(n10:n10+n1-1,1:n2); % 3D full prediction is the target

% Using 3D sparse prediction instead
%P=zeros(size(T));
%P(:,1:4:end)=T(:,1:4:end);

figure(1);imagesc(P);title('Input (2D or 3D sparse)');
figure(2);imagesc(T);title('Target (3D full prediction)');



if (0)
  P=P(:).';
  T=T(:).';
end

 net = newff(minmax(P),[15 50],{'tansig' 'purelin'});
 
 Y = sim(net,P);
 figure(1);
 
 net.trainParam.epochs = 10;
 net = train(net,P,T);
 Y = sim(net,P);
 

if (0)  
 P=reshape(P(:),n1,n2);
 T=reshape(T(:),n1,n2);
 Y=reshape(Y(:),n1,n2);
end

figure(3);imagesc(Y);title('3D prediction');


% Take a window not in the training set

nn2dto3d %load data
n1=50;
n2=30;
n10=261;

P=P(n10:n10+n1-1,n2+1:2*n2); % 2D prediction as input

T=T(n10:n10+n1-1,n2+1:2*n2); % 3D full prediction is the target
% Using 3D sparse prediction instead
%P=zeros(size(T));
%P(:,1:4:end)=T(:,1:4:end);

figure(4);imagesc(P);title('Input (2D or 3D sparse)');
figure(5);imagesc(T);title('Target (3D full prediction)');

 Y = sim(net,P);
 figure(6);imagesc(Y);title('3D prediction');