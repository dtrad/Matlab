%function [seqno, spiralx, spiraly, spiral_seq, scol, srow, My_odd, Mx_odd]=spiral(Mx,My); 
function [spiralx, spiraly, My_odd, Mx_odd]=spiral_ver1(Mx,My); 
% calculates mapping from spiral COV to/from row,column
%    program first creates mapping from sequential number to column,row
%              spiral_seq, spiralx, spiraly
%    then figures out mapping from column,row to sequential number
%                    [scol' srow' seqno']
%       can access seqno directly by using index
%               index=(icol-1)*(2*Nmax+1)+irow
%
%
% Input
% Mx  fold in x-direction
% My  fold in y-direction
%        Nmax=5;  %maximum number of rows or columns from the center
%        flg_up:  true if 1st step of spiral is up, false if first step is
%        to the right

%
%  later:  need to account for rectangular COV bins
%

Mxc=ceil(Mx/2);  Mxr=abs(Mxc-Mx/2);
Myc=ceil(My/2);  Myr=abs(Myc-My/2);

%TODO: generalize so can handle different types of asymmetrical spreads
My_odd=0; Mx_odd=0;
if ((Mxr>0)&(Myr>0)) %both odd
    flg_up=0
    Mx_odd=1; My_odd=1;
    %col=-.5; row=-.5;
    col=-.5; row=.5;
elseif(Mxr>0)  %works for asymetrical spread where most of spread is to left of shot line
    flg_up=0
    Mx_odd=1;
    col=-.5; row=.5;
elseif(Myr>0)  %works for asymetrical spread where most of spread is above reciever line
    flg_up=1
    My_odd=1;
    col=-.5; row=.5;
else
    flg_up=1
    col=-.5; row=-.5;
end

%create mapping between spiral and row column
%Nmax=max(ceil(Mxc/2),ceil(Myc/2));
Nmax=max(Mxc,Myc);
m=1;
spiralx(m)=col; spiraly(m)=row;
Nk=0;

if (flg_up)
    for j=1:Nmax+1;
        Nk=Nk+1
        for k=1:Nk
            row=row+1;  % 1up
            m=m+1; spiralx(m)=col; spiraly(m)=row;
        end

        for k=1:Nk
            col=col+1;  % 1right
            m=m+1; spiralx(m)=col; spiraly(m)=row;
        end

        Nk=Nk+1
        for k=1:Nk
            row=row-1;  %2 down
            m=m+1; spiralx(m)=col; spiraly(m)=row;
        end

        for k=1:Nk
            col=col-1;  %2 left
            m=m+1; spiralx(m)=col; spiraly(m)=row;
        end
    end
else
    for j=1:Nmax+1;
        Nk=Nk+1;
        for k=1:Nk;
            col=col+1;  % 1right
            m=m+1; spiralx(m)=col; spiraly(m)=row;
        end

        for k=1:Nk;
            row=row-1;  % 1down
            m=m+1; spiralx(m)=col; spiraly(m)=row;
        end

        Nk=Nk+1;
        for k=1:Nk
            col=col-1;  %2 left
            m=m+1; spiralx(m)=col; spiraly(m)=row;
        end

        for k=1:Nk;
            row=row+1;  %2 up
            m=m+1; spiralx(m)=col; spiraly(m)=row;
        end
     end
end

%I=find((abs(spiralx) <=  Nmax) & (abs(spiraly)<=Nmax)); 
%spiralx=spiralx(I);
%spiraly=spiraly(I);
%N=length(spiraly);
%spiral_seq=[1:N];

figure(121); plot(spiralx,spiraly,'-*')

%[srow,I]=sort(spiraly);
%tmp=spiralx(I);
%tmp2=spiral_seq(I);
%[scol,J]=sort(tmp);
%srow=srow(J);
%seqno=tmp2(J);
%[scol' srow' seqno']

return

