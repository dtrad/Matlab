function [ D ] = Get_data( FDFD, freq, fwave, vel_true, R )
            
U = FDFD(freq,fwave,vel_true);
D=zeros(size(R,1),size(U,2),size(U,3));
if length(freq)>1
    for jjj=1:length(freq)
        D(:,:,jjj) = R*U(:,:,jjj);
    end
else
    D = R*U;
end

D=full(D);

end

