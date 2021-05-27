function MM=free_precession_exponentiation(params,TR,TRFmax,order)
M=Bloch.Mfun(0,1e-3,params{:},TR,0);
if order==0
    M=M(1:5,1:5);
end
k=2;
m=size(M,1);
MM=zeros(m,m,2);
A=zeros(m);
A(:,1:5)=M;
i1=m/5;
for j=2:i1
    A((j-1)*5+(1:5),(j-1)*5+(1:5))=M(1:5,1:5);
end
MM(:,:,1)=eye(m);
EA=abs(A)*TRFmax;
ERR=EA;
while max(ERR(:))>1e-19
    if k>size(MM,3)
        MM=cat(3,MM,MM);
    end
    MM(:,:,k)=A*MM(:,:,k-1)/(k-1);
    ERR=EA*ERR/(k-1);
    k=k+1;
end
MM=MM(:,:,1:k-1);

end