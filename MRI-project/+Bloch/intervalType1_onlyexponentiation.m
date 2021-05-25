function EC=intervalType1_onlyexponentiation(alpha,TRF,params,TR,order)
M=Bloch.Mfun(alpha,TRF,params{:},TR,1)*TRF;
n=5;
if order==0
    EC=expm(M(1:5,1:5));
elseif order==1
    C=zeros(8*n);
    C(:,1:5)=M(:,1:5);
    for i=2:8
        C( (i-1)*5+(1:5),(i-1)*5+(1:5))=M(1:5,1:5);
    end
    EC=expm(C); 
    EC=EC(:,1:5);
elseif order==2
    Malpha=Bloch.dMdalphafun(alpha,TRF,params{:},TR,1)*TRF;
    MTRF=Bloch.dMdTRFfun(alpha,TRF,params{:},TR,1)*TRF^2;
    EC=zeros(8*5*3);
    EC(:,1:5)=cat(1,M,Malpha,MTRF);
    EC(8*5+(1:8*5),8*5+(1:5))=M;
    EC(8*5*2+(1:8*5),8*5*2+(1:5))=M;
    A=M(1:5,1:5);
    for i=1:8*3
        EC( (i-1)*5+(1:5),(i-1)*5+(1:5))=A;
    end
    A=Malpha(1:5,1:5);
    for i=1:8
        EC( 8*5+(i-1)*5+(1:5),(i-1)*5+(1:5))=A;
    end
    
    A=MTRF(1:5,1:5);
    for i=1:8
        EC( 8*5*2+(i-1)*5+(1:5),(i-1)*5+(1:5))=A;
    end
    M1=EC(1:40,1:40);
    EC=expm(EC);
    EC=EC(:,1:5);
    EC(81:120,:)=EC(81:120,:)+M1*EC(1:40,:);
    EC(81:120,:)=EC(81:120,:)/TRF;
end
end