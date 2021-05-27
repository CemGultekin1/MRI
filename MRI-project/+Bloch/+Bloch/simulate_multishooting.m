function [fingerprints,C,CTOT,grads]=simulate_multishooting(params,TR,TRFmax,sweep_phase,control,order,basis)
alphaseq=control(:,1);
TRFseq=control(:,2);
rfnum=length(alphaseq);
sweepnum=length(sweep_phase);
switch order
    case 0
        msize=5;
    case 1
        msize=40;
    case 2
        msize=120;
end
msize1=min(msize,40);
ERF=zeros(msize,5,rfnum);

pp=cell(rfnum,1);
for i=1:rfnum
   pp{i}=params; 
end
for k=1:rfnum
   ERF(:,:,k)=Bloch.intervalType1_onlyexponentiation(...
       alphaseq(k),TRFseq(k),pp{k},TR,order);
end

FP=Bloch.free_precession_exponentiation(params,TR,TRFmax,order);
EFP=zeros(min(msize,80),5,rfnum+2);
TRFseq=[0;TRFseq;0];
for k=1:rfnum+2
   A=zeros(size(FP,1));
   for i=1:size(FP,3)
       A=A+FP(:,:,i)*(TR/2-TRFseq(k)/2)^(i-1);
   end
   A=A(:,1:5);
   EFP(1:msize1,:,k)=A(1:msize1,:);
   if order==2
        M1=FP(:,:,2);
        EFP(msize1+1:end,:,k)=cat(1,M1*A*(-1/2));
   end
end

k1=msize1/5;

Rz=zeros(msize1,msize1,sweepnum);
for s=1:sweepnum
    phi = sweep_phase(s)+pi;
    Rz(:,:,s)=eye(msize1);
    for k=1:k1
        Rz((k-1)*5+(1:2),(k-1)*5+(1:2),s)=[cos(phi) -sin(phi);
                                        sin(phi)  cos(phi)];
    end
end

BB=Bloch.boundaryfun(params{:},TR,TRFmax);
BB0(1:msize1,1:5)=BB(1:msize1,1:5);
[MM,VV]=multishooting(rfnum,msize1,k1,EFP,ERF,BB0);
Y=repmat(VV,1,sweepnum);
for j=2:3:3*rfnum+2
    MM( msize1*(j-1)+(1:msize1),msize1*(j-1)+(1:msize1))=Rz(:,:,1)*MM( msize1*(j-1)+(1:msize1),msize1*(j-1)+(1:msize1));
end
Y(:,1)=MM\VV;
for s=2:sweepnum
    Rzs=Rz(:,:,s)*(Rz(:,:,s-1)^(-1));
    for j=2:3:3*rfnum+2
        MM( msize1*(j-1)+(1:msize1),msize1*(j-1)+(1:msize1))=Rzs*MM( msize1*(j-1)+(1:msize1),msize1*(j-1)+(1:msize1));
    end
    Y(:,s)=MM\VV;
end

Y=reshape(Y,msize1,(3*rfnum+3),sweepnum);
II=sort([(1:5:msize1),(2:5:msize1)]);
fingerprints=Y(II,2:3:end,:);
if order==0
    C=[];
    CTOT=[];
    grads=[];
    return;
end
%% Computing rel CRB
fng=fingerprints;
fng=reshape(fng,2,8,[]);% 2,8,rf+1,sweepnum
fng=permute(fng,[2,1,3]); % 8,2,rf+1,sweepnum
%phase_weights=reshape(phase_weights,1,[]);
%fng=reshape(fng,[],sweepnum).*phase_weights;
fng=reshape(fng,8,[]).'; % 2,rf+1,sweepnum,8
basis_flag=exist('basis','var');
if basis_flag
    basis_flag=~isempty(basis);
end
if basis_flag
    bfng=basis{1}*(basis{1}.'*fng);
    F=bfng.'*bfng;
else
    F=fng.'*fng;
end
I=eye(8);
%F=F+I;
FI=F\I(:,2:4);
C=diag(FI(2:4,:)); % Relative CRB values
CTOT=sum(C);
if order==1
    grads=[];
    return;
end


%% Computing dfdy
if basis_flag
    FF=FI*FI.';
    dfdy=-2*bfng*FF;
    QQ=basis{1};
    RR=basis{2};
    W=QQ*RR/(RR.'*RR);
    dfdu=dfdy*(fng.'*W)+fng*(dfdy.'*W)-QQ*((QQ.'*dfdy)*(fng.'*W)+(QQ.'*fng)*(dfdy.'*W));
    
else
    FF=FI*FI.';
    dfdy=-2*fng*FF;
end

dfdy=reshape(dfdy,2,(rfnum+1),sweepnum,8);
dfdy=permute(dfdy,[1,4,2,3]);
dfdy=reshape(dfdy,16,(rfnum+1),sweepnum);
 
%% Solving for adjoint
VP=zeros(40,3*rfnum+3,sweepnum);
VP(II,2:3:end,:)=dfdy;
VP=reshape(VP,[],sweepnum);
P=VP;
P(:,end)=MM.'\VP(:,end);
for s=1:sweepnum-1
    i= mod( s-2,sweepnum)+1;
    Rzs=Rz(:,:,s)*Rz(:,:,i)^(-1);
    for j=2:3:3*rfnum+2
        MM( msize1*(j-1)+(1:msize1),msize1*(j-1)+(1:msize1))=Rzs*MM( msize1*(j-1)+(1:msize1),msize1*(j-1)+(1:msize1));
    end
    P(:,s)=MM.'\VP(:,s);
end
P=reshape(P,msize1,(3*rfnum+3),sweepnum);
P=permute(P,[1,3,2]);
Y=permute(Y,[1,3,2]);
alphagrad=zeros(rfnum,1);
trfgrad=zeros(rfnum,1);

F=zeros(40);
for i=1:rfnum
    ii=3*i-1;
    F(:,1:5)=EFP(41:80,1:5,i+1); % dAdtrf
    for j=2:8
        F( (j-1)*5+(1:5),(j-1)*5+(1:5))=F(1:5,1:5);
    end
    pii=P(:,:,ii);
    yii=Y(:,:,ii);
    for s=1:sweepnum
        trfgrad(i)=trfgrad(i)-pii(:,s).'*Rz(:,:,s)*F*yii(:,s);
    end
    ii=3*i+1;
    pii=P(:,:,ii);
    yii=Y(:,:,ii);
    trfgrad(i)=trfgrad(i)-sum(diag(pii.'*F*yii));
    
    ii=3*i;
    F(:,1:5)=ERF(41:80,1:5,i); % dAdalpha
    for j=2:8
        F( (j-1)*5+(1:5),(j-1)*5+(1:5))=F(1:5,1:5);
    end
    pii=P(:,:,ii);
    yii=Y(:,:,ii);
    alphagrad(i)=alphagrad(i)-sum(diag(pii.'*F*yii));
    
    F(:,1:5)=ERF(81:120,1:5,i); % dAdtrf
    for j=2:8
        F( (j-1)*5+(1:5),(j-1)*5+(1:5))=F(1:5,1:5);
    end
    trfgrad(i)=trfgrad(i)-sum(diag(pii.'*F*yii));
    
end
grads=cat(2,alphagrad,trfgrad);
if basis_flag
   grads=cat(1,grads(:),dfdu(:)); 
end
end




function [MM,VV]=multishooting(rfnum,msize1,k1,EFP,ERF,BB0)
nmnz1=msize1*5*2-5*5;
%EFP=randn(msize,5,rfnum+2);
%ERF=randn(msize,5,rfnum);
innum=2*(rfnum+1)+rfnum;
bdnum=2;
trnum=innum+bdnum;
EE=zeros(nmnz1,trnum);

EE(:,1)=[reshape(EFP(1:msize1,:,1),[],1);reshape(repmat(EFP(1:5,1:5,1),1,k1-1),[],1)];
for i=1:rfnum
    EE(:,3*(i-1)+2)=[reshape(EFP(1:msize1,:,i+1),[],1);reshape(repmat(EFP(1:5,1:5,i+1),1,k1-1),[],1)];
    EE(:,3*(i-1)+3)=[reshape(ERF(1:msize1,:,i),[],1);reshape(repmat(ERF(1:5,1:5,i),1,k1-1),[],1)];
    EE(:,3*(i-1)+4)=EE(:,3*(i-1)+2);
end
EE(:,3*(rfnum-1)+5)=[reshape(EFP(1:msize1,:,end),[],1);reshape(repmat(EFP(1:5,1:5,end),1,k1-1),[],1)];

EE(:,end-1)=[BB0(:);reshape(repmat(BB0(1:5,1:5),1,k1-1),[],1)];
BB1=zeros(msize1,5);BB1(1:4,1:4)=-eye(4);
EE(:,end)=[BB1(:);reshape(repmat(BB1(1:5,1:5),1,k1-1),[],1)];
FF=-ones( msize1*innum,1);

SS=cat(1,EE(:),FF(:));

nmnz=trnum*nmnz1+innum*msize1;
II1=zeros(nmnz,1);
II2=zeros(nmnz,1);
[I2,I1]=meshgrid(1:5,1:msize1);
[D2,D1]=meshgrid(1:5,1:5);
D2=D2(:)+(1:msize1/5-1)*5;
D1=D1(:)+(1:msize1/5-1)*5;
I1=[I1(:);D1(:)];I2=[I2(:);D2(:)];
for i=1:innum
    II1( (i-1)*nmnz1+(1:nmnz1))=I1+(i-1)*msize1;
    II2( (i-1)*nmnz1+(1:nmnz1))=I2+(i-1)*msize1;
end
i=innum+1;
II1( (i-1)*nmnz1+(1:nmnz1))=I1+innum*msize1;
II2( (i-1)*nmnz1+(1:nmnz1))=I2;
i=innum+2;
II1( (i-1)*nmnz1+(1:nmnz1))=I1+innum*msize1;
II2( (i-1)*nmnz1+(1:nmnz1))=I2+innum*msize1;

II1(nmnz1*trnum+(1:innum*msize1))=(1:innum*msize1);
II2(nmnz1*trnum+(1:innum*msize1))=msize1+(1:innum*msize1);

MM=sparse(II1,II2,SS);
VV=zeros(size(MM,1),1);
VV(end-msize1+5)=1;
end

function [MM,VV]=multishooting1(rfnum,msize1,k1,EFP,ERF,BB0)
nmnz1=msize1*5*2-5*5;
%EFP=randn(msize,5,rfnum+2);
%ERF=randn(msize,5,rfnum);
EE=zeros(nmnz1,2*rfnum+1+2);

F1=EFP(1:msize1,:,1);F2=EFP(1:msize1,:,2);
FF=zeros(msize1,5);
FF(1:5,1:5)=F2(1:5,1:5)*F1(1:5,1:5);
for j=2:k1
FF( (j-1)*5+(1:5),1:5)=F2((j-1)*5+(1:5),1:5)*F1(1:5,1:5)+F2(1:5,1:5)*F1((j-1)*5+(1:5),1:5);
end
EE(:,1)=[FF(:);reshape(repmat(FF(1:5,1:5),1,k1-1),[],1)];
for i=1:rfnum
j=i+1;
EE(:,2*i)=[reshape(ERF(1:msize1,:,i),[],1);reshape(repmat(ERF(1:5,:,i),1,k1-1),[],1)];
F1=EFP(1:msize1,:,j);F2=EFP(1:msize1,:,j+1);
FF(1:5,1:5)=F2(1:5,1:5)*F1(1:5,1:5);
for j=2:msize1/5
FF( (j-1)*5+(1:5),1:5)=F2((j-1)*5+(1:5),1:5)*F1(1:5,1:5)+F2(1:5,1:5)*F1((j-1)*5+(1:5),1:5);
end
EE(:,2*i+1)=[FF(:);reshape(repmat(FF(1:5,1:5),1,k1-1),[],1)];
end
%BB0=randn(k1*5,5);
%BB1=randn(k1*5,5);
EE(:,end-1)=[BB0(:);reshape(repmat(BB0(1:5,1:5),1,k1-1),[],1)];
BB1=zeros(msize1,5);BB1(1:4,1:4)=-eye(4);
EE(:,end)=[BB1(:);reshape(repmat(BB1(1:5,1:5),1,k1-1),[],1)];
FF=-ones( msize1*(2*rfnum+1),1);

SS=cat(1,EE(:),FF(:));

nmnz=(2*rfnum+1+2)*(nmnz1)+(2*rfnum+1)*msize1;
II1=zeros(nmnz,1);
II2=zeros(nmnz,1);
[I2,I1]=meshgrid(1:5,1:msize1);
[D2,D1]=meshgrid(1:5,1:5);
D2=D2(:)+(1:msize1/5-1)*5;
D1=D1(:)+(1:msize1/5-1)*5;
I1=[I1(:);D1(:)];I2=[I2(:);D2(:)];
for i=1:2*rfnum+1
    II1( (i-1)*nmnz1+(1:nmnz1))=I1+(i-1)*msize1;
    II2( (i-1)*nmnz1+(1:nmnz1))=I2+(i-1)*msize1;
end
i=2*rfnum+1+1;
II1( (i-1)*nmnz1+(1:nmnz1))=I1+(i-1)*msize1;
II2( (i-1)*nmnz1+(1:nmnz1))=I2;
i=2*rfnum+1+2;
II1( (i-1)*nmnz1+(1:nmnz1))=I1+(i-2)*msize1;
II2( (i-1)*nmnz1+(1:nmnz1))=I2+(i-2)*msize1;

II1(nmnz1*(2*rfnum+3)+(1:(2*rfnum+1)*msize1))=(1:(2*rfnum+1)*msize1);
II2(nmnz1*(2*rfnum+3)+(1:(2*rfnum+1)*msize1))=msize1+(1:(2*rfnum+1)*msize1);

MM=sparse(II1,II2,SS);
VV=zeros(size(MM,1),1);
VV(end-msize1+5)=1;
end
