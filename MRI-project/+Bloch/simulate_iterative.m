function [fingerprints,y01,p01,C,CTOT,grads]=simulate_iterative(params,TR,TRFmax,sweep_phase,control,order,y0guess,p0guess,basis)
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
yguessflag=0;
if exist('y0guess','var')
    yguessflag=~isempty(y0guess);
end
if ~yguessflag
    y0guess=[0;0;1-params{1};params{1};1;...
                        0;0;1;-1;0;...
                        zeros(5*6,1)];

    y0guess=y0guess(1:msize1);
    y0guess=repmat(y0guess,1,sweepnum);
else
    y0guess=y0guess(1:msize1,:);
end
y00=y0guess;
err=inf;
maxiter=20;
iter=0;
BB=Bloch.boundaryfun(params{:},TR,TRFmax);
BB0=zeros(msize1);
BB0(1:msize1,1:5)=BB(1:msize1,1:5);
for i=2:k1
    BB0( (i-1)*5+(1:5),(i-1)*5+(1:5))=BB0(1:5,1:5);
end

while err>1e-12 && iter<maxiter
    [y01,fingerprints,dAydalpha,dAydtrf]=forward_pass(rfnum,sweepnum,EFP,ERF,Rz,y00,k1,order,1,BB0);
    err=max(abs(y00-y01),[],[1,2]);
    iter=iter+1;
    y00=y01;
end
if err>1e-12
    [MM,VV]=multishooting(rfnum,msize1,k1,EFP,ERF,BB0);
    Y=repmat(VV,1,sweepnum);
    for s=1:sweepnum
        MM1=MM;
        for j=1:2:2*rfnum+1
            MM1( msize1*(j-1)+(1:msize1),msize1*(j-1)+(1:msize1))=Rz(:,:,s)*MM1( msize1*(j-1)+(1:msize1),msize1*(j-1)+(1:msize1));
        end
        Y(:,s)=MM1\VV;
    end
    Y=reshape(Y,msize1,(2*rfnum+2),sweepnum);
    II=sort([(1:5:msize1),(2:5:msize1)]);
    fingerprints=Y(II,2:end-1,:);
end
if order==0
    C=[];
    CTOT=[];
    grads=[];
    p01=[];
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
if rcond(F)<1e-15
    disp('here');
end
FI=F\I(:,2:4);
C=diag(FI(2:4,:)); % Relative CRB values
CTOT=sum(C);
if order==1
    grads=[];
    p01=[];
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

%dfdy=reshape(dfdy,2*(rfnum+1),sweepnum,8);
%dfdy=dfdy.*phase_weights;
dfdy=reshape(dfdy,2,(rfnum+1),sweepnum,8);
dfdy=permute(dfdy,[1,4,3,2]);
dfdy=reshape(dfdy,16,sweepnum,(rfnum+1));

%{
%% Solving adjoint the long way
dfdy__=reshape(dfdy,16,rfnum+1);
dfdy_=zeros(40,rfnum+3);
dfdy_(1:5:end,2:end-1)=dfdy__(1:2:end,:);
dfdy_(2:5:end,2:end-1)=dfdy__(2:2:end,:);
ptrue=MM.'\dfdy_(:);
ptrue=reshape(ptrue,40,[]);
%}
%% Solving for adjoint
I=(1:80);I(5:5:end)=[];J=(1:4);
EFP=EFP(I,J,:);
I=(1:120);I(5:5:end)=[];J=(1:4);
ERF=ERF(I,J,:);
I=(1:40);I(5:5:end)=[];
dAydalpha=dAydalpha(I,:,:);
dAydtrf=dAydtrf(I,:,:);
Rz=Rz(I,I,:);
%Qs=Qs(I,I,:);
BB0=BB0(I,I);
pguessflag=0;
if exist('p0guess','var')
    pguessflag=~isempty(p0guess);
end
if ~pguessflag
    p0guess=zeros(32,sweepnum);
else
    p0guess=p0guess(1:msize1,:);
end
p00=p0guess;
err=1;
iter=0;
while err>1e-12 && iter<maxiter
    [p01,alphagrad,trfgrad]=backward_pass(rfnum,sweepnum,EFP,ERF,Rz,dfdy,p00,dAydalpha,dAydtrf,BB0);
    err=max(abs(p00-p01),[],1)/sqrt(sum(p01.^2,1));
    iter=iter+1;
    p00=p01;
    if iter>9
        disp('here');
    end
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
EE=zeros(nmnz1,2*rfnum+1+2);

F1=EFP(1:msize1,:,1);F2=EFP(1:msize1,:,2);
FF=zeros(msize1,5);
FF(1:5,1:5)=F2(1:5,1:5)*F1(1:5,1:5);
for j=2:k1
FF( (j-1)*5+(1:5),1:5)=F2((j-1)*5+(1:5),1:5)*F1(1:5,1:5)+F2(1:5,1:5)*F1((j-1)*5+(1:5),1:5);
end
EE(:,1)=[FF(:);reshape(repmat(FF(1:5,1:5),k1-1,1),[],1)];
for i=1:rfnum
j=i+1;
EE(:,2*i)=[reshape(ERF(1:msize1,:,i),[],1);reshape(repmat(ERF(1:5,:,i),k1-1,1),[],1)];
F1=EFP(1:msize1,:,j);F2=EFP(1:msize1,:,j+1);
FF(1:5,1:5)=F2(1:5,1:5)*F1(1:5,1:5);
for j=2:msize1/5
FF( (j-1)*5+(1:5),1:5)=F2((j-1)*5+(1:5),1:5)*F1(1:5,1:5)+F2(1:5,1:5)*F1((j-1)*5+(1:5),1:5);
end
EE(:,2*i+1)=[FF(:);reshape(repmat(FF(1:5,1:5),k1-1,1),[],1)];
end
%BB0=randn(k1*5,5);
%BB1=randn(k1*5,5);
EE(:,end-1)=[BB0(:);reshape(repmat(BB0(1:5,1:5),k1-1,1),[],1)];
BB1=zeros(msize1,5);BB1(1:4,1:4)=-eye(4);
EE(:,end)=[BB1(:);reshape(repmat(BB1(1:5,1:5),k1-1,1),[],1)];
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
function [y0,fingerprints,dAydalpha,dAydtrf]=forward_pass(rfnum,sweepnum,EFP,ERF,Rz,y0,k1,order,collect,BB0)
dAydtrf=[];
dAydalpha=[];
if collect && order==0
    fingerprints=zeros(2,rfnum+1,sweepnum);
elseif collect
    fingerprints=zeros(16,rfnum+1,sweepnum);
    if order==2
        dAydtrf=zeros(40,sweepnum,rfnum*3);
        dAydalpha=zeros(40,sweepnum,rfnum);
    end
end
y=y0;
%A=EFP(1:msize1,1:5,1);
for i=2:k1
  y( (i-1)*5+(1:5),:)=EFP((i-1)*5+(1:5),1:5,1)*y( 1:5,:)+EFP(1:5,1:5,1)*y( (i-1)*5+(1:5),:);
end
y(1:5,:)=EFP(1:5,1:5,1)*y(1:5,:);
if collect
    fingerprints(1:2:end,1,:)=y(1:5:end,:);
    fingerprints(2:2:end,1,:)=y(2:5:end,:);
end
for i=1:rfnum
  jj=i+1;
  if order==2 && collect
    for j=2:k1
      dAydtrf( (j-1)*5+(1:5),:,(i-1)*3+1)=EFP(40+(j-1)*5+(1:5),1:5,jj)*y( 1:5,:)+EFP(41:45,1:5,jj)*y( (j-1)*5+(1:5),:);
    end
    dAydtrf( 1:5,:,(i-1)*3+1)=EFP(41:45,1:5,jj)*y(1:5,:);
    for s=1:sweepnum
        dAydtrf(:,s,(i-1)*3+1)=Rz(:,:,s)*dAydtrf(:,s,(i-1)*3+1);
    end
  end
  for j=2:k1
      y( (j-1)*5+(1:5),:)=EFP((j-1)*5+(1:5),1:5,jj)*y( 1:5,:)+EFP(1:5,1:5,jj)*y( (j-1)*5+(1:5),:);
  end
  y(1:5,:)=EFP(1:5,1:5,jj)*y(1:5,:);
  for s=1:sweepnum
      y(:,s)=Rz(:,:,s)*y(:,s);
  end
  if order==2 && collect
        for j=2:k1
          dAydtrf( (j-1)*5+(1:5),:,(i-1)*3+2)=ERF(80+(j-1)*5+(1:5),1:5,i)*y( 1:5,:)+ERF(81:85,1:5,i)*y( (j-1)*5+(1:5),:);
        end
        dAydtrf( 1:5,:,(i-1)*3+2)=ERF(81:85,1:5,i)*y(1:5,:);
        for j=2:k1
          dAydalpha( (j-1)*5+(1:5),:,i)=ERF(40+(j-1)*5+(1:5),1:5,i)*y( 1:5,:)+ERF(41:45,1:5,i)*y( (j-1)*5+(1:5),:);
        end
        dAydalpha( 1:5,:,i)=ERF(41:45,1:5,i)*y(1:5,:);
  end
  for j=2:k1
      y((j-1)*5+(1:5),:)=ERF((j-1)*5+(1:5),1:5,i)*y(1:5,:)+ERF(1:5,1:5,i)*y((j-1)*5+(1:5),:);
  end
  y(1:5,:)=ERF(1:5,1:5,i)*y(1:5,:);
  
  
  if order==2 && collect
    for j=2:k1
      dAydtrf( (j-1)*5+(1:5),:,(i-1)*3+3)=EFP(40+(j-1)*5+(1:5),1:5,jj)*y( 1:5,:)+EFP(41:45,1:5,jj)*y( (j-1)*5+(1:5),:);
    end
    dAydtrf( 1:5,:,(i-1)*3+3)=EFP(41:45,1:5,jj)*y(1:5,:);
  end
  
  for j=2:k1
      y( (j-1)*5+(1:5),:)=EFP((j-1)*5+(1:5),1:5,jj)*y( 1:5,:)+EFP(1:5,1:5,jj)*y( (j-1)*5+(1:5),:);
  end
  y(1:5,:)=EFP(1:5,1:5,jj)*y(1:5,:);
  if collect
      fingerprints(1:2:end,jj,:)=y(1:5:end,:);
      fingerprints(2:2:end,jj,:)=y(2:5:end,:);
  end
end

jj=size(EFP,3);
for j=2:k1
  y( (j-1)*5+(1:5),:)=EFP((j-1)*5+(1:5),1:5,jj)*y( 1:5,:)+EFP(1:5,1:5,jj)*y( (j-1)*5+(1:5),:);
end
y(1:5,:)=EFP(1:5,1:5,jj)*y(1:5,:);

for s=1:sweepnum
   y(:,s)=Rz(:,:,s)*y(:,s);
end
y0=BB0\y;
end


function [p0,alphagrad,trfgrad]=backward_pass(rfnum,sweepnum,EFP,ERF,Rz,dfdy,p0,dAydalpha,dAydtrf,BB0)
trfgrad=zeros(rfnum,1);
alphagrad=zeros(rfnum,1);
p=-p0;
%{
p_n+1
for k=1:8
  p( (k-1)*4+(1:4),:)=-p( (k-1)*4+(1:4),:);
end
%}
% p_n
jj=rfnum+2;
for s=1:sweepnum
    p(:,s)=Rz(:,:,s).'*p(:,s);
end
p(1:4,:)=EFP(1:32,:,jj).'*p;
for k=2:8
    p( (k-1)*4+(1:4),:)=EFP(1:4,1:4,jj).'*p( (k-1)*4+(1:4),:);
end

for i=rfnum:-1:1
  
  jj=i+1;
  p(1:4:end,:)=p(1:4:end,:)-dfdy(1:2:end,:,jj);
  p(2:4:end,:)=p(2:4:end,:)-dfdy(2:2:end,:,jj);
  
  for s=1:sweepnum
      trfgrad(i)=trfgrad(i)-p(:,s).'*dAydtrf(:,s,(i-1)*3+3);
  end
  % p_n-1
  p(1:4,:)=EFP(1:32,:,jj).'*p;
  for k=2:8
      p( (k-1)*4+(1:4),:)=EFP(1:4,1:4,jj).'*p( (k-1)*4+(1:4),:);
  end
  
  for s=1:sweepnum
      trfgrad(i)=trfgrad(i)-p(:,s).'*dAydtrf(:,s,(i-1)*3+2);
      alphagrad(i)=alphagrad(i)-p(:,s).'*dAydalpha(:,s,i);
  end
  p(1:4,:)=ERF(1:32,:,i).'*p;
  for k=2:8
      p( (k-1)*4+(1:4),:)=ERF(1:4,1:4,i).'*p( (k-1)*4+(1:4),:);
  end
  
  for s=1:sweepnum
      trfgrad(i)=trfgrad(i)-p(:,s).'*dAydtrf(:,s,(i-1)*3+1);
  end
  for s=1:sweepnum
      p(:,s)=Rz(:,:,s).'*p(:,s);
  end
  
  p(1:4,:)=EFP(1:32,:,jj).'*p;
  for k=2:8
      p( (k-1)*4+(1:4),:)=EFP(1:4,1:4,jj).'*p( (k-1)*4+(1:4),:);
  end
end

jj=1;
p(1:4:end,:)=p(1:4:end,:)-dfdy(1:2:end,:,jj);
p(2:4:end,:)=p(2:4:end,:)-dfdy(2:2:end,:,jj);
p(1:4,:)=EFP(1:32,:,jj).'*p;
for k=2:8
  p( (k-1)*4+(1:4),:)=EFP(1:4,1:4,jj).'*p( (k-1)*4+(1:4),:);
end
p0=-BB0.'\p;
end