%{
clear;
clc;
m0s=.1;
T1=1.6;
T2f=0.065;
R=30;
T2s=60e-6;
B0=0;
B1=1;
TR=4.5e-3;
alpha=0.1;
TRF=1e-3;
TRFmax=1e-3;
load('+BlochSim/defaultControl.mat');
rfnum=666;
alphaseq= interp1q(linspace(0,1,100).',defaultControl(:,1),linspace(0,1,rfnum).')*2;
TRFseq= interp1q(linspace(0,1,100).',defaultControl(:,2),linspace(0,1,rfnum).');
control=[alphaseq,TRFseq];
sweep_phase=linspace(-pi/2,pi/2,266);
params={m0s,T1,T2f,R,T2s,B0,B1};
order=2;
%}

function [fingerprints,C,CTOT,grads]=simulate(params,TR,TRFmax,sweep_phase,control,order,basis,extras)
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
A=repmat(EFP(1:msize1,:,1),1,sweepnum);

Rz=zeros(msize1,msize1,sweepnum);
for s=1:sweepnum
    phi = sweep_phase(s)+pi;
    Rz(:,:,s)=eye(msize1);
    for k=1:k1
        Rz((k-1)*5+(1:2),(k-1)*5+(1:2),s)=[cos(phi) -sin(phi);
                                        sin(phi)  cos(phi)];
    end
end

for i=1:rfnum
  jj=i+1;
  D1=A(1:5,:);
  D2=EFP(1:5,1:5,jj);
  A(1:5,:)=D2*D1;
  for j=2:k1
      A((j-1)*5+(1:5),:)=D2*A((j-1)*5+(1:5),:)+EFP((j-1)*5+(1:5),1:5,jj)*D1;
  end
  for s=1:sweepnum
    A(:,(s-1)*5+(1:5))=Rz(:,:,s)*A(:,(s-1)*5+(1:5));
  end
  D1=A(1:5,:);
  D2=ERF(1:5,1:5,i);
  A(1:5,:)=D2*D1;
  for j=2:k1
      A((j-1)*5+(1:5),:)=D2*A((j-1)*5+(1:5),:)+ERF((j-1)*5+(1:5),1:5,i)*D1;
  end
  D1=A(1:5,:);
  D2=EFP(1:5,1:5,jj);
  A(1:5,:)=D2*D1;
  for j=2:k1
      A((j-1)*5+(1:5),:)=D2*A((j-1)*5+(1:5),:)+EFP((j-1)*5+(1:5),1:5,jj)*D1;
  end
end
jj=size(EFP,3);
D1=A(1:5,:);
D2=EFP(1:5,1:5,jj);
A(1:5,:)=D2*D1;
for j=2:k1
  A((j-1)*5+(1:5),:)=D2*A((j-1)*5+(1:5),:)+EFP((j-1)*5+(1:5),1:5,jj)*D1;
end
for s=1:sweepnum
    A(:,(s-1)*5+(1:5))=Rz(:,:,s)*A(:,(s-1)*5+(1:5));
end
%A(1:5:end,:)=-A(1:5:end,:);
%A(2:5:end,:)=-A(2:5:end,:);

BB=Bloch.boundaryfun(params{:},TR,TRFmax);
BB0=BB(1:msize1,1:5);
BB1=BB(41:45,1:5);
%BB2=zeros(msize1,1);
%BB2(1:5)=BB(46,:);

y0s=zeros(msize1,sweepnum);
Qs=zeros(msize1,msize1,sweepnum);
for s=1:sweepnum
    MA=A(:,(s-1)*5+(1:5));
    D1=MA(1:5,:);
    D2=BB1(1:5,1:5);
    MA(1:5,:)=D2*D1;
    for j=2:k1
      MA((j-1)*5+(1:5),:)=D2*MA((j-1)*5+(1:5),:);
    end
    Q=[MA+BB0,zeros(msize1,(k1-1)*5)];
    for j=2:k1
        Q( (j-1)*5+(1:5),(j-1)*5+(1:5))=Q(1:5,1:5);
    end
    Qs(:,:,s)=Q;
    QA=zeros(msize1/5*4);
    bb=zeros(msize1/5*4,1);
    for i=1:msize1/5
        QA( (i-1)*4+(1:4),(i-1)*4+(1:4))=Q((i-1)*5+(1:4),(i-1)*5+(1:4));
        bb((i-1)*4+(1:4))=Q((i-1)*5+(1:4),5);
    end
    %rcond(QA(1:4,1:4))
    if rcond(QA(1:4,1:4))<1e-6
        disp(['B0 value: ', num2str(params{end-1})]);
    end
    yy0=-QA\bb;
    for i=1:msize1/5
        y0s((i-1)*5+(1:4),s)=yy0((i-1)*4+(1:4));
    end
    y0s(5,s)=1;
    %y0s_=Q\BB2;
end
%{
%% Full matrix build
MM=zeros( (rfnum+2+1)*40);
ii=1;
EE=EFP(1:40,1:5,1);
MM( (ii-1)*40+(1:40),(ii-1)*40+(1:5))=EE;
for j=2:8
    MM((ii-1)*40+(j-1)*5+(1:5),(ii-1)*40+(j-1)*5+(1:5))=EE(1:5,1:5);
end
MM( (ii-1)*40+(1:40),(ii)*40+(1:40))=-eye(40);
for i=1:rfnum
    ii=i+1;
    MM( (ii-1)*40+(1:40),(ii)*40+(1:40))=-eye(40);
    EE=EFP(1:40,1:5,i+1);
    UU=zeros(40);
    UU(:,1:5)=EE;
    for j=2:8
        UU((j-1)*5+(1:5),(j-1)*5+(1:5))=EE(1:5,1:5);
    end
    UM=UU;
    EE=ERF(1:40,1:5,i);
    EE(:,1:2)=-EE(:,1:2);
    UU=zeros(40);
    UU(:,1:5)=EE;
    for j=2:8
        UU((j-1)*5+(1:5),(j-1)*5+(1:5))=EE(1:5,1:5);
    end
    UM=UU*UM;
    EE=EFP(1:40,1:5,i+1);
    UU=zeros(40);
    UU(:,1:5)=EE;
    for j=2:8
        UU((j-1)*5+(1:5),(j-1)*5+(1:5))=EE(1:5,1:5);
    end
    UM=UU*UM;
    MM( (ii-1)*40+(1:40),(ii-1)*40+(1:40))=UM;
end
ii=rfnum+2;
EE=EFP(1:40,1:5,end);
EE(:,1:2)=-EE(:,1:2);
MM( (ii-1)*40+(1:40),(ii-1)*40+(1:5))=EE;
for j=2:8
    MM((ii-1)*40+(j-1)*5+(1:5),(ii-1)*40+(j-1)*5+(1:5))=EE(1:5,1:5);
end
MM( (ii-1)*40+(1:40),(ii)*40+(1:40))=-eye(40);

MM(end-40+1:end,1:5)=BB0;
for j=2:8
    MM(end-40+(j-1)*5+(1:5),(j-1)*5+(1:5))=BB0(1:5,1:5);
end
for j=1:8
    MM(end-40+(j-1)*5+(1:5),end-40+(j-1)*5+(1:5))=BB1(1:5,1:5);
end
VV=zeros(size(MM,1),1);
VV(end-40+1:end)=BB2;
YY=reshape(MM\VV,40,[]);
ftrue=YY(1:5:end,:);
ftrue=[ftrue(:,2:end-1)];
%}
%% Forward Pass
if order==0
    fingerprints=zeros(2,rfnum+1,sweepnum);
else
    fingerprints=zeros(16,rfnum+1,sweepnum);
    if order==2
        dAydtrf=zeros(40,sweepnum,rfnum*3);
        dAydalpha=zeros(40,sweepnum,rfnum);
    end
end
y=y0s;
%A=EFP(1:msize1,1:5,1);
for i=2:k1
  y( (i-1)*5+(1:5),:)=EFP((i-1)*5+(1:5),1:5,1)*y( 1:5,:)+EFP(1:5,1:5,1)*y( (i-1)*5+(1:5),:);
end
y(1:5,:)=EFP(1:5,1:5,1)*y(1:5,:);

fingerprints(1:2:end,1,:)=y(1:5:end,:);
fingerprints(2:2:end,1,:)=y(2:5:end,:);
for i=1:rfnum
  jj=i+1;
  if order==2
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
  if order==2
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
  
  
  if order==2
    for j=2:k1
      dAydtrf( (j-1)*5+(1:5),:,(i-1)*3+3)=EFP(40+(j-1)*5+(1:5),1:5,jj)*y( 1:5,:)+EFP(41:45,1:5,jj)*y( (j-1)*5+(1:5),:);
    end
    dAydtrf( 1:5,:,(i-1)*3+3)=EFP(41:45,1:5,jj)*y(1:5,:);
  end
  
  for j=2:k1
      y( (j-1)*5+(1:5),:)=EFP((j-1)*5+(1:5),1:5,jj)*y( 1:5,:)+EFP(1:5,1:5,jj)*y( (j-1)*5+(1:5),:);
  end
  y(1:5,:)=EFP(1:5,1:5,jj)*y(1:5,:);
  fingerprints(1:2:end,jj,:)=y(1:5:end,:);
  fingerprints(2:2:end,jj,:)=y(2:5:end,:);
end
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

if contains(extras,'T1pT2f')
    FI=F\I(:,3:4);
    C=diag(FI(3:4,:));
elseif contains(extras, 'R-only')
    FI=F\I(:,5);
    C=diag(FI(5,:));
elseif contains(extras, 'default')
    FI=F\I(:,2:4);
    C=diag(FI(2:4,:));
end


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
Qs=Qs(I,I,:);
BB1=BB1(1:4,1:4);

pz=zeros(32,sweepnum);
%pz_track=zeros(40,sweepnum,rfnum);
for i=rfnum:-1:1
  jj=i+1;
  pz(1:4:end,:)=pz(1:4:end,:)-dfdy(1:2:end,:,jj);
  pz(2:4:end,:)=pz(2:4:end,:)-dfdy(2:2:end,:,jj);
  pz(1:4,:)=EFP(1:32,:,jj).'*pz;
  for k=2:8
      pz( (k-1)*4+(1:4),:)=EFP(1:4,1:4,jj).'*pz( (k-1)*4+(1:4),:);
  end
  
  
 pz(1:4,:)=ERF(1:32,:,i).'*pz;
  for k=2:8
      pz( (k-1)*4+(1:4),:)=ERF(1:4,1:4,i).'*pz( (k-1)*4+(1:4),:);
  end
  
  for s=1:sweepnum
      pz(:,s)=Rz(:,:,s).'*pz(:,s);
  end
  
  pz(1:4,:)=EFP(1:32,:,jj).'*pz;
  for k=2:8
      pz( (k-1)*4+(1:4),:)=EFP(1:4,1:4,jj).'*pz( (k-1)*4+(1:4),:);
  end
  %pz_track(:,:,i)=pz;
end

jj=1;
pz(1:4:end,:)=pz(1:4:end,:)-dfdy(1:2:end,:,jj);
pz(2:4:end,:)=pz(2:4:end,:)-dfdy(2:2:end,:,jj);
pz(1:4,:)=EFP(1:32,:,jj).'*pz;
for k=2:8
  pz( (k-1)*4+(1:4),:)=EFP(1:4,1:4,jj).'*pz( (k-1)*4+(1:4),:);
end

p=zeros(msize1/5*4,sweepnum);
for s=1:sweepnum
    Q=Qs(:,:,s);
    p(:,s)=-(Q.')\pz(:,s);
end

%% Computing derivatives
trfgrad=zeros(rfnum,1);
alphagrad=zeros(rfnum,1);

% p_n+1
for k=1:8
  p( (k-1)*4+(1:4),:)=BB1.'*p( (k-1)*4+(1:4),:);
end

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
%{
jj=1;
p(1:5:end,:)=p(1:5:end,:)-dfdy(1:2:end,:,jj);
p(2:5:end,:)=p(2:5:end,:)-dfdy(2:2:end,:,jj);

for s=1:sweepnum
  trfgrad(i)=trfgrad(i)-p(:,s).'*dAydtrf(:,s,(i-1)*3+3);
end
%}
grads=cat(2,alphagrad,trfgrad);
if basis_flag
   grads=cat(1,grads(:),dfdu(:)); 
end
end