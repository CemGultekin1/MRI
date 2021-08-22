function [x,  fval, exitflag,output]=fmincon_optimize_basis(readfolder,filename,basis_size,T,extras)
%warning('off', 'MATLAB:nearlySingularMatrix');
%readfolder=varargin{1};
%filename=varargin{2};
%basis_size=varargin{3};
sweepnum=2;
S=dir([readfolder,'/fin*']);
TR=3.5e-3;
rfnum=ceil(T/TR);
history=struct('fval',[],'crb',[],'basis',[],'iter',[],'funccount',[],'stepsize',[],'firstorderopt',[]);
save(filename,'history');
%{
S=dir([readfolder,'/fin*']);
history=struct('fval',[],'crb',[],'basis',[],'iter',[],'funccount',[],'stepsize',[],'firstorderopt',[]);
save(filename,'history');
load([S(1).folder,'/',S(1).name],'data');
rfnum=size(data.fpg,1)/6-1;

splinenum=ceil(5/666*rfnum);
SB=Bloch.spline_basis(linspace(0,1,rfnum+1).',4*splinenum,'C0-spline');
[SB,~]=qr(SB,0);
sb2=size(SB,2);
%}


splinenum=ceil(10/666*rfnum);%ceil(50/666*rfnum);
SB=Bloch.spline_basis(linspace(0,1,rfnum+1).',ceil(2*splinenum/5),'C0-spline');
[SB,~]=qr(SB,0);



options = optimoptions('fmincon',...
        'Display','iter',...
        'MaxIter',3e3,...
        'MaxFunEvals',inf,...
        'Algorithm','interior-point',...
        'CheckGradients',true,...%'HessianFcn',@(x,lambda) hessianfcn(x,lambda,sim2,smc,alpha,trfmax,fixedtrf),...
        'SubproblemAlgorithm', 'cg',...%'factorization',...
        'SpecifyObjectiveGradient', true,...
        'SpecifyConstraintGradient',true,...
        'HonorBounds', false,...
        'CheckGradients', false,...
        'UseParallel', true,...
        'Output',@(x,optimValues,state) myoutfun(x,optimValues,state,sweepnum,rfnum,filename,SB),...
        'StepTolerance', 1e-9,...
        'ConstraintTolerance',1e-3,...
        'SpecifyConstraintGradient',true);%'ScaleProblem',true);


basis=randn(2*sweepnum*(rfnum+1),basis_size)/sqrt(2*sweepnum*(rfnum+1));
basis=reshape(basis,2,rfnum+1,sweepnum,[]);% 2,rfnum+1,3,bs
basis=permute(basis,[2,1,3,4]);  % rfnum+1,2,3,bs
x0=SB.'*reshape(basis,rfnum+1,[]);  % sb2,2,3,bs
costfun=@(x) fmincon_costfun(x,S,sweepnum,rfnum,filename,SB,extras);
nnlnr=@(x) nonlinear_constraint(x,S,sweepnum,rfnum,SB); 
[x, fval, exitflag,output,lambda,grad,hessian] = fmincon(costfun,x0,[],[],[],[],[],[],nnlnr,options); %
save([filename,'_finalized'],'x','fval','exitflag','output');
end
function [C,ceq,f,geq]=nonlinear_constraint(x0,S,sweepnum,rfnum,SB)
sbnum=size(SB,2);
bs=length(x0(:))/(sbnum)/2/sweepnum;
x0=reshape(x0,[],bs);
[QQ,RR]=qr(x0,0);
QQ=SB*reshape(QQ,sbnum,[]);  % rfnum+1,2,sweepnum,bs
QQ=reshape(QQ,rfnum+1,2,sweepnum,[]);
QQ=permute(QQ,[2,1,3,4]);
QQ=reshape(QQ,[],size(QQ,4)); % 2*rfnum+1*3,bs
C=0;
f=zeros(sbnum*2*sweepnum,bs);
totsnum=0;
for i=1:length(S)
    load([S(i).folder,'/',S(i).name],'fingerprints');
    fpg=fingerprints;
    snum=size(fpg,2);% (2,8,rf+1,sweepnum),snum
    %fpg=reshape(fpg,[],8,snum); 
    fpg=reshape(fpg,[],snum); 
    totsnum=totsnum+snum;
    CC=zeros(snum,1);
    ff=zeros(sbnum*2*sweepnum,bs,snum);
    parfor j=1:snum
        %fpgj=fpg(:,:,j);
        fingerprints=fpg(:,j);
        
        
        
         % (2,8,rf+1,sweepnum)
        fng=fingerprints;% 2,8,rf+1,sweepnum
        fng=reshape(fng,2,8,[]);% 2,8,rf+1,sweepnum
        fng=permute(fng,[2,1,3]); % 8,2,rf+1,sweepnum
        %phase_weights=reshape(phase_weights,1,[]);
        %fng=reshape(fng,[],sweepnum).*phase_weights;
        fpgj=reshape(fng,8,[]).'; % 2,rf+1,sweepnum,8

        fng=fpgj(:,1);
        bfng=QQ*(QQ.'*fng);
        
        CC(j)=-norm(bfng)/norm(fng);
        
        dfdy=bfng;
        W=QQ*RR/(RR.'*RR);
        dfdu=dfdy*(fng.'*W)+fng*(dfdy.'*W)-QQ*((QQ.'*dfdy)*(fng.'*W)+(QQ.'*fng)*(dfdy.'*W));
        dfdu=-dfdu/norm(fng)/norm(bfng);
        
        dfdu=reshape(dfdu,2,rfnum+1,sweepnum,[]);
        dfdu=permute(dfdu,[2,1,3,4]);
        dfdu=SB.'*reshape(dfdu,rfnum+1,[]);
        ff(:,:,j)=reshape(dfdu,[],bs);
    end
    f=f+sum(ff,3);
    C=C+sum(CC);
end
f=f(:)/totsnum;
C=C/totsnum+0.9;
ceq=0;
geq=f.*0;
end
function [C,f]=fmincon_costfun(x0,S,sweepnum,rfnum,filename,SB,extras)
sbnum=size(SB,2);
bs=length(x0(:))/(sbnum)/2/sweepnum;
x0=reshape(x0,[],bs);
[QQ,RR]=qr(x0,0);
QQ=SB*reshape(QQ,sbnum,[]);  % rfnum+1,2,sweepnum,bs
QQ=reshape(QQ,rfnum+1,2,sweepnum,[]);
QQ=permute(QQ,[2,1,3,4]);
QQ=reshape(QQ,[],size(QQ,4)); % 2*rfnum+1*sweepnum,bs

C=0;
f=zeros(sbnum*2*sweepnum,bs);
totsnum=0;
i=1;
load([S(i).folder,'/',S(i).name],'fingerprints');
fpg=fingerprints;
snum=size(fpg,2);%/8;
CRBVEC=zeros(8,length(S)*snum);
for i=1:length(S)
    if i>1
        load([S(i).folder,'/',S(i).name],'fingerprints');
        fpg=fingerprints;
        snum=size(fpg,2);%/8;
    end
    %fpg=reshape(fpg,[],8,snum);
    totsnum=totsnum+snum;
    CC=zeros(snum,1);
    ff=zeros(sbnum*2*sweepnum,bs,snum);
    %for j=1:snum
    for j=1:snum
        fingerprints=fpg(:,j);
        
        % (2,8,rf+1,sweepnum)
        fng=fingerprints;% 2,8,rf+1,sweepnum
        fng=reshape(fng,2,8,[]);% 2,8,rf+1,sweepnum
        fng=permute(fng,[2,1,3]); % 8,2,rf+1,sweepnum
        %phase_weights=reshape(phase_weights,1,[]);
        %fng=reshape(fng,[],sweepnum).*phase_weights;
        fpgj=reshape(fng,8,[]).'; % 2,rf+1,sweepnum,8
        bfng=QQ*(QQ.'*fpgj);
        
        F=bfng.'*bfng;
        I=eye(8);
        CRBVEC(:,(i-1)*snum+j)=diag(inv(F));
        FI=F\I(:,2:4);
        CC(j)=sum(diag(FI(2:4,:))); % Relative CRB values

        if contains(extras,'T1pT2f')
            FI=F\I(:,3:4);
            C_=diag(FI(3:4,:));
            CC(j)=sum(C_);
        elseif contains(extras, 'R-only')
            FI=F\I(:,5);
            C_=diag(FI(5,:));
            CC(j)=sum(C_);
        elseif contains(extras, 'default')
            FI=F\I(:,2:4);
            C_=diag(FI(2:4,:));
            CC(j)=sum(C_);
        end
        
        
        FF=FI*FI.';
        dfdy=-2*bfng*FF;
        W=QQ*RR/(RR.'*RR);
        dfdu=dfdy*(fpgj.'*W)+fpgj*(dfdy.'*W)-QQ*((QQ.'*dfdy)*(fpgj.'*W)+(QQ.'*fpgj)*(dfdy.'*W));
        
        dfdu=reshape(dfdu,2,rfnum+1,sweepnum,[]);
        dfdu=permute(dfdu,[2,1,3,4]);
        dfdu=SB.'*reshape(dfdu,rfnum+1,[]);
        ff(:,:,j)=reshape(dfdu,[],bs);
    end
    f=f+sum(ff,3);
    C=C+sum(CC);
end
load(filename,'history');
history.crb=CRBVEC;
save(filename,'history');
f=f/totsnum;
C=C/totsnum;
end

function stop = myoutfun(x0,optimValues,state,sweepnum,rfnum,filename,SB),...
stop = false;
load(filename,'history');

basis_visu=SB*reshape(x0,size(SB,2),[]);
basis=reshape(basis_visu,rfnum+1,2,sweepnum,[]);
basis=permute(basis,[2,1,3,4]);
basis=reshape(basis,[],size(basis,4)); % 2 x rfnum+1 x sweepnum x bs
history.basis=basis;

history.iter=[history.iter; optimValues.iteration];
history.funccount=[history.funccount; optimValues.funccount];
history.fval = [history.fval; optimValues.fval];
%history.crb = [history.crb; optimValues.fval];
history.stepsize=[history.stepsize;optimValues.stepsize];
history.firstorderopt=[history.firstorderopt;optimValues.firstorderopt];

clf;
basis_visu=reshape(basis_visu,(rfnum+1)*2*sweepnum,[]);
plot(basis_visu);
pause(0.05);

switch state
   case 'done'
       stop=true;
   otherwise
end
save(filename,'history');
end
