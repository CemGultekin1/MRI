function [x,  fval, exitflag,output]=fmincon_optimize(varargin)
%warning('off', 'MATLAB:nearlySingularMatrix');
readfile=varargin{1};
filename=varargin{2};
maxiter=varargin{3};
I=find(strcmp(varargin,'control'));
x0=varargin{I+1};
x0=reshape(x0,[],2);
rfnum=size(x0,1);

I=find(strcmp(varargin,'TR'));
TR=varargin{I+1};

I=find(strcmp(varargin,'TRFmax'));
trfmax=varargin{I+1};

I=find(strcmp(varargin,'samples'));
samples=varargin{I+1};

I=find(strcmp(varargin,'basis'));
if ~isempty(I)
    basis=varargin{I+1};
    basis_size=size(basis,2);
    basis_opt=1;
else
    basis_opt=0;
end

I=find(strcmp(varargin,'SARMAX'));
if ~isempty(I)
    SARMAX=varargin{I+1};
else
    SARMAX=inf;
end

I=find(strcmp(varargin,'extras'));
if ~isempty(I)
    extras=varargin{I+1};
else
    extras='';
end
if ~isempty(readfile)
    S=dir([readfile,'*']);
else
    S=[];
end
if isempty(S)
    history=struct('fval',[],'crb',[],'control',[],'basis',[],'iter',[],'funccount',[],'stepsize',[],'firstorderopt',[]);
else
    load(S(1).name,'history');
    if size(history.control,2)>=2
        x0=history.control;
        if basis_opt
            if size(history.basis,2)==basis_size
                basis=history.basis;
            end
        end
    else
        history=struct('fval',[],'crb',[],'control',[],'basis',[],'iter',[],'funccount',[],'stepsize',[],'firstorderopt',[]);
    end
end
save(filename,'history');
splinenum=ceil(10/666*rfnum);%ceil(50/666*rfnum);
sweepnum=2;
S=Bloch.spline_basis(linspace(0,1,rfnum).',splinenum,'hann');
[S,~]=qr(S,0);
splinenum2=size(S,2);

SB=Bloch.spline_basis(linspace(0,1,rfnum+1).',ceil(2*splinenum/5),'C0-spline');
[SB,~]=qr(SB,0);

B1=4500;

%lb=[zeros(rfnum,1);5e-5*ones(rfnum,1)].*0;
%ub=[pi*ones(rfnum,1);5e-4*ones(rfnum,1)].*0+1;
%D=eye(rfnum);
%A=[D,-eye(rfnum)*B1];
A=[S*pi,-eye(rfnum)*B1*(trfmax-5e-5)];
Alb=[-S,zeros(rfnum)];
Aub=[S,zeros(rfnum);];
A=cat(1,A,Alb,Aub);
%xsc=[pi*ones(rfnum,1);(trfmax-5e-5)*ones(rfnum,1)];
%xoff=[zeros(rfnum,1);5e-5*ones(rfnum,1)];
%b=-A*xoff;
%A=A*diag(xsc);
b=cat(1,B1*ones(rfnum,1)*5e-5,zeros(rfnum,1),ones(rfnum,1));
lb=cat(1,-inf*ones(splinenum2,1),zeros(rfnum,1));
ub=cat(1,inf*ones(splinenum2,1),ones(rfnum,1));
if basis_opt
    %fng_len=size(basis,1);
    fng_len=size(SB,2)*4;
    A=cat(2,A,zeros(size(A,1),fng_len*basis_size));
    lb=cat(1,lb,-inf*ones(fng_len*basis_size,1));
    ub=cat(1,ub,inf*ones(fng_len*basis_size,1));
end

options = optimoptions('fmincon',...
        'Display','iter',...
        'MaxIter',maxiter,...
        'MaxFunEvals',inf,...
        'Algorithm','interior-point',...
        'CheckGradients',true,...%'HessianFcn',@(x,lambda) hessianfcn(x,lambda,sim2,smc,alpha,trfmax,fixedtrf),...
        'SubproblemAlgorithm', 'cg',...
        'SpecifyObjectiveGradient', true,...
        'SpecifyConstraintGradient',true,...
        'HonorBounds', false,...
        'CheckGradients', false,...
        'UseParallel', true,...
        'Output',@(x,optimValues,state) myoutfun(x,optimValues,state,filename,rfnum,sweepnum,trfmax,S,SB),...
        'StepTolerance', 1e-9,...
        'ConstraintTolerance',1e-3,...
        'SpecifyConstraintGradient',true);%'ScaleProblem',true);
x0=reshape(x0,[],2);

x0(:,1)=x0(:,1)/pi;
x0(:,2)=(x0(:,2)-5e-5)/(trfmax-5e-5);
x0=cat(1,S.'*x0(:,1),x0(:,2));
if basis_opt
     % 2,rfnum+1,sweepnum,8
     basis=reshape(basis,2,rfnum+1,sweepnum,[]);
     basis=permute(basis,[2,1,3,4]);
     basis=SB.'*reshape(basis,rfnum+1,[]);
     x0=cat(1,x0,basis(:));
end
costfun=@(x) fmincon_costfun(x,rfnum,sweepnum,samples,TR,trfmax,S,SB,extras);
nnlnr=@(x) sarmax_const(x,rfnum,samples,TR,trfmax,S,SARMAX);
[x, fval, exitflag,output,lambda,grad,hessian] = fmincon(costfun,x0,A,b,[],[],lb,ub,nnlnr,options); %
basis=reshape(x(rfnum+splinenum*2+1:end),(rfnum+1)*2*sweepnum,[]);
x=x(1:rfnum+splinenum*2);
x=cat(2,S*x(1:splinenum*2),x(splinenum*2+1:end));
x(:,1)=x(:,1)*pi;
x(:,2)=x(:,2)*(trfmax-5e-5)+5e-5;

if basis_opt
    save([filename,'_finalized'],'x','basis','fval','exitflag','output');
else
    save([filename,'_finalized'],'x','fval','exitflag','output');
end
end
function [c,ceq,g,geq]=sarmax_const(x0,rfnum,samples,TR,trfmax,S,SARMAX)
splnum=size(S,2);
x=x0(1:rfnum+splnum);
x=cat(2,S*x(1:splnum),x(splnum+1:end));
x(:,1)=x(:,1)*pi;
x(:,2)=x(:,2)*(trfmax-5e-5)+5e-5;
T=(rfnum+1)*TR;
c=sum( x(:,1).^2./x(:,2))/T/SARMAX-1;
g=zeros(rfnum,2);
g(:,1)=2/T/SARMAX*x(:,1)./x(:,2);
g(:,2)=-1/T/SARMAX*x(:,1).^2./x(:,2).^2;
g(:,1)=g(:,1)*pi;
g(:,2)=g(:,2)*(trfmax-5e-5);
g=cat(1,S.'*g(:,1),g(:,2),zeros(length(x0)-splnum-rfnum,1));
ceq=0;
geq=g.*0;
end
function [C,G]=fmincon_costfun(x0,rfnum,sweepnum,samples,TR,trfmax,S,SB,extras)
splnum=size(S,2);
x=x0(1:rfnum+splnum);
x=cat(2,S*x(1:splnum),x(splnum+1:end));
x(:,1)=x(:,1)*pi;
x(:,2)=x(:,2)*(trfmax-5e-5)+5e-5;
basis=SB*reshape(x0(rfnum+splnum+1:end),size(SB,2),[]);
basis=reshape(basis,rfnum+1,2,sweepnum,[]);
basis=permute(basis,[2,1,3,4]);
basis=reshape(basis,[],size(basis,4));
g=0;c=0;
switch nargout
    case 1
        %[C,~]=Bloch.parallel_simulation_multishooting(samples,x,TR,trfmax,1,basis);
        [C,~]=Bloch.parallel_simulation(samples,x,TR,trfmax,1,basis,extras);
    case 2
        %[C,G]=Bloch.parallel_simulation_multishooting(samples,x,TR,trfmax,2,basis);
        [C,G]=Bloch.parallel_simulation(samples,x,TR,trfmax,2,basis,extras);
end
if nargout>=2
    G(1:rfnum)=G(1:rfnum)*pi;
    G(rfnum+1:2*rfnum)=G(rfnum+1:2*rfnum)*(trfmax-5e-5);
    G(1:rfnum)=G(1:rfnum)+g;
    gbasis=G(2*rfnum+1:end);
    gbasis=reshape(gbasis,2,rfnum+1,sweepnum,[]);
    gbasis=permute(gbasis,[2,1,3,4]);
    gbasis=SB.'*reshape(gbasis,rfnum+1,[]);
    G=cat(1,S.'*G(1:rfnum),G(rfnum+1:2*rfnum),gbasis(:));
end
C=C+c;
end
function [c,g]=ptotvar(x,delta)
c=sum(sqrt((x(2:end)-x(1:end-1)).^2+delta^2)-delta);
g=zeros(length(x),1);
for i=1:length(x)
    if i>1
        di=sqrt((x(i)-x(i-1))^2+delta^2);
        gi=(x(i)-x(i-1))/di;
    else
        gi=0; 
    end
    if i<length(x)
        di1=sqrt((x(i)-x(i+1))^2+delta^2);
        gi1=(x(i)-x(i+1))/di1;
    else
        gi1=0;
    end
    g(i)=gi+gi1;
end

end
function stop = myoutfun(x0,optimValues,state,filename,rfnum,sweepnum,trfmax,S,SB)
stop = false;
load(filename,'history');

splnum=size(S,2);
x=x0(1:rfnum+splnum);
x=cat(2,S*x(1:splnum),x(splnum+1:end));
x(:,1)=x(:,1)*pi;
x(:,2)=x(:,2)*(trfmax-5e-5)+5e-5;
if length(x0)>rfnum+splnum
    basis_visu=SB*reshape(x0(rfnum+splnum+1:end),size(SB,2),[]);
    basis=reshape(basis_visu,rfnum+1,2,sweepnum,[]);
    basis=permute(basis,[2,1,3,4]);
    basis=reshape(basis,[],size(basis,4));
    history.basis=basis;
end

history.control = x;
history.iter=[history.iter; optimValues.iteration];
history.funccount=[history.funccount; optimValues.funccount];
history.fval = [history.fval; optimValues.fval];
history.crb = [history.crb; optimValues.fval];
history.stepsize=[history.stepsize;optimValues.stepsize];
history.firstorderopt=[history.firstorderopt;optimValues.firstorderopt];

clf;
subplot(2,2,1);
plot(x(:,1));
subplot(2,3,4);
semilogy(x(:,2));
subplot(1,2,2);
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
