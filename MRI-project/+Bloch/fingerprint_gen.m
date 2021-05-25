function data=fingerprint_gen(sample_num1,sample_num2,control,TR,TRFmax)
if ~exist('TR','var')
    TR=3.5e-3;
    TRFmax=5e-4;
end
sweep_phase=[-0.3,0,0.3]*pi;
sweepnum=length(sweep_phase);
rng(0);
samples=Bloch.sample_parameter_space(sample_num2,TR);
samples={samples{sample_num1+1:sample_num2}};
sample_num=sample_num2-sample_num1;
order=1;
C=zeros(3,sample_num);
rfnum=size(control,1);
fpglen=2*8*(rfnum+1)*sweepnum;
fng=zeros(fpglen,sample_num);

parfor i=1:sample_num
     %[ff,C(:,i),~]=Bloch.simulate_multishooting(samples{i},TR,TRFmax,sweep_phase,control,order);
     [ff,C(:,i),~]=Bloch.simulate(samples{i},TR,TRFmax,sweep_phase,control,order);
     fng(:,i)=reshape(ff,[],1);
end
fng=reshape(fng,2,8,[],sample_num);% 2,8,rf+1,sweepnum,sample_num
fng=permute(fng,[1,3,2,4]); % 2,rf+1,sweepnum,8,sample_num
fng=reshape(fng,[],8*sample_num);
data=struct('CRB',C,'fpg',fng,'smpl',[sample_num1+1,sample_num2]);
%{
[u,s,v]=svd(fng,0);
I=find(diag(s)<s(1,1)/1e2,1);
u=u(:,1:I);
v=v(:,1:I);
s=s(1:I,1:I);

data=struct('CRB',C,'fpg_u',u,'fpg_v',v,'fpg_s',s,'smpl',[sample_num1+1,sample_num2]);
%}
end