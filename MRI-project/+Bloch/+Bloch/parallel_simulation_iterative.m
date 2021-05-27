function [c,g,y01,p01]=parallel_simulation_iterative(samples,control,TR,TRFmax,order,y01,p01,basis)
if ~exist('y01','var')
    y01=cell(length(samples),1);
elseif length(y01)<length(samples)
    y01=cell(length(samples),1);
end

if ~exist('p01','var')
    p01=cell(length(samples),1);
elseif length(p01)<length(samples)
    p01=cell(length(samples),1);
end


sweep_phase=[-0.3,0,0.3]*pi;
c=zeros(1,length(samples));

if ~isempty(basis)
    rfnum=size(control,1);
    fnglen=(rfnum+1)*2*3;
    basis=reshape(basis,fnglen,[]);
    basis_size=size(basis,2);
    [QQ,RR]=qr(basis,0);
    basis={QQ,RR};
end
if order==2
    if ~isempty(basis)
        totlen=basis_size*fnglen+size(control,1)*2;
        g=zeros(totlen,length(samples));
    else
        g=zeros(size(control,1)*2,length(samples));
    end
    parfor i=1:length(samples)
        %[~,y01,~,~,CTOTh(i),~]
        [~,y01{i},p01{i},~,c(i),grads]=Bloch.simulate_iterative(samples{i},TR,TRFmax,sweep_phase,control,order,y01{i},p01{i},basis);
        g(:,i)=grads(:);
    end
else
    g=[];
    p01=[];
    parfor i=1:length(samples)
        [~,y01{i},~,~,c(i),~]=Bloch.simulate_iterative(samples{i},TR,TRFmax,sweep_phase,control,order,y01{i},p01{i},basis);
    end
end
c=mean(c);
if nargout==2
    g=mean(g,2);
    %g=reshape(g,[],2);
end
end