function [c,g]=parallel_simulation_multishooting(samples,control,TR,TRFmax,order,basis)

sweep_phase=[-0.3,0.3]*pi;
sweepnum=length(sweep_phase);
c=zeros(1,length(samples));

if ~isempty(basis)
    rfnum=size(control,1);
    fnglen=(rfnum+1)*2*sweepnum;
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
        [~,~,c(i),grads]=Bloch.simulate_multishooting(samples{i},TR,TRFmax,sweep_phase,control,order,basis);
        g(:,i)=grads(:);
    end
else
    parfor i=1:length(samples)
        [~,~,c(i)]=Bloch.simulate_multishooting(samples{i},TR,TRFmax,sweep_phase,control,order,basis);
    end
end
c=mean(c);
if order==2
    g=mean(g,2);
    %g=reshape(g,[],2);
end
end