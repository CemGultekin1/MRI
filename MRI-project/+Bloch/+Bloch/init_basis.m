function U=init_basis(sample_num,TR,TRFmax,control,basis_size)
control=reshape(control,[],2);
rfnum=size(control,1);
samples=Bloch.sample_parameter_space(sample_num,TR);
fng=zeros((rfnum+1)*2*3*8,sample_num);
order=1;
sweep_phase=[-0.3,0,0.3]*pi;
parfor i=1:sample_num
    [fps,~,~]=Bloch.simulate_multishooting(samples{i},TR,TRFmax,sweep_phase,control,order);
    fps=reshape(fps,2,8,[],3); % 2,8,rfnum+1,3
    fps=permute(fps,[1,3,4,2]); % 2,rfnum+1,3,8
    fps=reshape(fps,[],8);% 2,rfnum+1,3,8
    fng(:,i)=fps(:); % 2,rfnum+1,3,8
end
fng=reshape(fng,(rfnum+1)*2*3,[]);
[Q,~]=qr(fng,0);
[U,~,~]=svd(Q.'*fng,0);
U=Q*U(:,1:basis_size);
end