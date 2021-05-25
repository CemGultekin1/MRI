function server_bloch_svd_basis(ID,max_bss)
writeaddress=['2020_ID_',num2str(ID),'_data/fin*'];
S=dir(writeaddress);
m=600;
rng(0);
snum=0;
for i=1:length(S)
   load([S(i).folder,'/',S(i).name],'data');
   fpg=data.fpg(:,1:8:end);
   snum=size(fpg,2);
   fpg=reshape(fpg,[],3,snum);
   fpg=fpg(:,[1,3],:);
   fpg=reshape(fpg,[],snum);
   if ~exist('X','var')
       X=randn(size(fpg,2),m);
       Y=zeros(size(fpg,1),m);
   else
       X=randn(size(fpg,2),m);
   end
   Y=Y+fpg*X;
   snum=snum+size(fpg,2);
end
[Q,~]=qr(Y,0);
X=zeros(m,snum);
for i=1:length(S)
    load([S(i).folder,'/',S(i).name],'data');
    fpg=data.fpg(:,1:8:end);
    snum=size(fpg,2);
    fpg=reshape(fpg,[],3,snum);
    fpg=fpg(:,[1,3],:);
    fpg=reshape(fpg,[],snum);
    X(:,data.smpl(1):data.smpl(2))=Q.'*fpg;
end
[U,~,~]=svd(X,0);
U=Q*U;
U=U(:,1:max_bss);
basis=U;
save(['2020_ID_',num2str(ID),'_data/svd_basis_',num2str(max_bss)],'basis');
end