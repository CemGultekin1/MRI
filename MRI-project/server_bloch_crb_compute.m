function server_bloch_crb_compute(ID,bsssize)
basis_address=['2020_ID_',num2str(ID),'_data/*basis*'];
data_address=['2020_ID_',num2str(ID),'_data/fin*'];

SB=dir(basis_address);
S=dir(data_address);
snum=0;
for i=1:length(S)
   load([S(i).folder,'/',S(i).name],'data');
   snum=snum+size(data.fpg,2)/8;
end
for k=1:length(SB)
    svdflag=0;
    crbflag=contains(SB(k).name,'CRB') || contains(SB(k).name,'sw2');
    if contains(SB(k).name,'svd') && ~crbflag
        load([SB(k).folder,'/',SB(k).name],'basis');
        
        svdflag=1;
    elseif contains(SB(k).name,'finalized') && ~crbflag
        basis=[];
    elseif ~crbflag
        load([SB(k).folder,'/',SB(k).name],'history');
        basis=history.basis;
        
    else
        basis=[];
    end
    if size(basis,2)<13
        basis=[];
    end
    if ~isempty(basis)
        if svdflag
            for ll=1:length(bsssize)
                bs=bsssize(ll);
                u=basis(:,1:bs);
                [u,~]=qr(u,0);
                CRB=zeros(4,snum);
                energy=zeros(5,snum);
                for i=1:length(S)
                   load([S(i).folder,'/',S(i).name],'data');
                   fpg=data.fpg;
                   fpg=reshape(fpg,size(fpg,1),8,[]);
                   snum=size(fpg,3);
                   fpg=reshape(fpg,[],3,8*snum);
                   fpg=fpg(:,[1,3],:);
                   fpg=reshape(fpg,[],8*snum);
                   %disp(['size u (',num2str(size(u,1)),',',num2str(size(u,2)),')  size fpg ',num2str(size(fpg,1))]);
                   ufpg=u*(u.'*fpg);
                   ufpg=reshape(ufpg,size(ufpg,1),8,[]);
                   fpg=reshape(fpg,size(fpg,1),8,[]);
                   II=eye(8);II=II(:,1:4);
                   CC1=zeros(4,size(ufpg,3));
                   EE1=zeros(5,size(ufpg,3));
                   for j=1:size(ufpg,3)
                       ff=ufpg(:,:,j);
                       C=(ff.'*ff)\II;
                       CC1(:,j)=diag(C(1:4,1:4));
                       
                       f=fpg(:,:,j);
                       C=(f.'*f)\II;
                       C=diag(C(1:4,1:4));
                       EE1(2:5,j)=sqrt(C./CC1(:,j));
                       EE1(1,j)=norm(ff(:,1))/norm(f(:,1));
                   end
                   CRB(:,data.smpl(1):data.smpl(2))=CC1;
                   energy(:,data.smpl(1):data.smpl(2))=EE1;
                end
                save([SB(1).folder,'/svd_basis_',num2str(size(u,2)),'_CRB.mat'],'CRB','energy');
            end
            CRB=zeros(3,snum);
            energy=ones(4,snum);
            for i=1:length(S)
               load([S(i).folder,'/',S(i).name],'data');
               CRB(:,data.smpl(1):data.smpl(2))=data.CRB;
            end
            save([SB(1).folder,'/no_basis_CRB.mat'],'CRB','energy');
        else
            
            u=basis;
            [u,~]=qr(u,0);
            CRB=zeros(4,snum);
            energy=zeros(5,snum);
            for i=1:length(S)
               load([S(i).folder,'/',S(i).name],'data');
               fpg=data.fpg;
               fpg=reshape(fpg,size(fpg,1),8,[]);
               snum=size(fpg,3);
               fpg=reshape(fpg,[],3,8*snum);
               fpg=fpg(:,[1,3],:);
               fpg=reshape(fpg,[],8*snum);
               %disp(['size u (',num2str(size(u,1)),',',num2str(size(u,2)),')  size fpg ',num2str(size(fpg,1))]);
               ufpg=u*(u.'*fpg);
               fpg=reshape(fpg,size(fpg,1),8,[]);
               ufpg=reshape(ufpg,size(ufpg,1),8,[]);
               II=eye(8);II=II(:,1:4);
               CC1=zeros(4,size(ufpg,3));
               EE1=zeros(5,size(ufpg,3));
               for j=1:size(ufpg,3)
                   ff=ufpg(:,:,j);
                   C=(ff.'*ff)\II;
                   CC1(:,j)=diag(C(1:4,1:4));%diag(C(2:4,1:3));
                   f=fpg(:,:,j);
                   C=(f.'*f)\II;
                   C=diag(C(1:4,1:4));
                   EE1(2:5,j)=sqrt(C./CC1(:,j));
                   EE1(1,j)=norm(ff(:,1))/norm(f(:,1));
               end
               CRB(:,data.smpl(1):data.smpl(2))=CC1;
               energy(:,data.smpl(1):data.smpl(2))=EE1;
            end
            save([SB(1).folder,'/optimized_basis_',num2str(size(u,2)),'_CRB.mat'],'CRB','energy');
        end
    end
end

end