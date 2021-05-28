function C=statistics(sample_num1,sample_num2,control,basis,extras)
TR=3.5e-3;
TRFmax=5e-4;
sweep_phase=[-0.3,0.3]*pi;
rng(0);
samples=Bloch.sample_parameter_space(sample_num2,TR,extras);
samples={samples{sample_num1+1:sample_num2}};
sample_num=sample_num2-sample_num1;
order=1;
CRBs=cell(1);
%{
if pastver==1
    load('+Bloch/past_versions.mat','past_versions');
    control=reshape(control,[],2);
    CRBs=cell(length(past_versions)+1,1);
    for j=1:length(past_versions)
        control1=past_versions(j).control;
        TR1=past_versions(j).TR;
        alphaseq=control1(1:end-1,1)+control1(2:end,1);
        control1=cat(2,alphaseq,control1(1:end-1,2));
        switch TR1
            case 4.5e-3
                TRFmax1=1e-3;
            case 3.5e-3
                TRFmax1=5e-4;
        end
        C=zeros(3,sample_num);
        disp(['          computing version index: ' , num2str(j),'/',num2str(length(past_versions)), '...']);
        parfor i=1:sample_num
             [~,C(:,i),~]=Bloch.simulate_multishooting(samples{i},TR1,TRFmax1,sweep_phase,control1,order,basis);
        end
        CRBs{j}=C;
    end
else
    CRBs=cell(1,1);
end
C=zeros(3,sample_num);
if pastver
    disp(['          computing version index: ' , num2str(length(past_versions)),'/',num2str(length(past_versions)), '...']);
end
%}


parfor i=1:sample_num
    [~,C(:,i)]=Bloch.simulate(samples{i},TR,TRFmax,sweep_phase,control,order,basis,extras);
end
%{
if pastver
    for i=1:length(CRBvals)-1
        control1=past_versions(i).control;
        TR1=past_versions(j).TR;
        rfnum1=size(control1,1);
        CRBvals(i).name=past_versions(i).name;
        CRBvals(i).CRB=CRBvals(i).CRB*(rfnum1+1)*TR1;
    end
end
%}
C=C*(size(control,1)+1)*TR;
end