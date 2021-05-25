function [T,bss,SARMAX,extras]=server_opt_ID(ID)
SARMAX_v5=2.8278e+05*2;
% primary concern:
extras='';
if ID<=6 && 1<=ID
    T=4;
    ID2=ceil(ID/3);
    bs=[9,13];
    SARV=1e-2*[60,75,90];
    bss=bs(ID2);
    ID3=mod(ID-1,3)+1;
    SARMAX=SARV(ID3)*SARMAX_v5;
    Tstr=strrep(num2str(T),'.','p');
    sarstr=num2str(SARMAX,'%.2e');
    sarstr=strrep(strrep(sarstr,'.','p'),'+','');
    readfile=['2021_01_12_T',Tstr,'_b',num2str(bss),'_SAR',sarstr,'.mat'];
    if exist(readfile,'file')==2
        load(readfile,'history');
        control=history.control;
    else
        control=[];
    end
    TR=3.5e-3;
    TRFmax=5e-4;
    tag=['ID-',num2str(ID)];%['2021_01_12_T',Tstr,'_b',num2str(bss),'_SAR',sarstr];
elseif ID<=22
    SARMAX=0.6*SARMAX_v5;
    bss=9;
    tt=1:.5:9;
    [~,I]=sort(abs(tt-4.5),'ascend');
    tt=tt(I);
    II=find(mod(tt,1)==0.5);
    I=find(mod(tt,1)==0);
    tt=cat(2,tt(I),tt(II));
    tt(tt==4)=[];
    T=tt(ID-6);
    
    Tstr=strrep(num2str(T),'.','p');
    sarstr=num2str(SARMAX,'%.2e');
    sarstr=strrep(strrep(sarstr,'.','p'),'+','');
    readfile=['2021_01_12_T',Tstr,'_b',num2str(bss),'_SAR',sarstr,'.mat'];
    if exist(readfile,'file')==2
        load(readfile,'history');
        control=history.control;
    else
        control=[];
    end
     TR=3.5e-3;
     TRFmax=5e-4;
     tag=['ID-',num2str(ID)];%['2021_01_12_T',Tstr,'_b',num2str(bss),'_SAR',sarstr];
elseif ID>=23 && ID<=26
    ID2=ID-22;
    load('+Bloch/past_versions','past_versions');
    TR=past_versions(ID2).TR;
    if TR>3.5e-3
        TRFmax=1e-3;
    else
        TRFmax=5e-4;
    end
    control=past_versions(ID2).control;
    control=[control(1:end-1,1)+control(2:end,1),control(1:end-1,2)];
    rfnum=size(control,1);
    T=TR*(rfnum+1);
    bss=[];
    SARMAX=[];
    tag=past_versions(ID2).name;
elseif 37>=ID  && ID>=27
    SARMAX=0.6*SARMAX_v5;
    bss=9;
    T=(ID-26)+9;
    
    Tstr=strrep(num2str(T),'.','p');
    sarstr=num2str(SARMAX,'%.2e');
    sarstr=strrep(strrep(sarstr,'.','p'),'+','');
    readfile=['2021_01_16_T',Tstr,'_b',num2str(bss),'_SAR',sarstr,'.mat'];
    if exist(readfile,'file')==2
        load(readfile,'history');
        control=history.control;
    else
        control=[];
    end
     TR=3.5e-3;
     TRFmax=5e-4;
     tag=['ID-',num2str(ID)];%['2021_01_16_T',Tstr,'_b',num2str(bss),'_SAR',sarstr];
elseif 38<=ID
    SARMAX=0.6*SARMAX_v5;
    bss=9;
    T=4;
    switch ID
        case 38
            extras='default m0s=0.1 tag=default';
        case 39
            extras='T1+T2f m0s=0.1 tag=T1T2f';
        case 40
            extras='R-only tag=R';
        case 41
            extras='T1pT2f m0s=0.1 tag=T1pT2f';
    end
end
end