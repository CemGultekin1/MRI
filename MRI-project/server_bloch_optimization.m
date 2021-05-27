function server_bloch_optimization(T,basis_size,SARMAX,sample_num,cpu_num,readfile,writefile,extras)
%T=4;
maxiter=5e3;
if cpu_num>0
    parpool(cpu_num);
end
TR=3.5e-3;
TRFmax=5e-4;
%sample_num=500;
%smc=50;
rfnum=ceil(T/TR);
rng(0);
alphaseq=(rand(rfnum,1)*pi+pi)*0.2;
TRFseq=rand(rfnum,1)*(TRFmax-5e-5)+5e-5;
X=[alphaseq,TRFseq];
Tstr=strrep(num2str(T),'.','p');
samples=Bloch.sample_parameter_space(sample_num,TR,extras);
sweepnum=2;
basis=randn((rfnum+1)*2*sweepnum,basis_size)/sqrt((rfnum+1)*2*sweepnum);
Bloch.fmincon_optimize(...
            readfile,writefile,...
            maxiter,...
            'SARMAX',SARMAX,...
            'control',X,...
            'TR',TR,...
            'TRFmax',TRFmax,...
            'samples',samples,...
            'basis',basis,...
            'extras',extras);
end

