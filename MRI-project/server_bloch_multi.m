function server_bloch_multi(ID,cpu_num)
[T,bss,SARMAX,writefile,extras]=server_opt_ID(ID);
readfile=[];

sample_num=500;
server_bloch_optimization(T,bss,SARMAX,sample_num,cpu_num,readfile,writefile,extras)
end