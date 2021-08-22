function params=sample_parameter_space(n,TR,extras)
params=cell(n,1);
defaultm0s=.1;%0.3;%
defaultT1=1.6;%0.7810;%
defaultT2f=.065;%.065;%
defaultR=30;%30;30;%
defaultT2s=60e-6;%10e-6;10e-6;%
defaultB0=0;
defaultB1=1;
params{1}={defaultm0s,defaultT1,defaultT2f,defaultR,defaultT2s,defaultB0,defaultB1};


for i=2:n
    switch i<inf%i<=600 || i>1101
        case true
            m0s = -1;
            while (m0s<0)
                m0sflag=1;
                if ~isempty(extras)
                    if contains(extras,'m0s')
                       m0s=0.1;%str2num(extras(dd+4:dd+6));
                       m0sflag=0;
                    end
                end
                if m0sflag
                    m0s = 0.1 + 0.1 * randn();
                end
            end

            logT1 = -1;
            while (logT1<log10(.3) || log10(5)<logT1)
                logT1 = 0.2*randn() + log10(1.5);
            end

            logT2f = -100;
            while (10^logT2f < 0.01 || logT2f > logT1)
                logT2f = log10(10^logT1 * 0.1) + 0.2*randn();
            end

            T1 = 10^logT1;
            T2f = 10^logT2f;

            R = 30+30*randn();
            while (R<15 || R>100)
                R = 30+30*randn();
            end

            T2s = 150*10^(-6) + 150*10^(-6)*randn();
            while(T2s < 20*10^(-6) || T2s > 500*10^(-6))
                T2s =  150*10^(-6) + 150*10^(-6)*randn();
            end
            Bvary='gaussian';
            switch Bvary
                case 'pio5'
                    B0 = randn() * pi/TR / 5;
                    B1 = 1 + randn() * 0.1;
                case 'pio2'
                    B0 = randn() * pi/TR / 2;
                    B1 = 1 + randn() * 0.1;
                case 'fixed'
                    B0 = 0;
                    B1 = 1;
                case 'gaussian'
                    B0 = randn()*pi/TR;
                    B1 = 0.9+randn()*0.15;
                case 'uniform'
                    B0 = -2*pi/TR+4*pi/TR*rand(1);
                    B1 = 0.9+randn()*0.15;
                case 'uniform-b1cutoff'
                    B0 = -2*pi/TR+4*pi/TR*rand(1);
                    %B1 = 0.9+randn()*0.15;
                    B1=0;
                    while B1<0.6 ||B1>1.3
                        B1 = 1+randn()*0.2;
                    end
                case 'uniform-2'
                    B0 = -1.4*pi/TR +((-0.6+1.4)*pi/TR)*rand(1); 
                    B1 = 0.9+randn()*0.15;
                case 'uniform-3'
                    B0 = 0.6*pi/TR+((1.4-0.6)*pi/TR)*rand(1);
                    B1 = 0.9+randn()*0.15;
                case 'uniform-B1fix'
                    B0 = -2*pi/TR+4*pi/TR*rand(1);
                    B1 = 1;
                otherwise
                    sprintf('please secify B0 and B1')

            end
        case false
            B0 = -2*pi/TR+4*pi/TR*( (i-601)/500);
            m0s=defaultm0s;
            T1=defaultT1;
            T2f=defaultT2f;
            R=defaultR;
            T2s=defaultT2s;
            B1=defaultB1;
    end
    params{i}={m0s,T1,T2f,R,T2s,B0,B1};
end


end