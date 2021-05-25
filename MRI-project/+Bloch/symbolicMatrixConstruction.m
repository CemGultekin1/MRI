function symbolicMatrixConstruction()
syms m0s T1 T2f R T2s_ TR alpha0 TRF B0 B1 TRFmax RFflag

alpha=alpha0*B1*RFflag;
wz = B0;
wy = alpha / TRF;
pol = @(x) 0.4922812113644589 + 0.4734209035397098*x + 0.21429448654701067*x^2 + 0.061682063760768846*x^3 + 0.003405803678011424*x^4 - 0.013520010778091904*x^5 - 0.0042308073556300155*x^6 + 0.002510379993273342*x^7 + 0.0007252948004778164*x^8 - 0.0003617855549952223*x^9 - 5.339997364891677e-5*x^10 + 3.462515593522486e-5*x^11 + 5.75925759532023e-7*x^12 - 1.919567796635778e-6*x^13 + 1.5974668626567822e-7*x^14 + 4.950447003497e-8*x^15 - 9.135634780140215e-9*x^16 - 8.805182533780258e-11*x^17 + 1.4407177783963282e-10*x^18 - 1.337098474521444e-11*x^19 + 4.0352556351291685e-13*x^20;
R_rf = T2s_ * (alpha/TRF)^2 * pol(log(TRF/T2s_));
m0f=1-m0s;
A=[
    -1/T2f,       -wz,            wy,                          0,                0
        wz,    -1/T2f,             0,                          0,                0
       -wy,         0,    -1/T1-R*m0s,                     R*m0f,           m0f/T1
         0,         0,          R*m0s,          -1/T1-R_rf-R*m0f,           m0s/T1
         0,         0,             0,                          0,                0];


ms=size(A,1);
M=[A];
params='M  ';
for i=2:8
    if i==2
        D=diff(A,m0s)*0.1;
        params=[params,':  m0s  '];
    elseif i==3
        D=diff(A,T1)*T1;
        params=[params,':  T1  '];
    elseif i==4
        D=diff(A,T2f)*T2f;
        params=[params,':  T2f  '];
    elseif i==5
        D=diff(A,R)*R;
         params=[params,':  R  '];
    elseif i==6
        D=diff(A,T2s_)*T2s_;
        params=[params,':  T2s'];
    elseif i==7
        D=diff(A,B0);
        params=[params,':  B0'];
    elseif i==8
        D=diff(A,B1);
        params=[params,':  B1'];
    else
        error('coding error');
    end
    M=[M;D];
end
disp(['ODE states=  ',params]);

matlabFunction(M,'File','+Bloch/Mfun','Vars',[alpha0,TRF,m0s,T1,T2f,R,T2s_,B0,B1,TR,RFflag]);
dMdalpha=diff(M,alpha0);
dMdTRF=diff(M,TRF);
matlabFunction(dMdalpha,'File','+Bloch/dMdalphafun','Vars',[alpha0,TRF,m0s,T1,T2f,R,T2s_,B0,B1,TR,RFflag]);
matlabFunction(dMdTRF,'File','+Bloch/dMdTRFfun','Vars',[alpha0,TRF,m0s,T1,T2f,R,T2s_,B0,B1,TR,RFflag]);

attenuation= (B1 * pi)^2 * T2s_ / TRFmax;
boundary0=diag(1./[sin(pi/2*B1)^2,-sin(pi/2*B1)^2,cos(pi*B1),exp(- attenuation),1]);
boundary1=diag([-1,-1,-1,-1,0]);
boundaryv=zeros(5,1);
boundaryv(5)=1;
M=[boundary0];
for i=2:8
    if i==2
        D=diff(boundary0,m0s)*0.1;
        params=[params,':  m0s  '];
    elseif i==3
        D=diff(boundary0,T1)*T1;
        params=[params,':  T1  '];
    elseif i==4
        D=diff(boundary0,T2f)*T2f;
        params=[params,':  T2f  '];
    elseif i==5
        D=diff(boundary0,R)*R;
         params=[params,':  R  '];
    elseif i==6
        D=diff(boundary0,T2s_)*T2s_;
        params=[params,':  T2s'];
    elseif i==7
        D=diff(boundary0,B0);
        params=[params,':  B0'];
    elseif i==8
        D=diff(boundary0,B1);
        params=[params,':  B1'];
    else
        error('coding error');
    end
    M=[M;D];
end
M=[M;boundary1;boundaryv.'];
matlabFunction(M,'File','+Bloch/boundaryfun','Vars',[m0s,T1,T2f,R,T2s_,B0,B1,TR,TRFmax]);
end