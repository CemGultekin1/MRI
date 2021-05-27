function S=spline_basis(grid,nedge,type)
grid=grid(:);
mg=min(grid);
Mg=max(grid);
grid=(grid-mg)/(Mg-mg);
sedge=linspace(0,1,nedge).';
ngrid=length(grid);
switch type
    case 'C0-spline'
        mm=3;
        S=zeros(ngrid,(nedge-1)*mm+1);
        for i=2:nedge-1
            I=grid>=sedge(i-1) & grid<=sedge(i+1);
            xx=(grid(I)-sedge(i-1))/(sedge(i+1)-sedge(i-1))*2-1;
            yy=(1-abs(xx));
            S(I,(i-1)*mm+1)=yy;
        end
        I=grid>=sedge(1) & grid<=sedge(2);
        xx=(grid(I)-sedge(1))/(sedge(2)-sedge(1));
        yy=(1-abs(xx));
        S(I,1)=yy;

        I=grid>=sedge(end-1) & grid<=sedge(end);
        xx=(grid(I)-sedge(end-1))/(sedge(end)-sedge(end-1))-1;
        yy=(1-abs(xx));
        S(I,end)=yy;
        for i=1:nedge-1
            I=grid>=sedge(i) & grid<=sedge(i+1);
            xx=(grid(I)-sedge(i))/(sedge(i+1)-sedge(i));
            yy=(1-xx).*2.*xx;
            for j=2:mm
                ct=chebyshevT(j-2,2*xx-1);
                S(I,(i-1)*mm+j)=yy.*ct;
            end
        end
    case 'hann'
        mm=1;
        S=zeros(ngrid,(nedge-1)*mm+1);
        for i=2:nedge-1
            I=grid>=sedge(i-1) & grid<=sedge(i+1);
            xx=(grid(I)-sedge(i-1))/(sedge(i+1)-sedge(i-1))*2-1;
            yy=cos(pi/2*xx).^2;
            S(I,(i-1)*mm+1)=yy;
        end
        I=grid>=sedge(1) & grid<=sedge(2);
        xx=(grid(I)-sedge(1))/(sedge(2)-sedge(1));
        yy=cos(pi/2*xx).^2;
        S(I,1)=yy;

        I=grid>=sedge(end-1) & grid<=sedge(end);
        xx=(grid(I)-sedge(end-1))/(sedge(end)-sedge(end-1))-1;
        yy=cos(pi/2*xx).^2;
        S(I,end)=yy;
end
%{
S=zeros(ngrid,nedge*2);
for i=2:nedge-1
    I=grid>=sedge(i-1) & grid<=sedge(i+1);
    xx=(grid(I)-sedge(i-1))/(sedge(i+1)-sedge(i-1))*2-1;
    yy=(xx+1).^2.*(xx-1).^2.*[(xx.^2+1),xx];
    S(I,(i-1)*2+(1:2))=yy;
end
I=grid>=sedge(1) & grid<=sedge(2);
xx=(grid(I)-sedge(1))/(sedge(2)-sedge(1));
yy=(xx+1).^2.*(xx-1).^2.*[(xx.^2+1),xx];
S(I,1:2)=yy;

I=grid>=sedge(end-1) & grid<=sedge(end);
xx=(grid(I)-sedge(end-1))/(sedge(end)-sedge(end-1))-1;
yy=(xx+1).^2.*(xx-1).^2.*[(xx.^2+1),xx];
S(I,end-2+(1:2))=yy;
%}

%{
Y=randn(ngrid,1);
Y=exp(-grid);
coeff=S\Y;
plot(grid,Y);hold on;plot(grid,S*coeff);
%}
end