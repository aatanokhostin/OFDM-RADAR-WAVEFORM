function [AMB]=cross_AMB(t,f,x1,x2)
%% fd Calculation
dfd=(max(f)-min(f))./(length(f)-1);
fd=min(f):dfd:max(f);

%% Convolution Matrices Calculation

x1_t=zeros(1,length(x1));
x1_t(1,1:end)=conj(x1(end:-1:1));

x2_t=zeros(length(fd),length(x2));
for i=1:length(fd)
    x2_t(i,:)=x2.*exp(j.*2.*pi.*fd(i).*t);
end

%% Ambiguity Function Calculation

AMB=zeros(length(fd),2.*length(x1)-1);
for i=1:length(fd)
    AMB(i,:)=conv(x2_t(i,:),x1_t);
end

return