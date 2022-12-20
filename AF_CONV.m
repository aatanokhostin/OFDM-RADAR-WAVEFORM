function [AMB]=AF_CONV(t,f,y_t)
%% fd Calculation
dfd=(max(f)-min(f))./(length(f)-1);
fd=min(f):dfd:max(f);

%% Convolution Matrices Calculation

y1_t=zeros(1,length(y_t));
y1_t(1,1:end)=conj(y_t(end:-1:1));

y2_t=zeros(length(fd),length(y_t));
for i=1:length(fd)
    y2_t(i,:)=y_t.*exp(j.*2.*pi.*fd(i).*t);
end

%% Ambiguity Function Calculation

AMB=zeros(length(fd),2.*length(y_t)-1);
for i=1:length(fd)
    AMB(i,:)=conv(y2_t(i,:),y1_t);
end