function [r]=aperacfsiso_time(x)

N=length(x);
r=zeros(1,2.*N-1);

for k=0:N-1
    x_shift=zeros(1,N);
    x_shift(1,1:N-k)=x(1+k:end);
    r(k+N)=(x*x_shift')./(x*x');
end

r(1:(N-1)) = conj(r(end:-1:(N+1)));
return