function [r ISL PSL]=aperacfsiso(x)

N=length(x);
r=zeros(1,2.*N-1);

for k=0:N-1
    x_shift=zeros(1,N);
    x_shift(1,1:N-k)=x(1+k:end);
    r(k+N)=(x*x_shift')./(x*x');
end

r(1:(N-1)) = conj(r(end:-1:(N+1)));

if nargout>=2
    ISL=sum(abs(r).^2)-(x*x').^2;
    if nargout==3
        PSL=(max(abs(r(1:N-1))).^2)./((x*x').^2);
    end
end
return