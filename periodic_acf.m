function [val]=periodic_acf(x_n)
N=length(x_n);
val=zeros(1,2.*N-1);

for k=0:N-1
    x_n_1=zeros(1,N);
    x_n_1(1,1:end-k)=x_n(1,1+k:end);
    x_n_1(1,end-k+1:1:end)=x_n(1,1:k);
    val(1,N+k:end)=x_n*x_n_1'./N;
end
val(1,1:N-1)=conj(val(1,end:-1:N+1));