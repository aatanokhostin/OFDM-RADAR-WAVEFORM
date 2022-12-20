function [val]=P4_code(N)
n=1:N;
val=pi.*((n-1).^2)./N-pi.*(n-1);
return