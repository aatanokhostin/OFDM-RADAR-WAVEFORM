function [val]=PMEPR(code)
N=length(code);
Inst_Power=abs(code).^2;
Peak=max(Inst_Power);
MEAN=mean(Inst_Power);
val=Peak./MEAN;
return