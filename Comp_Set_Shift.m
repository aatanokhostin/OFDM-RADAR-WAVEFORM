function [Comp_Set]=Comp_Set_Shift(code)
N=length(code);
Comp_Set=zeros(N);
for i=0:N-1
    Comp_Set(i+1,1:end-i)=code(1,1+i:end);
    Comp_Set(i+1,end-i+1:end)=code(1,1:i);
end