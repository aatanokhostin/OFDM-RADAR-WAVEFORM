%% Complementary Set Of Ideal Periodic Autocorrelation Phase Sequence
clc
close all
clear all
N=5;
Int_Seq=P4_code(N);
Comp_Set=zeros(N);
for i=0:N-1
    Comp_Set(i+1,1:end-i)=Int_Seq(1,1+i:end);
    Comp_Set(i+1,end-i+1:end)=Int_Seq(1,1:i);
end
Result=fopen('Complementary Set.txt','w');
Comp_Set_Deg=Comp_Set.*180./pi;
for i=1:N
    fprintf(Result,'%d\r',Comp_Set_Deg(i,:));
    fprintf(Result,'\n');
end
fclose(Result);
xlswrite('exportdata.xlsx',Comp_Set_Deg);