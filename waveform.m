function [val]=waveform(tc,Rate,code)
N=length(code);
Samples=Rate*N;
t=linspace(0,N*tc,Samples);
val=zeros(1,Samples);
for i=0:N-1
    val=val+(stepfun(t,(i*tc))-stepfun(t,(i+1).*tc)).*exp(j.*code(i+1));
end
return
