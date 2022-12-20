function [val]=OFDM_Waveform(tb,Rate,signal)
fs=1./tb;
M=size(signal,1);
Samples=size(signal,2);
t=linspace(0,Samples.*tb./Rate,Samples);
val=zeros(1,Samples);
for i=1:M
    val(1,:)=val(1,:)+exp(j.*2.*pi.*fs.*t.*((M+1)./2-i)).*signal(i,:);
end
return