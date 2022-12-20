function [val]=Comp_Train_Waveform(M,Extention,code)

Samples=Extention*size(code,2);
New_Arrange=zeros(M,size(code,2)./M);

for i=1:M
    New_Arrange(i,:)=code(1,1+(i-1)*size(code,2)./M:i*size(code,2)./M);
end

val=zeros(1,Samples);

for j=1:M
    val(1,1+(j-1)*Extention*size(code,2)./M:(size(code,2)./M)+(j-1)*Extention*size(code,2)./M)=New_Arrange(j,:);
end

return