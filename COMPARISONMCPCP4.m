%% Comparison Between MCPC And P4 Waveform
clc
clear all
close all

%% MCPC And P4 Codes Length
M=5; 
N=M^2;

 %% P4 Phase Code Sequence & Waveform
 P4_Phase=P4_code(N);
 
 figure();
 subplot(2,1,1); stem(mod((P4_Phase.*180./pi),360),'linewidth',2.5);
 xlabel('Samples'); ylabel('Phase In degree'); axis tight
 boldify;
 
 tc=1;
 Rate=25;
 Samples=Rate*N;
 t=linspace(0,N*tc,Samples); 
 x_t=waveform(tc,Rate,P4_Phase);
 subplot(2,1,2); plot(t,real(x_t));
 xlabel('Time In Seconds'); ylabel('Real Part Of Waveform'); axis tight
 boldify;

 %% P4 Phase Code Periodic Auto-Correlation
 P4_Seq=exp(j.*P4_Phase);
 Periodic_ACF=periodic_acf(P4_Seq);
 
 Shift=-length(P4_Seq)+1:length(P4_Seq)-1;
 
 figure();
 plot(Shift,abs(Periodic_ACF));
 xlabel('Shift'); ylabel('Magnitude Of Periodic AutoCorrelation'); axis tight
 boldify;
 
 %% Produce Complementary Set Of P4 Codes
 P4_MCPC=P4_code(M);
 Comp_Set=Comp_Set_Shift(P4_MCPC);
 Comp_Set_Deg=Comp_Set.*180./pi;
 
 figure();
 for i=1:size(Comp_Set,1)
     plot(Comp_Set_Deg(i,:)); hold on
 end
 xlabel('Phase Sample'); ylabel('Complementary Set Phases In Degree'); axis tight
 legend('Phase-seq1','Phase-seq2','Phase-seq3','Phase-seq4','Phase-seq5');
 boldify;
 
 %% Ideal Aperiodic Auto-Correlation
 MCPC_Comp_Seq=exp(j.*Comp_Set);
 MCPC_Comp_Seq_ACF=zeros(size(MCPC_Comp_Seq,1),2*size(MCPC_Comp_Seq,2)-1);
 
 for i=1:size(MCPC_Comp_Seq_ACF,1)
     MCPC_Comp_Seq_ACF(i,:)=aperacfsiso(MCPC_Comp_Seq(i,:));
 end
 
 Z_P=sum(MCPC_Comp_Seq_ACF);
 
 Comp_Shift=-size(MCPC_Comp_Seq,2)+1:size(MCPC_Comp_Seq,2)-1;
 
 figure();
 plot(Comp_Shift,abs(Z_P));
 xlabel('Shift'); ylabel('Magnitude Of Aperiodic Auto-Correlation'); axis tight
 boldify;
 
 %% P4 Phase Code Aperiodeic-Autocorrelation
 Aperiod_ACF=aperacfsiso(P4_Seq);
 
 figure();
 plot(Shift,abs(Aperiod_ACF));
 xlabel('Shift'); ylabel('Magnitude Of P4-Phase Code Aperiodic-Auto-correaltion'); axis tight
 boldify;
 
 %% P4 Phase coded Waveform PSD
 PSD_P4=abs(fft(x_t)).^2;
%  figure();
%  plot(PSD_P4);
%  xlabel('Frequency In Hertz*Sampling Rate'); ylabel('PSD Of Signal'); axis tight
%  boldify;
 
 figure();
 Frequency=linspace(0,Samples./Rate,Samples);
 plot(Frequency,PSD_P4);
 xlabel('Frequency In Hertz'); ylabel('PSD Of Signal'); axis tight
 boldify;
 
 
 %% MCPC Waveform Design
 tb=M*tc;
 fs=1./tb;
 Rate_1=100;
 Samples_1=Rate_1*M;
 Comp_Set_Waveforms=zeros(size(Comp_Set,1),Samples_1);
 for i=1:size(Comp_Set,1)
     Comp_Set_Waveforms(i,:)=waveform(tb,Rate_1,Comp_Set(i,:));
 end

 y_t=OFDM_Waveform(tb,Rate_1,Comp_Set_Waveforms);
 t_OFDM=linspace(0,M*tb,Samples_1);
 
 figure();
 plot(t_OFDM,real(y_t));
 xlabel('Time In Seconds'); ylabel('Real Part Of OFDM Waveform');
 boldify;
 
 %% OFDM Waveform Autocorrelation
 OFDM_ACF=aperacfsiso_time(y_t);
 Tau_OFDM=linspace(-max(t_OFDM)+1,max(t_OFDM)-1,2*length(t_OFDM)-1);
 
 figure();
 plot(Tau_OFDM,abs(OFDM_ACF));
 xlabel('Time in Seconds'); ylabel('Magnitude Of P4-OFDM AutoCorrelation'); axis tight
 boldify;
 
 %% OFDM With P4 Cyclic Shift Phase Code PSD
 PSD_OFDM=abs(fft(y_t)).^2;
 Frequency_1=linspace(0,M*Samples_1./Rate_1,Samples_1);
 figure();
 plot(Frequency_1,abs(PSD_OFDM));
 xlabel('Frequency In Hertz'); ylabel('PSD Of OFDM Waveform'); axis tight
 boldify;

 %% P4 Phase Code Waveform Ambiguity Function
 f_P4=linspace(-1/tc,1/tc,Samples);
 P4_d_delay=((max(t)-min(t))./(length(t)-1));
 P4_d_fd=(max(f_P4)-min(f_P4))./(length(f_P4)-1);
 P4_delay=-max(t):P4_d_delay:max(t);
 P4_fd=-max(f_P4):P4_d_fd:max(f_P4);
 P4_AMB=AF_CONV(t,f_P4,x_t);
 
 figure();
 mesh(P4_delay,P4_fd,abs(P4_AMB));
 xlabel('Delay in second'); ylabel('Doppler in Hz'); zlabel('Ambiguity Function'); axis tight
 boldify;
 
 %% P4-OFDM Ambiguity Function With Cyclic Order
 f_OFDM=linspace(-1/tb,1/tb,Samples_1);
 d_delay=(max(t_OFDM)-min(t_OFDM))./(length(t_OFDM)-1);
 d_fd=(max(f_OFDM)-min(f_OFDM))./(length(f_OFDM)-1);
 fd=min(f_OFDM):d_fd:max(f_OFDM);
 delay=-max(t_OFDM):d_delay:max(t_OFDM);
 P4_OFDM=AF_CONV(t_OFDM,f_OFDM,y_t);
 
 figure();
 mesh(delay,fd,abs(P4_OFDM));
 xlabel('Delay in second'); ylabel('Doppler in Hz'); zlabel('Ambiguity Function'); axis tight
 boldify;
  
 %% Waveform With Optimal Sequence Permutation On Subcarriers [3 5 2 1 4] For Sidelobe Reuction
Permutation=[3 5 2 1 4];
Permutation_Set_Waveforms=zeros(size(Comp_Set_Waveforms,1),size(Comp_Set_Waveforms,2));

for i=1:size(Comp_Set_Waveforms,1)
    Permutation_Set_Waveforms(i,:)=Comp_Set_Waveforms(Permutation(i),:);
end

z_t=OFDM_Waveform(tb,Rate_1,Permutation_Set_Waveforms);

figure();
plot(t_OFDM,real(z_t)); hold on
plot(t_OFDM,real(y_t),'-.');  
xlabel('Time In Seconds'); ylabel('Real Part Of waveform'); axis tight;
legend('[3 5 2 1 4]','[1 2 3 4 5]');
boldify;

%% OFDM Waveform Autocorrelation 
Permutation_OFDM_ACF=aperacfsiso_time(z_t);

figure();
plot(Tau_OFDM,abs(Permutation_OFDM_ACF));
hold on;
plot(Tau_OFDM,abs(OFDM_ACF),'-.');
xlabel('Time in Seconds'); ylabel('Magnitude OFDM ACF'); axis tight
legend('[3 5 2 1 4]','[1 2 3 4 5]');
boldify;

%% OFDM Ambiguity Function With Optimum Permutation
Permutation_OFDM_AMB=AF_CONV(t_OFDM,f_OFDM,z_t);

figure();
mesh(delay,fd,abs(Permutation_OFDM_AMB));
xlabel('Delay in second'); ylabel('Doppler in Hz'); zlabel('Ambiguity Function'); axis tight
boldify;

%% PMEPR Comparison Between Cyclic, Simple and Optimum Order
% Cyclic For PMEPR Reduction 
Cyclic_Permutation=[3 4 5 1 2];
Cyclic_Set_Waveforms=zeros(size(Comp_Set_Waveforms,1),size(Comp_Set_Waveforms,2));

for i=1:length(Cyclic_Permutation)
    Cyclic_Set_Waveforms(i,:)=Comp_Set_Waveforms(Cyclic_Permutation(i),:);
end

w_t=OFDM_Waveform(tb,Rate_1,Cyclic_Set_Waveforms);

 figure();
 plot(t_OFDM,real(z_t),t_OFDM,real(w_t));
 xlabel('Time In Seconds'); ylabel('Real Part Of Different OFDM Waveforms'); axis tight
 boldify;
 legend('Optimum SLL Distribution','Cyclic Distribution');
 
 PAPR_Optimum_SLL=PMEPR(z_t)
 PAPR_Cyclic=PMEPR(w_t)
 
 %% Train Of Complementary MCPC Of Permutation [3 5 2 1 4]
 Comp_Set_Time=Comp_Set_Shift(Permutation).';
 Comp_Waveforms_Time=zeros(size(Comp_Set_Time,1),size(Comp_Set_Time,2)*Samples_1);
 
 for i=1:size(Comp_Set_Time,1)
     for j=1:size(Comp_Set_Time,2)
         Comp_Waveforms_Time(i,(j-1).*Samples_1+1:j*Samples_1)=Comp_Set_Waveforms(Comp_Set_Time(i,j),:);
     end
 end
 
 Comp_Train_Code=OFDM_Waveform(tb,Rate_1,Comp_Waveforms_Time);
 
 Extention=2;
 t_Train=linspace(0,Extention*(M^2)*tb,Extention*size(Comp_Set_Time,2)*Samples_1);
 
 u_t=Comp_Train_Waveform(M,Extention,Comp_Train_Code);
 
 figure();
 plot(t_Train,real(u_t));
 xlabel('Time In Seconds'); ylabel('Real Part Of Complementary Train Waveform'); axis tight
 boldify;
 
 %% Train Of Complementary MCPC Auto_Correlation
 
 Comp_Train_ACF=aperacfsiso_time(u_t);
 Tau_Comp_Train=linspace(-max(t_Train)+1,max(t_Train)-1,2*length(t_Train)-1);
 
 figure();
 plot(Tau_Comp_Train,abs(Comp_Train_ACF));
 xlabel('Time in Seconds'); ylabel('Magnitude Complementary Train ACF'); axis tight
 boldify;
 
 %% Train Of Complementary MCPC Of Permutation [3 5 2 1 4] With Hamming Frequency Weighted
 a0=0.53826;
 a1=0.46164;
 alfa=0.5;
 
 Frq_Weighted_Waveform=Hamming_OFDM_Waveform(a0,a1,alfa,tb,Rate_1,Comp_Waveforms_Time);
 v_t=Comp_Train_Waveform(M,Extention,Frq_Weighted_Waveform);
 
 figure();
 plot(t_Train,real(v_t));
 xlabel('Time In Seconds'); ylabel('Real Part Of Hamming Weighted Complementary Train Waveform'); axis tight
 boldify;
 
 %% Train Of Complementary MCPC Of Permutation [3 5 2 1 4] With Hamming Frequency Weighted ACF
  
 Frq_Weighted_Waveform_ACF=aperacfsiso_time(v_t);
 Tau_Comp_Train=linspace(-max(t_Train)+1,max(t_Train)-1,2*length(t_Train)-1);
 
 figure();
 plot(Tau_Comp_Train,abs(Frq_Weighted_Waveform_ACF));
 xlabel('Time in Seconds'); ylabel('Magnitude Hamming Weighted Complementary Train ACF'); axis tight
 boldify;
 
 %% Cross-Ambiguity Function Between Cyclic And Optimum Sequence Order
 Cross=cross_AMB(t_OFDM,f_OFDM,z_t,w_t);
 
 figure();
 subplot(2,1,1); mesh(delay,fd,abs(Cross));
 xlabel('Delay in second'); ylabel('Doppler in Hz'); zlabel('Cross-Ambiguity Function'); axis tight
 boldify;
 
 subplot(2,1,2); mesh(delay,fd,abs(Permutation_OFDM_AMB));
 xlabel('Delay in second'); ylabel('Doppler in Hz'); zlabel('Ambiguity Function'); axis tight
 boldify;
 