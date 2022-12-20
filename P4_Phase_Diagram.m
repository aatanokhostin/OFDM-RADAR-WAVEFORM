%% Phase Sequence Of a P4 Code
M=input('Insert Number Of Code Length:');
m=1:M;
Phi=pi.*((m-1).^2)./M-pi.*(m-1);
plot(m,Phi*180./pi,'k','linewidth',1.5);
xlabel('mth Phase Of Sequence','fontsize',12,'fontweight','bold');
ylabel('Phase In Degree','fontsize',12,'fontweight','bold');
axis tight;
grid on;