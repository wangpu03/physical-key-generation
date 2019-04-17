%% simulation of the key generation between two Backscatter devices
clc;
clear;
addpath(genpath('C:\Code file\Matlab\physical-key-generation\MIToolbox-master\matlab'));

%Pt = 110;           %transmission power of the RF source
rho = 0.2;          %power split ratio for information receiver, range is [0,1]

d_f1 = 5.1;
d_f2 = 6;
d_h12 = 3;

%fig.1 
%Pt = [100,316.2,1000,3162,10000];
log_pt = 25:2.5:50;
Pt = 10.^(log_pt/10);
c_real = zeros(1,length(Pt));
c_theo = zeros(1,length(Pt));
e_real = zeros(1,length(Pt));
e_theo = zeros(1,length(Pt));
for i = 1:length(Pt)
    [c_real(i),c_theo(i),e_real(i), e_theo(i)] = key_gen_func_02(d_f1,d_f2,d_h12,rho, Pt(i));
end

figure(1);
yyaxis left
plot(log_pt,c_real,'-o','DisplayName','Real MI','LineWidth',1);
hold on;
plot(log_pt,c_theo,'-*','DisplayName','Theo MI','LineWidth',1);
%plot(Pt,c_sk);
yyaxis right
plot(log_pt,e_real,'-->','DisplayName','Real EH','LineWidth',1);
plot(log_pt,e_theo,'--<','DisplayName','Theo EH','LineWidth',1);
%xlim([25 50])
hold off;
grid on;
xlabel('SNR(dB)');
yyaxis left
ylabel ('Mutual Informaiton');
yyaxis right
ylabel ('Energy harvesting');
legend;


% figure(2);
% yyaxis left
% plot(Pt,c_real,'-o','DisplayName','Real MI','LineWidth',1);
% hold on;
% plot(Pt,c_theo,'-*','DisplayName','Theo MI','LineWidth',1);
% %plot(Pt,c_sk);
% yyaxis right
% plot(Pt,e_real,'-->','DisplayName','Real EH','LineWidth',1);
% plot(Pt,e_theo,'--<','DisplayName','Theo EH','LineWidth',1);
% hold off;
% grid on;
% xlabel('SNR(dB)');
% yyaxis left
% ylabel ('Mutual Informaiton');
% yyaxis right
% ylabel ('Energy harvesting');
% 
% legend;
