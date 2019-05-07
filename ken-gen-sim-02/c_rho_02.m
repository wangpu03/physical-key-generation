%% simulation of the key generation between two Backscatter devices
clc;
clear;
addpath(genpath('C:\Code file\Matlab\physical-key-generation\MIToolbox-master\matlab'));

log_pt = 35;
Pt = 10.^(log_pt/10);                %transmission power of the RF source
%rho = 0.2;             %power split ratio for information receiver, range is [0,1]

d_f1 = 5.1;
d_f2 = 6;
d_h12 = 3;

rho = 0:0.1:1;
c_real = zeros(1,length(rho));
c_theo = zeros(1,length(rho));
e_real = zeros(1,length(rho));
e_theo = zeros(1,length(rho));

for i = 1:length(rho)
    [c_real(i),c_theo(i),e_real(i), e_theo(i)] = key_gen_func_02(d_f1,d_f2,d_h12,rho(i), Pt);
end

figure(1);
yyaxis left
plot(rho,c_real,'-o','DisplayName','Real MI','LineWidth',1);
hold on;
plot(rho,c_theo,'-*','DisplayName','Theo MI','LineWidth',1);
%plot(Pt,c_sk);
yyaxis right
plot(rho,e_real,'-->','DisplayName','Real EH','LineWidth',1);
plot(rho,e_theo,'--<','DisplayName','Theo EH','LineWidth',1);
%xlim([25 50])
hold off;
grid on;
xlabel('rho');
yyaxis left
ylabel ('Mutual Informaiton');
yyaxis right
ylabel ('Energy harvesting');
legend;

%% need to run it gain and plot with different SNR

lambda1 = 0.7;
G1 = (1-lambda1)*(c_real)/max(c_real)+lambda1*(e_real/max(e_real));
lambda3 = 0.4;
G3 = (1-lambda3)*(c_real)/max(c_real)+lambda3*(e_real/max(e_real));
lambda5 = 0.1;
G5 = (1-lambda5)*(c_real)/max(c_real)+lambda5*(e_real/max(e_real));

figure(2);
plot(rho,G1,'-o','DisplayName','lambda = 0.7','LineWidth',1);
hold on;
plot(rho,G3,'-*','DisplayName','lambda = 0.5','LineWidth',1);
plot(rho,G5,'-d','DisplayName','lambda = 0.3','LineWidth',1);
grid on;
hold off;
title('Pt = 25')
xlabel('rho');
ylabel ('Mutual Informaiton');
legend;