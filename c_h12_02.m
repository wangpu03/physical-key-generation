%% simulation of the key generation between two Backscatter devices
clc;
clear;
addpath(genpath('C:\Code file\Matlab\physical-key-generation\MIToolbox-master\matlab'));

%Pt = 110;           %transmission power of the RF source
rho = 0.2;          %power split ratio for information receiver, range is [0,1]

d_f1 = 5.1;
d_f2 = 6;
%d_h12 = 3;

%fig.1 
%Pt = [100,316.2,1000,3162,10000];
log_pt = 35;
Pt = 10.^(log_pt/10);
d_h12 = 2:0.4:6;
c_real = zeros(1,length(d_h12));
c_theo = zeros(1,length(d_h12));
e_real = zeros(1,length(d_h12));
e_theo = zeros(1,length(d_h12));

for i = 1:length(d_h12)
    [c_real(i),c_theo(i),e_real(i), e_theo(i)] = key_gen_func_02(d_f1,d_f2,d_h12(i),rho, Pt);
end

figure(1);
yyaxis left
plot(d_h12,c_real,'-o','DisplayName','Real MI','LineWidth',1);
hold on;
plot(d_h12,c_theo,'-*','DisplayName','Theo MI','LineWidth',1);
%plot(Pt,c_sk);
ylim([0 3.3])
yyaxis right
plot(d_h12,e_real,'-->','DisplayName','Real EH','LineWidth',1);
plot(d_h12,e_theo,'--<','DisplayName','Theo EH','LineWidth',1);
ylim([50 155])

hold off;
grid on;
xlabel('the distance between two BDs');
yyaxis left
ylabel ('Mutual Informaiton');
yyaxis right
ylabel ('Energy harvesting');
legend;

