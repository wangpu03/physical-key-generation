% plot_mi_eth.m

clc;
clear;

d = [10 7 3];
alpha = 0.3 + 1i*0.4;

lamda = 3;
h = 5*d.^(-lamda/2);

rho1 =0:0.01:1;
rho2 =0:0.01:1;

%%
snr = 30;
pt = 1;
P = 4*pt*abs(alpha)*h(1)*h(2)*h(3);
pw_noise = pt/10^(snr/10);

E_th = 0:0.001:0.1;
mi_max = zeros(1,length(E_th));
mi_diag_max = zeros(1,length(E_th));
for j = 1:length(E_th)
    mi = zeros(length(rho1),length(rho2));
    for i = 1:length(rho1)
        E1 = 0.5 * pt * (1-rho1(i))* h(1)^2;
        if E1 <  E_th(j)
            break;
        end
    end
    num_rho1 = i-1;
    
    for i = 1:length(rho2)
        E2 = 0.5 * pt * (1-rho2(i))* h(2)^2;
        if E2 <  E_th(j)
            break;
        end
    end
    num_rho2 = i-1;
    
    for i = 1:num_rho1
        for k = 1:num_rho2
            mi(i,k) = 1/2*log10(1+rho1(i)^2*rho2(k)^2*P/(rho1(i)^2*pw_noise+rho2(k)^2*pw_noise+pw_noise*pw_noise/P));
        end
    end
    
    mi_max (j) = max(max(mi)');
    mi_diag_max(j) = max(diag(mi));
end


%%
snr = 25;
pt = 1;
P = 4*pt*abs(alpha)*h(1)*h(2)*h(3);
pw_noise = pt/10^(snr/10);

E_th = 0:0.001:0.1;
mi_25_max = zeros(1,length(E_th));
mi_diag_25_max = zeros(1,length(E_th));
for j = 1:length(E_th)
    mi = zeros(length(rho1),length(rho2));
    for i = 1:length(rho1)
        E1 = 0.5 * pt * (1-rho1(i))* h(1)^2;
        if E1 <  E_th(j)
            break;
        end
    end
    num_rho1 = i-1;
    
    for i = 1:length(rho2)
        E2 = 0.5 * pt * (1-rho2(i))* h(2)^2;
        if E2 <  E_th(j)
            break;
        end
    end
    num_rho2 = i-1;
    
    for i = 1:num_rho1
        for k = 1:num_rho2
            mi(i,k) = 1/2*log10(1+rho1(i)^2*rho2(k)^2*P/(rho1(i)^2*pw_noise+rho2(k)^2*pw_noise+pw_noise*pw_noise/P));
        end
    end
    
    mi_25_max (j) = max(max(mi)');
    %mi_diag_35_max(j) = max(diag(mi));
end

%%
snr = 30;
pt = 1.5;
P = 4*pt*abs(alpha)*h(1)*h(2)*h(3);
pw_noise = 1/10^(snr/10);

E_th = 0:0.001:0.1;
mi_2_max = zeros(1,length(E_th));
mi_diag_2_max = zeros(1,length(E_th));
for j = 1:length(E_th)
    mi = zeros(length(rho1),length(rho2));
    for i = 1:length(rho1)
        E1 = 0.5 * pt * (1-rho1(i))* h(1)^2;
        if E1 <  E_th(j)
            break;
        end
    end
    num_rho1 = i-1;
    
    for i = 1:length(rho2)
        E2 = 0.5 * pt * (1-rho2(i))* h(2)^2;
        if E2 <  E_th(j)
            break;
        end
    end
    num_rho2 = i-1;
    
    for i = 1:num_rho1
        for k = 1:num_rho2
            mi(i,k) = 1/2*log10(1+rho1(i)^2*rho2(k)^2*P/(rho1(i)^2*pw_noise+rho2(k)^2*pw_noise+pw_noise*pw_noise/P));
        end
    end
    
    mi_2_max (j) = max(max(mi)');
    %mi_diag_35_max(j) = max(diag(mi));
end

num_max = sum(mi_max ~= 0);
num_diag_max = sum(mi_diag_max ~= 0);
num_25_max = sum(mi_25_max ~= 0);
num_2_max = sum(mi_2_max ~= 0);
% 
plot([E_th(1:num_max),E_th(num_max)], [mi_max(1:num_max),0],'r-*','LineWidth',1);
hold on;
plot(E_th(1:num_diag_max), mi_diag_max(1:num_diag_max),'b-s','LineWidth',1);
hold on;
plot(E_th(1:num_25_max), mi_25_max(1:num_25_max),'k-v','LineWidth',1);
hold on;
plot([E_th(1:num_2_max),E_th(num_2_max)], [mi_2_max(1:num_2_max),0],'k-o','LineWidth',1);
hold on;
% % plot(E_th, mi_diag_max_2,'LineWidth',1.5);

area1 = area([0,E_th(4)],[mi_max(4),mi_max(4)], 'LineStyle','none');
plot([0,E_th(4)],[mi_max(4),mi_max(4)],'LineStyle','--','LineWidth',1.5,'Color','#0072BD');
plot([E_th(4),E_th(4)],[0,mi_max(4)],'LineStyle','--','LineWidth',1.5,'Color','#0072BD');

area2 = area([0,E_th(9)],[mi_max(9),mi_max(9)], 'LineStyle','none');
plot([0,E_th(9)],[mi_max(9),mi_max(9)],'LineStyle','--','LineWidth',1.5,'Color','#77AC30');
plot([E_th(9),E_th(9)],[0,mi_max(9)],'LineStyle','--','LineWidth',1.5,'Color','#77AC30');

area3 = area([0,E_th(12)],[mi_max(12),mi_max(12)], 'LineStyle','none');
plot([0,E_th(12)],[mi_max(12),mi_max(12)],'LineStyle','--','LineWidth',1.5,'Color','#A2142F');
plot([E_th(12),E_th(12)],[0,mi_max(12)],'LineStyle','--','LineWidth',1.5,'Color','#A2142F');

area1.FaceColor = '#0072BD';
area2.FaceColor = '#77AC30';
area3.FaceColor = '#A2142F';
area1.FaceAlpha = 0.1;
area2.FaceAlpha = 0.1;
area3.FaceAlpha = 0.1;
grid on;

dim = [0.003/0.0127 0.725/1.05 0.0500 0.0500];
a1 = annotation('textbox',dim,'String','K1');
a1.Color = 'red';
a1.FontSize = 12;
a1.LineStyle = 'none';

dim = [0.003/0.0068 0.725/1.476 0.0500 0.0500];
a2 = annotation('textbox',dim,'String','K2');
a2.Color = 'red';
a2.FontSize = 12;
a2.LineStyle = 'none';

dim = [0.003/0.0054 0.725/3.1 0.0500 0.0500];
a3 = annotation('textbox',dim,'String','K3');
a3.Color = 'red';
a3.FontSize = 12;
a3.LineStyle = 'none';


axis([0 0.02 0 1.05]);
ylabel('每比特互信息（bits）','Fontname','<宋体>');
xlabel('能量收集门限','Fontname','<宋体>');
legend('信号功率 = 1, 噪声功率 =0.001','信号功率 = 1, 噪声功率 =0.001 (\rho_1 = \rho_2)','信号功率 = 1, 噪声功率 =0.003','信号功率 = 1.5, 噪声功率 =0.001','Fontname','<宋体>')

%mesh(rho1,rho2,mi);





