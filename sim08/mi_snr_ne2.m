clc;
clear;
datestr(now)
addpath('E:\Github\physical-key-generation\MIToolbox-master\matlab');

% modulation methods: BPSK, QPSK,16QAM, 32QAM,64QAM
mod_method = 'QPSK';
% calculate modulation order from modulation method
mod_methods = {'BPSK', 'QPSK','8PSK','16QAM', '32QAM','64QAM'};
mod_order = find(ismember(mod_methods, mod_method));

% 该数据可以生成一次之后可以永久使用，保证每次调制的数据一致
% rand_ints_gen = randi(2,10000,1)-1;
% save data_input.txt -ascii rand_ints_gen
rand_ints = load("data_input.txt");

num_mc = 5000;

%
% nfft = 256, K =1
n_fft = 256;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=1;
c_flag = 1;

% generate the data 
rand_ints_1256 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_1256, mod_method, n_fft, n_cp, c_flag);



%% snr=0
cor = [0.9999    0.9991    0.9964    0.9857    0.9118    0.6701    0.0001, 0, 0,0];
v1 = zeros(num_mc, length(cor));
v2 = zeros(num_mc, length(cor));
snr = 0;
for i = 1:length(cor)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr_ne(data_ofdm, n_ofdm,n_cp, K, snr, cor(i));
    [v1_trad(j,i),v2_trad(j,i)] = mi_snr_ne_trad2(data_ofdm, n_ofdm,n_cp, K, snr, cor(i));
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);
v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);
for i = 1:length(cor)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

mean_v1_trad = repmat(mean(v1_trad), num_mc,1);
mean_v2_trad = repmat(mean(v2_trad), num_mc,1);
v1_bit_trad = double(v1_trad>= mean_v1_trad);
v2_bit_trad = double(v2_trad>= mean_v2_trad);
for i = 1:length(cor)
    mi_bit_trad(i) = mi(v1_bit_trad(:,i),v2_bit_trad(:,i));
end

v1_snr0 = v1;
v2_snr0 = v2;
v1_bit_snr0 = v1_bit;
v2_bit_snr0 = v2_bit;
mi_bit_snr0 = mi_bit;

v1_snr0_trad = v1_trad;
v2_snr0_trad = v2_trad;
v1_bit_snr0_trad = v1_bit_trad;
v2_bit_snr0_trad = v2_bit_trad;
mi_bit_snr0_trad = mi_bit_trad;

figure(1);
plot(mi_bit_snr0,'r-x','LineWidth',1.5);
hold on;
plot(mi_bit_snr0_trad,'r--x','LineWidth',1.5);


%% snr=30
cor = [0.9999    0.9991    0.9964    0.9857    0.9118    0.6701    0.0001, 0, 0, 0];
v1 = zeros(num_mc, length(cor));
v2 = zeros(num_mc, length(cor));
snr = 30;
for i = 1:length(cor)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr_ne(data_ofdm, n_ofdm,n_cp, K, snr, cor(i));
    [v1_trad(j,i),v2_trad(j,i)] = mi_snr_ne_trad2(data_ofdm, n_ofdm,n_cp, K, snr, cor(i));
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);
v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);
for i = 1:length(cor)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

mean_v1_trad = repmat(mean(v1_trad), num_mc,1);
mean_v2_trad = repmat(mean(v2_trad), num_mc,1);
v1_bit_trad = double(v1_trad>= mean_v1_trad);
v2_bit_trad = double(v2_trad>= mean_v2_trad);
for i = 1:length(cor)
    mi_bit_trad(i) = mi(v1_bit_trad(:,i),v2_bit_trad(:,i));
end

v1_snr30 = v1;
v2_snr30 = v2;
v1_bit_snr30 = v1_bit;
v2_bit_snr30 = v2_bit;
mi_bit_snr30 = mi_bit;

v1_snr30_trad = v1_trad;
v2_snr30_trad = v2_trad;
v1_bit_snr30_trad = v1_bit_trad;
v2_bit_snr30_trad = v2_bit_trad;
mi_bit_snr30_trad = mi_bit_trad;


%% snr=0
cor2 = [1.0000    1.0000    1.0000    0.9999    0.9991    0.9964    0.9856    0.9116    0.6699   0.0002];
v1 = zeros(num_mc, length(cor));
v2 = zeros(num_mc, length(cor));
snr = 0;
for i = 1:length(cor)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr_ne(data_ofdm, n_ofdm,n_cp, K, snr, cor2(i));
    [v1_trad(j,i),v2_trad(j,i)] = mi_snr_ne_trad2(data_ofdm, n_ofdm,n_cp, K, snr, cor2(i));
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);
v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);
for i = 1:length(cor)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

mean_v1_trad = repmat(mean(v1_trad), num_mc,1);
mean_v2_trad = repmat(mean(v2_trad), num_mc,1);
v1_bit_trad = double(v1_trad>= mean_v1_trad);
v2_bit_trad = double(v2_trad>= mean_v2_trad);
for i = 1:length(cor)
    mi_bit_trad(i) = mi(v1_bit_trad(:,i),v2_bit_trad(:,i));
end

v1_snr0_cor2 = v1;
v2_snr0_cor2 = v2;
v1_bit_snr0_cor2 = v1_bit;
v2_bit_snr0_cor2 = v2_bit;
mi_bit_snr0_cor2 = mi_bit;

v1_snr0_trad_cor2 = v1_trad;
v2_snr0_trad_cor2 = v2_trad;
v1_bit_snr0_trad_cor2 = v1_bit_trad;
v2_bit_snr0_trad_cor2 = v2_bit_trad;
mi_bit_snr0_trad_cor2 = mi_bit_trad;


%% snr=30
cor2 = [1.0000    1.0000    1.0000    0.9999    0.9991    0.9964    0.9856    0.9116    0.6699   0.0002];
v1 = zeros(num_mc, length(cor));
v2 = zeros(num_mc, length(cor));
snr = 30;
for i = 1:length(cor)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr_ne(data_ofdm, n_ofdm,n_cp, K, snr, cor2(i));
    [v1_trad(j,i),v2_trad(j,i)] = mi_snr_ne_trad2(data_ofdm, n_ofdm,n_cp, K, snr, cor2(i));
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);
v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);
for i = 1:length(cor)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

mean_v1_trad = repmat(mean(v1_trad), num_mc,1);
mean_v2_trad = repmat(mean(v2_trad), num_mc,1);
v1_bit_trad = double(v1_trad>= mean_v1_trad);
v2_bit_trad = double(v2_trad>= mean_v2_trad);
for i = 1:length(cor)
    mi_bit_trad(i) = mi(v1_bit_trad(:,i),v2_bit_trad(:,i));
end

v1_snr30 = v1;
v2_snr30_cor2 = v2;
v1_bit_snr30_cor2 = v1_bit;
v2_bit_snr30_cor2 = v2_bit;
mi_bit_snr30_cor2 = mi_bit;

v1_snr30_trad_cor2 = v1_trad;
v2_snr30_trad_cor2 = v2_trad;
v1_bit_snr30_trad_cor2 = v1_bit_trad;
v2_bit_snr30_trad_cor2 = v2_bit_trad;
mi_bit_snr30_trad_cor2 = mi_bit_trad;




time = [10^-4,2.5*10^-4,5*10^-4,10^-3,2.5 *10^-3,5 *10^-3, 10^-2, 2.5*10^-2, 5*10^-2,10^-1];
%num = floor([10^-5,2.5*10^-5,5*10^-5,10^-4,2.5 *10^-4,5 *10^-4,10^-3,2.5*10^-3,5*10^-3,10^-2]*10^6*3.828);
% plot(mi_cp_bit_cor1,'-d','LineWidth',1.5,'Color', '#77AC30');
% axis([-0.1 0.95 0.0 0.06]);

% 
figure(1);
plot(mi_bit_snr0,'r-s','LineWidth',1.5);
hold on;
plot(mi_bit_snr0_trad,'r--s','LineWidth',1.5);
plot(mi_bit_snr30,'k-o','LineWidth',1.5);
plot(mi_bit_snr30_trad,'k--o','LineWidth',1.5);
hold off;
grid on;
axis([0.5 10.5 -0.1 1]);
%time = [10^-4,2.5*10^-4,5*10^-4,10^-3,2.5 *10^-3,5 *10^-3, 10^-2, 2.5*10^-2, 5*10^-2,10^-1];
set(gca,'XTickLabel',{'0.0001','0.00025', '0.0005', '0.001', '0.0025', '0.005', '0.01', '0.025', '0.05','0.1'});
xtickangle(50)
legend('snr = 0,    三角信道','snr = 0,    内向信道','snr = 30,  三角信道','snr = 30,  内向信道','Fontname','<宋体>');
ylabel('每比特互信息（bits）','Fontname','<宋体>');
xlabel('测量间隔 (s)','Fontname','<宋体>');


%%
figure(2);
plot(mi_bit_snr0_cor2,'r-s','LineWidth',1.5);
hold on;
plot(mi_bit_snr0_trad_cor2,'r--s','LineWidth',1.5);
plot(mi_bit_snr30_cor2,'k-o','LineWidth',1.5);
plot(mi_bit_snr30_trad_cor2,'k--o','LineWidth',1.5);
hold off;
grid on;
axis([0.5 10.5 -0.1 1]);
set(gca,'XTickLabel',{'0.0001','0.00025', '0.0005', '0.001', '0.0025', '0.005', '0.01', '0.025', '0.05','0.1'});
xtickangle(50)
legend('snr = 0,    三角信道','snr = 0,    内向信道','snr = 30,  三角信道','snr = 30,  内向信道','Fontname','<宋体>');
ylabel('每比特互信息（bits）','Fontname','<宋体>');
xlabel('测量间隔 (s)','Fontname','<宋体>');

figure(3);
plot(mi_bit_snr0_cor2,'k-s','LineWidth',1.5);
hold on;
plot(mi_bit_snr0_trad_cor2,'k--s','LineWidth',1.5);
plot(mi_bit_snr0,'r-s','LineWidth',1.5);
plot(mi_bit_snr0_trad,'r--s','LineWidth',1.5);
hold off;
grid on;
axis([0.5 10.5 -0.1 1]);
set(gca,'XTickLabel',{'0.0001','0.00025', '0.0005', '0.001', '0.0025', '0.005', '0.01', '0.025', '0.05','0.1'});
xtickangle(50)
ylabel('Mutual Information per bit');
xlabel('Time delay (s)');
legend('fd = 10Hz,    triangle','fd = 10Hz,    inward','fd = 100Hz,  triangle','fd = 100Hz,  inward');

