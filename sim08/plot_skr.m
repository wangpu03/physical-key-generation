clc;
clear;

addpath('E:\Github\physical-key-generation\MIToolbox-master\matlab');

% modulation methods: BPSK, QPSK,16QAM, 32QAM,64QAM
mod_method = 'QPSK';
% calculate modulation order from modulation method
mod_methods = {'BPSK', 'QPSK','8PSK','16QAM', '32QAM','64QAM'};
mod_order = find(ismember(mod_methods, mod_method));

%
% nfft = 256, K =1
n_fft = 256;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K = 1;
c_flag = 1;
% generate the data 
% rand_ints_gen = randi(2,n_fft*mod_order*K,1)-1;
% save data_input.txt -ascii rand_ints_gen
rand_ints = load("data_input.txt");
rand_ints_1256 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_1256, mod_method, n_fft, n_cp, c_flag);

num_sim = 50000;
snr = 0:3:30;
v1 = zeros(num_sim, length(snr));
v2 = zeros(num_sim, length(snr));
ve1 = zeros(num_sim, length(snr));
ve2 = zeros(num_sim, length(snr));

%%
cor = 0;
d_e1 = 1;  % distance between device 1 and eavesdropper
%[v_y1_cp, v_y2_cp, v_ye1_cp, v_ye2_cp] = skr_snr(data_ofdm, n_ofdm, n_cp ,n_frame, cor, d_e1, snr)
for i = 1:length(snr)
for j = 1:num_sim
    [v1(j,i),v2(j,i),ve1(j,i),ve2(j,i)] = skr_snr(data_ofdm, n_ofdm, n_cp ,K, cor, d_e1, snr(i));
end
end

mean_v1 = repmat(mean(v1), num_sim,1);
mean_v2 = repmat(mean(v2), num_sim,1);
mean_ve1 = repmat(mean(ve1), num_sim,1);
mean_ve2 = repmat(mean(ve2), num_sim,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);
ve1_bit = double(ve1>= mean_ve1);
ve2_bit = double(ve2>= mean_ve2);

cmi_ve2 = zeros(1,length(snr));   %%发射信号和直射信号相乘，得到第一类窃听方法
cmi_ve1 = zeros(1,length(snr));   %%发射信号和直射信号相乘，得到第二类窃听方法
mi_bit = zeros(1,length(snr));
for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

for i = 1:length(snr)
    cmi_ve2(i) = cmi(v1_bit(:,i),v2_bit(:,i),ve2_bit(:,i));
end

for i = 1:length(snr)
    cmi_ve1(i) = cmi(v1_bit(:,i),v2_bit(:,i),ve1_bit(:,i));
end

v1_0 = v1;
v2_0 = v2;
ve1_0 = ve1;
ve2_0 = ve2;
mi_bit_0 = mi_bit;
cmi_ve2_0 = cmi_ve2;
cmi_ve1_0 = cmi_ve1;

%%
cor = 0.6;
d_e1 = 1;
%[v_y1_cp, v_y2_cp, v_ye1_cp, v_ye2_cp] = skr_snr(data_ofdm, n_ofdm, n_cp ,n_frame, cor, d_e1, snr)
for i = 1:length(snr)
for j = 1:num_sim
    [v1(j,i),v2(j,i),ve1(j,i),ve2(j,i)] = skr_snr(data_ofdm, n_ofdm, n_cp ,K, cor, d_e1, snr(i));
end
end

mean_v1 = repmat(mean(v1), num_sim,1);
mean_v2 = repmat(mean(v2), num_sim,1);
mean_ve1 = repmat(mean(ve1), num_sim,1);
mean_ve2 = repmat(mean(ve2), num_sim,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);
ve1_bit = double(ve1>= mean_ve1);
ve2_bit = double(ve2>= mean_ve2);

cmi_ve2 = zeros(1,length(snr));
cmi_ve1 = zeros(1,length(snr));
mi_bit = zeros(1,length(snr));
for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end
for i = 1:length(snr)
    cmi_ve2(i) = cmi(v1_bit(:,i),v2_bit(:,i),ve2_bit(:,i));
end

for i = 1:length(snr)
    cmi_ve1(i) = cmi(v1_bit(:,i),v2_bit(:,i),ve1_bit(:,i));
end

v1_6 = v1;
v2_6 = v2;
ve1_6 = ve1;
ve2_6 = ve2;
mi_bit_6 = mi_bit;
cmi_ve2_6 = cmi_ve2;
cmi_ve1_6 = cmi_ve1;


%%  在窃听时，只考虑设备之间信道的随机性，其他信道保持一致，作为传统信道生成模型
cor = 0;
d_e1 = 1;
%[v_y1_cp, v_y2_cp, v_ye1_cp, v_ye2_cp] = skr_snr(data_ofdm, n_ofdm, n_cp ,n_frame, cor, d_e1, snr)
for i = 1:length(snr)
for j = 1:num_sim
    [v1(j,i),v2(j,i),ve1(j,i),ve2(j,i)] = skr_snr_trad(data_ofdm, n_ofdm, n_cp ,K, cor, d_e1, snr(i));
end
end

mean_v1 = repmat(mean(v1), num_sim,1);
mean_v2 = repmat(mean(v2), num_sim,1);
mean_ve1 = repmat(mean(ve1), num_sim,1);
mean_ve2 = repmat(mean(ve2), num_sim,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);
ve1_bit = double(ve1>= mean_ve1);
ve2_bit = double(ve2>= mean_ve2);

cmi_ve2 = zeros(1,length(snr));
cmi_ve1 = zeros(1,length(snr));
mi_bit = zeros(1,length(snr));
for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end
for i = 1:length(snr)
    cmi_ve2(i) = cmi(v1_bit(:,i),v2_bit(:,i),ve2_bit(:,i));
end

for i = 1:length(snr)
    cmi_ve1(i) = cmi(v1_bit(:,i),v2_bit(:,i),ve1_bit(:,i));
end

v1_trad_0 = v1;
v2_trad_0 = v2;
ve1_trad_0 = ve1;
ve2_trad_0 = ve2;
mi_trad_bit_0 = mi_bit;
cmi_trad_ve2_0 = cmi_ve2;
cmi_trad_ve1_0 = cmi_ve1;

%%  在窃听时，只考虑设备之间信道的随机性，其他信道保持一致
cor = 0.6;
d_e1 = 1;
%[v_y1_cp, v_y2_cp, v_ye1_cp, v_ye2_cp] = skr_snr(data_ofdm, n_ofdm, n_cp ,n_frame, cor, d_e1, snr)
for i = 1:length(snr)
for j = 1:num_sim
    [v1(j,i),v2(j,i),ve1(j,i),ve2(j,i)] = skr_snr_trad(data_ofdm, n_ofdm, n_cp ,K, cor, d_e1, snr(i));
end
end

mean_v1 = repmat(mean(v1), num_sim,1);
mean_v2 = repmat(mean(v2), num_sim,1);
mean_ve1 = repmat(mean(ve1), num_sim,1);
mean_ve2 = repmat(mean(ve2), num_sim,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);
ve1_bit = double(ve1>= mean_ve1);
ve2_bit = double(ve2>= mean_ve2);

cmi_ve2 = zeros(1,length(snr));
cmi_ve1 = zeros(1,length(snr));
mi_bit = zeros(1,length(snr));
for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end
for i = 1:length(snr)
    cmi_ve2(i) = cmi(v1_bit(:,i),v2_bit(:,i),ve2_bit(:,i));
end

for i = 1:length(snr)
    cmi_ve1(i) = cmi(v1_bit(:,i),v2_bit(:,i),ve1_bit(:,i));
end

v1_trad_6 = v1;
v2_trad_6 = v2;
ve1_trad_6 = ve1;
ve2_trad_6 = ve2;
mi_trad_bit_6 = mi_bit;
cmi_trad_ve2_6 = cmi_ve2;
cmi_trad_ve1_6 = cmi_ve1;
%%
plot(snr,cmi_ve2_0,'r-o','LineWidth',1.5);   %%第一类窃听方法后的安全速率
hold on;
plot(snr,cmi_ve1_0,'-v','LineWidth',1.5,'Color', '#0072BD');
plot(snr,cmi_trad_ve2_0,'k-d','LineWidth',1.5);
plot(snr,cmi_ve2_6,'r--o','LineWidth',1.5);
plot(snr,cmi_ve1_6,'--v','LineWidth',1.5,'Color', '#0072BD');
plot(snr,cmi_trad_ve2_6,'k--d','LineWidth',1.5);
hold off;
grid on;
axis([0 30 0 0.95])
legend('cor = 0, V|Ve^1','cor = 0, V|Ve^2','cor = 0, V_h|Ve','cor = 0.6, V|Ve^1','cor = 0.6, V|Ve^2','cor = 0.6, V_h|Ve');
ylabel('Secret Key Rate (bit)');
xlabel('SNR (dB)');

plot(snr,mi_bit_0-cmi_ve2_0,'r-o','LineWidth',1.5);
hold on;
plot(snr,mi_bit_0-cmi_ve1_0,'-v','LineWidth',1.5,'Color', '#0072BD');
plot(snr,mi_bit_0- cmi_trad_ve2_0,'k-d','LineWidth',1.5);
plot(snr,mi_bit_6-cmi_ve2_6,'r--o','LineWidth',1.5);
plot(snr,mi_bit_6-cmi_ve1_6,'--v','LineWidth',1.5,'Color', '#0072BD');
plot(snr,mi_bit_6- cmi_trad_ve2_6,'k--d','LineWidth',1.5);
grid on;
axis([0 30 0 0.8])
legend('cor = 0, V|Ve^1','cor = 0, V|Ve^2','cor = 0, V_h|Ve','cor = 0.6, V|Ve^1','cor = 0.6, V|Ve^2','cor = 0.6, V_h|Ve');
ylabel('Leak Information (bit)');
xlabel('SNR (dB)');
hold off;
