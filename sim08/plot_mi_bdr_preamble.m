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

num_mc = 50000;

snr = 0:3:30;
v1 = zeros(num_mc, length(snr));
v2 = zeros(num_mc, length(snr));

%%
% nfft = 256, K =1
n_fft = 256;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=1;
c_flag = 1;
% generate the data 
rand_ints_1256 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_1256, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr(data_ofdm, n_ofdm,n_cp, K, snr(i));
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_cp = v1;
v2_cp = v2;
v1_bit_cp = v1_bit;
v2_bit_cp = v2_bit;
mi_bit_cp = mi_bit;


%%
% nfft = 256, K =1， ZP=2
n_fft = 256;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=1;
c_flag = 2;
% generate the data 
rand_ints_1256 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_1256, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr(data_ofdm, n_ofdm,n_cp, K, snr(i));
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_zp = v1;
v2_zp = v2;

v1_bit_zp = v1_bit;
v2_bit_zp = v2_bit;
mi_bit_zp = mi_bit;

%%
% nfft = 256, K =1， OP=3
n_fft = 256;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=1;
c_flag = 3;
% generate the data 
rand_ints_1256 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_1256, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr(data_ofdm, n_ofdm,n_cp, K, snr(i));
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_op = v1;
v2_op = v2;

v1_bit_op = v1_bit;
v2_bit_op = v2_bit;
mi_bit_op = mi_bit;

%%
% nfft = 256, K =1， preamble=4
n_fft = 256;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=1;
c_flag = 4;
% generate the data 
rand_ints_1256 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_1256, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr(data_ofdm, n_ofdm,n_cp, K, snr(i));
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_pre = v1;
v2_pre = v2;

v1_bit_pre = v1_bit;
v2_bit_pre = v2_bit;
mi_bit_pre = mi_bit;


%%
snr = 0:3:30;
plot(snr, mi_bit_cp,'-d','LineWidth',1.5,'Color', '#77AC30');
hold on;
plot(snr, mi_bit_zp,'r-x','LineWidth',1.5);
hold on;
plot(snr, mi_bit_op,'k-o','LineWidth',1.5);
hold on;
plot(snr, mi_bit_pre,'-s','LineWidth',1.5,'Color', '#0072BD');
grid on;
axis([0 30 0.3 0.95]);
legend('CP','ZP','OP','1-0 bits');
%legend('K = 1, N_{FFT} = 256','K = 3, N_{FFT} = 256','K = 1, N_{FFT} = 256','K = 3, N_{FFT} = 256','K = 1, N_{Guard} = 256','K = 1, N_{Guard} = 256','K = 1, N_{Guard}= 256,Noiseless');
ylabel('每比特互信息（bits）','Fontname','<宋体>');
xlabel('SNR (dB)');

%% 计算BDR
bdr_cp = 1 - sum(v1_bit_cp == v2_bit_cp)/num_mc;
bdr_zp = 1 - sum(v1_bit_zp == v2_bit_zp)/num_mc;
bdr_op = 1 - sum(v1_bit_op == v2_bit_op)/num_mc;
bdr_pre = 1 - sum(v1_bit_pre == v2_bit_pre)/num_mc;

%%
snr = 0:3:30;
plot(snr, bdr_cp,'-d','LineWidth',1.5,'Color', '#77AC30');
hold on;
plot(snr, bdr_zp,'r-x','LineWidth',1.5);
hold on;
plot(snr, bdr_op,'k-o','LineWidth',1.5);
hold on;
plot(snr, bdr_pre,'-s','LineWidth',1.5,'Color', '#0072BD');
grid on;
axis([0 30 0.0 0.17]);
legend('CP','ZP','OP','1-0 bits');
%legend('K = 1, N_{FFT} = 256','K = 3, N_{FFT} = 256','K = 1, N_{FFT} = 256','K = 3, N_{FFT} = 256','K = 1, N_{Guard} = 256','K = 1, N_{Guard} = 256','K = 1, N_{Guard}= 256,Noiseless');
ylabel('互异比特比率','Fontname','<宋体>');
xlabel('SNR (dB)');
datestr(now)