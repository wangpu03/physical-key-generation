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
% nfft = 256, K =1， cor = 1;
n_fft = 256;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=1;
c_flag = 1;
cor = 1;
% generate the data 
rand_ints_1256 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_1256, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr_ne(data_ofdm, n_ofdm,n_cp, K, snr(i), cor);
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_cp_cor1 = v1;
v2_cp_cor1 = v2;
v1_cp_bit_cor1 = v1_bit;
v2_cp_bit_cor1 = v2_bit;
mi_cp_bit_cor1 = mi_bit;

%%
% nfft = 256, K =1， cor = 0.999;
n_fft = 256;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=1;
c_flag = 1;
cor = 0.999;
% generate the data 
rand_ints_1256 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_1256, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr_ne(data_ofdm, n_ofdm,n_cp, K, snr(i), cor);
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_cp_cor9 = v1;
v2_cp_cor9 = v2;
v1_cp_bit_cor9 = v1_bit;
v2_cp_bit_cor9= v2_bit;
mi_cp_bit_cor9 = mi_bit;
%%
% nfft = 256, K =1， cor = 0.998;
n_fft = 256;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=1;
c_flag = 1;
cor = 0.996;
% generate the data 
rand_ints_1256 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_1256, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr_ne(data_ofdm, n_ofdm,n_cp, K, snr(i), cor);
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_cp_cor8 = v1;
v2_cp_cor8 = v2;
v1_cp_bit_cor8 = v1_bit;
v2_cp_bit_cor8= v2_bit;
mi_cp_bit_cor8 = mi_bit;

%%
% nfft = 256, K =1， cor = 1; preamble
n_fft = 256;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=1;
c_flag = 4;
cor = 1;
% generate the data 
rand_ints_1256 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_1256, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr_ne_trad(data_ofdm, n_ofdm,n_cp, K, snr(i), cor);
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_pre_cor1 = v1;
v2_pre_cor1 = v2;
v1_pre_bit_cor1 = v1_bit;
v2_pre_bit_cor1 = v2_bit;
mi_pre_bit_cor1 = mi_bit;

%%
% nfft = 256, K =1， cor = 0.999; preamble
n_fft = 256;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=1;
c_flag = 4;
cor = 0.999;
% generate the data 
rand_ints_1256 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_1256, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr_ne_trad(data_ofdm, n_ofdm,n_cp, K, snr(i), cor);
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_pre_cor9 = v1;
v2_pre_cor9 = v2;
v1_pre_bit_cor9 = v1_bit;
v2_pre_bit_cor9= v2_bit;
mi_pre_bit_cor9 = mi_bit;
%%
% nfft = 256, K =1， cor = 0.998; preamble
n_fft = 256;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=1;
c_flag = 4;
cor = 0.996;
% generate the data 
rand_ints_1256 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_1256, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr_ne_trad(data_ofdm, n_ofdm,n_cp, K, snr(i), cor);
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_pre_cor8 = v1;
v2_pre_cor8 = v2;
v1_pre_bit_cor8 = v1_bit;
v2_pre_bit_cor8= v2_bit;
mi_pre_bit_cor8 = mi_bit;


% %%
% % nfft = 256, K =1， cor = 1; preamble
% n_fft = 256;
% n_cp = n_fft/4;    % size of cyclic prefix extension
% n_ofdm = n_fft + n_cp;
% K=1;
% c_flag = 4;
% cor = 1;
% % generate the data 
% rand_ints_1256 = rand_ints(1:n_fft*mod_order*K,:);
% data_ofdm = ofdm_module(rand_ints_1256, mod_method, n_fft, n_cp, c_flag);
% 
% for i = 1:length(snr)
% for j = 1:num_mc
%     [v1(j,i),v2(j,i)] = mi_snr_ne(data_ofdm, n_ofdm,n_cp, K, snr(i), cor);
% end
% end
% 
% mean_v1 = repmat(mean(v1), num_mc,1);
% mean_v2 = repmat(mean(v2), num_mc,1);
% 
% v1_bit = double(v1>= mean_v1);
% v2_bit = double(v2>= mean_v2);
% 
% for i = 1:length(snr)
%     mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
% end
% 
% v1_pre_cor1 = v1;
% v2_pre_cor1 = v2;
% v1_pre_bit_cor1 = v1_bit;
% v2_pre_bit_cor1 = v2_bit;
% mi_pre_bit_cor1 = mi_bit;
% 
% %%
% % nfft = 256, K =1， cor = 0.999; preamble
% n_fft = 256;
% n_cp = n_fft/4;    % size of cyclic prefix extension
% n_ofdm = n_fft + n_cp;
% K=1;
% c_flag = 4;
% cor = 0.999;
% % generate the data 
% rand_ints_1256 = rand_ints(1:n_fft*mod_order*K,:);
% data_ofdm = ofdm_module(rand_ints_1256, mod_method, n_fft, n_cp, c_flag);
% 
% for i = 1:length(snr)
% for j = 1:num_mc
%     [v1(j,i),v2(j,i)] = mi_snr_ne(data_ofdm, n_ofdm,n_cp, K, snr(i), cor);
% end
% end
% 
% mean_v1 = repmat(mean(v1), num_mc,1);
% mean_v2 = repmat(mean(v2), num_mc,1);
% 
% v1_bit = double(v1>= mean_v1);
% v2_bit = double(v2>= mean_v2);
% 
% for i = 1:length(snr)
%     mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
% end
% 
% v1_pre_cor9 = v1;
% v2_pre_cor9 = v2;
% v1_pre_bit_cor9 = v1_bit;
% v2_pre_bit_cor9= v2_bit;
% mi_pre_bit_cor9 = mi_bit;
% %%
% % nfft = 256, K =1， cor = 0.998; preamble
% n_fft = 256;
% n_cp = n_fft/4;    % size of cyclic prefix extension
% n_ofdm = n_fft + n_cp;
% K=1;
% c_flag = 4;
% cor = 0.998;
% % generate the data 
% rand_ints_1256 = rand_ints(1:n_fft*mod_order*K,:);
% data_ofdm = ofdm_module(rand_ints_1256, mod_method, n_fft, n_cp, c_flag);
% 
% for i = 1:length(snr)
% for j = 1:num_mc
%     [v1(j,i),v2(j,i)] = mi_snr_ne(data_ofdm, n_ofdm,n_cp, K, snr(i), cor);
% end
% end
% 
% mean_v1 = repmat(mean(v1), num_mc,1);
% mean_v2 = repmat(mean(v2), num_mc,1);
% 
% v1_bit = double(v1>= mean_v1);
% v2_bit = double(v2>= mean_v2);
% 
% for i = 1:length(snr)
%     mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
% end
% 
% v1_pre_cor8 = v1;
% v2_pre_cor8 = v2;
% v1_pre_bit_cor8 = v1_bit;
% v2_pre_bit_cor8= v2_bit;
% mi_pre_bit_cor8 = mi_bit;
% 
snr = 0:3:30;
plot(snr, mi_cp_bit_cor1,'-d','LineWidth',1.5,'Color', '#77AC30');
hold on;
plot(snr, mi_cp_bit_cor9,'r-x','LineWidth',1.5);
plot(snr, mi_cp_bit_cor8,'k-o','LineWidth',1.5);
plot(snr, mi_pre_bit_cor1,'--d','LineWidth',1.5,'Color', '#77AC30');
plot(snr, mi_pre_bit_cor9,'r--x','LineWidth',1.5);
plot(snr, mi_pre_bit_cor8,'k--o','LineWidth',1.5);
hold off;
grid on;
axis([0 30 0.5 0.95]);
legend('cor = 1,    triangle','cor = 0.9, triangle','cor = 0.8, triangle','cor = 1,    inward','cor = 0.9, inward','cor = 0.8, inward');
ylabel('Mutual Information per bit');
xlabel('SNR (dB)');



%% 计算BDR
bdr_cp_cor1 = 1 - sum(v1_cp_bit_cor1 == v2_cp_bit_cor1)/num_mc;
bdr_cp_cor9 = 1 - sum(v1_cp_bit_cor9 == v2_cp_bit_cor9)/num_mc;
bdr_cp_cor8 = 1 - sum(v1_cp_bit_cor8 == v2_cp_bit_cor8)/num_mc;
bdr_pre_cor1 = 1 - sum(v1_pre_bit_cor1 == v2_pre_bit_cor1)/num_mc;
bdr_pre_cor9 = 1 - sum(v1_pre_bit_cor9 == v2_pre_bit_cor9)/num_mc;
bdr_pre_cor8 = 1 - sum(v1_pre_bit_cor8 == v2_pre_bit_cor8)/num_mc;


%%
snr = 0:3:30;
plot(snr, bdr_cp_cor1,'-d','LineWidth',1.5,'Color', '#77AC30');
hold on;
plot(snr, bdr_cp_cor9,'r-x','LineWidth',1.5);
plot(snr, bdr_cp_cor8,'k-o','LineWidth',1.5);
plot(snr, bdr_pre_cor1,'--d','LineWidth',1.5,'Color', '#77AC30');
plot(snr, bdr_pre_cor9,'r--x','LineWidth',1.5);
plot(snr, bdr_pre_cor8,'k--o','LineWidth',1.5);
hold off;
grid on;
axis([0 30 0.0 0.06]);
legend('cor = 1,    triangle','cor = 0.9, triangle','cor = 0.8, triangle','cor = 1,    inward','cor = 0.9, inward','cor = 0.8, inward');
ylabel('Bit Disagreement Ratio');
xlabel('SNR (dB)');