% consists of all precedures of the OFDM 
% 
clear;
clc;

%% simulation parameters

n_fft = 128;        % IFFT/FFT size
n_cp = n_fft/4;    % size of cyclic prefix extension
n_frame = 3;        % the acount of the OFDM frame
n_ofdm = n_fft + n_cp;

% modulation methods: BPSK, QPSK,16QAM, 32QAM,64QAM
mod_method = 'QPSK';
% calculate modulation order from modulation method
mod_methods = {'BPSK', 'QPSK','8PSK','16QAM', '32QAM','64QAM'};
mod_order = find(ismember(mod_methods, mod_method));
% generate the data 
rand_ints_gen = randi(2,n_fft*mod_order*n_frame,1)-1;
save data_input.txt -ascii rand_ints_gen
rand_ints = load("data_input.txt");

data_ofdm = ofdm_module(rand_ints, mod_method, n_fft, n_cp);
size(data_ofdm)

d = [8, 7, 3];
taps = [8, 7, 3];

h1 = ray_model(d(1),taps(1));
h2 = ray_model(d(2),taps(2));

h_12 = ray_model(d(3),taps(3));

cor_h1_he = 0;
cor_h12_he1 = 0;

he = eave_ray_model(h1,0,d(1));
h_e1 = eave_ray_model(h_12,0,d(3));

pt = 1;
snr = 15;
alpha = 0.3 + 1i*4;
pw_noise = pt/10^(snr/10);
y1_d = ofdm_trans(data_ofdm, h1, pw_noise);
back1 = ofdm_back(y1_d, alpha, n_ofdm, n_frame);

y2_d = ofdm_trans(data_ofdm, h2, pw_noise);
back2 = ofdm_back(y2_d, alpha, n_ofdm, n_frame);

y1_b = ofdm_trans(back2, h_12, pw_noise);
y2_b = ofdm_trans(back1, h_12, pw_noise);

%% combinaiton of two signals
% user 1
% x_s_noise_fading1 and x_s_noise_fading2_h12
sym_rem = mod(n_ofdm-mod(length(x_s_noise_fading1),n_ofdm),n_ofdm);
padding = repmat(0+0i,sym_rem,1);
x_s_noise_fading1_padded = [x_s_noise_fading1;padding];

sym_rem = mod(n_ofdm-mod(length(x_s_noise_fading2_bc_h12),n_ofdm),n_ofdm);
padding = repmat(0+0i,sym_rem,1);
x_s_noise_fading2_bc_h12_padded = [x_s_noise_fading2_bc_h12;padding];

x_s_noise_bd1 = x_s_noise_fading1_padded+x_s_noise_fading2_bc_h12_padded;
x_s_noise_bd1_pwr = mean(abs(x_s_noise_bd1).^2);

% user 2
%x_s_noise_fading2 and x_s_noise_fading1_h21
sym_rem = mod(n_ofdm-mod(length(x_s_noise_fading2),n_ofdm),n_ofdm);
padding = repmat(0+0i,sym_rem,1);
x_s_noise_fading2_padded = [x_s_noise_fading2;padding];

sym_rem = mod(n_ofdm-mod(length(x_s_noise_fading1_bc_h21),n_ofdm),n_ofdm);
padding = repmat(0+0i,sym_rem,1);
x_s_noise_fading1_bc_h21_padded = [x_s_noise_fading1_bc_h21;padding];

x_s_noise_bd2 = x_s_noise_fading2_padded+x_s_noise_fading1_bc_h21_padded;

x_s_noise_bd2_pwr = mean(abs(x_s_noise_bd2).^2);




