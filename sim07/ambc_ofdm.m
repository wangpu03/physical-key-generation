% consists of all precedures of the OFDM 
% 
clear;
clc;

%% simulation parameters

n_fft = 128;        % IFFT/FFT size
n_cp = n_fft/4;    % size of cyclic prefix extension
n_frame = 3;        % the acount of the OFDM frame

% modulation methods: BPSK, QPSK,16QAM, 32QAM,64QAM
mod_method = 'QPSK';
% calculate modulation order from modulation method
mod_methods = {'BPSK', 'QPSK','8PSK','16QAM', '32QAM','64QAM'};
mod_order = find(ismember(mod_methods, mod_method));
% generate the data 
rand_ints_gen = randi(2,n_fft*mod_order*n_frame,1)-1;
save data_input.txt -ascii rand_ints_gen
rand_ints = load("data_input.txt");
rand_bits = dec2bin(rand_ints(:));

data_ofdm = ofdm_module(rand_bits, mod_method, n_fft, n_cp);

d = [8, 6, 3];
taps = [8, 6, 3];

h1 = ray_model(d(1),taps(1));
h2 = ray_model(d(2),taps(2));

h_12 = ray_model(d(3),taps(3));

cor_h1_he = 0;
cor_h12_he1 = 0;

he = eave_ray_model(h1,0,d(1));
h_e1 = eave_ray_model(h_12,0,d(3));



