%% simulation of the signals received by BDs over OFDM carrier with raylaigh and exponential model
% in this file, we generate OFDM signal and adopt the exponential model to
% represent indoor channel.
% each channel tap is modeled by an independent complex Gaussian random
% variable with its average power that follows the exponential PDP, while
% taking the time index of each channel tap by the integer multiples of
% sampling periods.

% the direct OFDM signal
% the backscatter OFDM signal with backscatter modulation
%  *(AP)
%  |  \
%  |   \                %h=d^(-r/2)  r is path-loss exponet
%  4    5               %            d is the distance between two device
%  |     \              %5-0.0894, 6-0.068, 3-0.1925
%  *---3--*
%  BD1     BD2

clear;
clc;

%% the path loss of different channel
% simulate the path loss from the channel gain model
% we can adopt different channg gain model
f1_cg = 0.0894*5;
f2_cg = 0.0680*5;
h12_cg = 0.1925*5;

% backscatter coefficient
%alpha = sqrt(0.25);
alpha = 0.3+0.4i;

%% generate the OFDM signals
% simulation parameters
% modulation methods: BPSK, QPSK,16QAM, 32QAM,64QAM
mod_method = 'QPSK';
n_fft = 64;% IFFT/FFT size
n_cpe = 32; % size of cyclic prefix extension
n_ofdm = n_fft+n_cpe;% the size of a OFDM frame
n_frame = 2; % the acount of the OFDM frame
snr =10; % target SNR (dB)

% number of channel taps (1 = one channel) 
n_taps1 = 7;        %f1
n_taps2 = 6;        %f1
n_taps_h12 = 3;        %h12

% the maximum of the channel spread
n_L1 = max(n_taps1,n_taps2)+n_taps_h12-1;

% channel estimaition method: none, LS
ch_est_method = 'LS';

% option to save plot to file
save_file = 0;

% calculate modulation order from modulation method
mod_methods = {'BPSK', 'QPSK','8PSK','16QAM', '32QAM','64QAM'};
mod_order = find(ismember(mod_methods, mod_method));


% generate a vector whose elements are 0 and 1.
% rand_ints_gen = randi(2,n_fft*mod_order*n_frame,1)-1; 
% save data_input.txt -ascii rand_ints_gen

% input data to binary stream
rand_ints = load("data_input.txt");
rand_bits = dec2bin(rand_ints(:));
im_bin = rand_bits(:);

% binary stream to symbols 
% Parse binary stream into mod_order bit symbols
% pads input signal to appropriate length
sym_rem = mod(mod_order-mod(length(im_bin),mod_order),mod_order);
padding = repmat('0',sym_rem,1);
im_bin_padded = [im_bin;padding];
cons_data = reshape(im_bin_padded, mod_order,length(im_bin_padded)/mod_order)';
cons_sym_id = bin2dec(cons_data);

% symbol modulation, mainly generate teh symbol book
% BPSK
if mod_order == 1
    mod_ind = 2^(mod_order-1);
    n = 0:pi/mod_ind:2*pi-pi/mod_ind;
    in_phase = cos(n);
    quadrature = sin(n);
    symbol_book = (in_phase+quadrature*1i)';
end

% phase shift keying about unit circle
if mod_order == 2 || mod_order ==3
    mod_ind = 2^(mod_order-1);
    n = 0:pi/mod_ind:2*pi-pi/mod_ind;
    in_phase = cos(n+pi/4);
    quadrature = sin(n+pi/4);
    symbol_book = (in_phase+quadrature*1i)';
end

% 16QAM, 64QAM modulaiton
if mod_order ==4 || mod_order ==6
     mod_ind = 2^(mod_order-2);
     in_phase = repmat(linspace(-1,1,mod_ind),mod_ind,1);
     quadrature =  repmat(linspace(-1,1,mod_ind)',1,mod_ind);
     symbol_book = in_phase(:)+quadrature(:)*1i;
end

% 32 QAM modulation
% generate 6x6 constellation and removes corners
if mod_order ==5
    mod_ind = 6;
    in_phase = repmat(linspace(-1,1,mod_ind),mod_ind,1);
    quadrature =  repmat(linspace(-1,1,mod_ind)',1,mod_ind);
    symbol_book = in_phase(:)+quadrature(:)*1i;
    symbol_book = symbol_book([2:5 7:30 32:35]);
end

% modulate data according to symbol_book
X = symbol_book(cons_sym_id+1);

% use IFFT to move to time domain
%  and pad input signal to appropriate length
fft_rem = mod(n_fft-mod(length(X),n_fft),n_fft);
X_padded = [X;zeros(fft_rem,1)];
X_blocks = reshape(X_padded, n_fft,length(X_padded)/n_fft);
x = ifft(X_blocks);

% add cyclic prefix extension and shift from parallel to serial
x_cpe = [x(end-n_cpe+1:end,:);x];
x_s = x_cpe(:);
%*********************************************************************

%% add noise
%add awgn
% calculation data power
%x_s = repmat(1+1i,n_ofdm,1);

x_s = x_s/sqrt(mean(abs(x_s).^2));
mean(x_s)
data_pwr = mean(abs(x_s).^2);
data_pwr2 = mean(x_s.*conj(x_s));

% add noise to channel 
noise_pwr = data_pwr/10^(snr/10);
noise = normrnd(0,sqrt(noise_pwr/2),size(x_s))+normrnd(0,sqrt(noise_pwr/2),size(x_s))*1i;
%x_s_noise = x_s+noise;

% without noise
x_s_noise = x_s;

% measure SNR
snr_meas = 10*log10(mean(abs(x_s).^2)/mean(abs(noise).^2));
data_x_s_noise_pwr2 = mean(x_s_noise.*conj(x_s_noise));

%*********************************************************************

%% apply fading channel
% there is a big difference with or without the parameter 'same'
% with the 'same' returns only the central part of the convolution and its
% length is the same with the original 

% signal pass through the channel f1
pow_f1 = exp(-(0:n_taps1-1));
pow_f1 = f1_cg*pow_f1/norm(pow_f1);

h_f1=(randn(1,1)+1i*randn(1,1)).*sqrt(pow_f1/2);
data_h_f1_pwr = mean(h_f1.*conj(h_f1));

x_s_noise_fading1 = conv(x_s_noise,h_f1);
data_x_s_noise_fading1_pwr = mean(abs(x_s_noise_fading1).^2);
data_x_s_noise_fading1_pwr02 = mean(x_s_noise_fading1.*conj(x_s_noise_fading1));

% signal pass through the channel f2
pow_f2 = exp(-(0:n_taps2-1));
pow_f2 = f2_cg*pow_f2/norm(pow_f2);

h_f2=(randn(1,1)+1i*randn(1,1)).*sqrt(pow_f2/2);
data_h_f2_pwr = mean(h_f2.*conj(h_f2));

x_s_noise_fading2 = conv(x_s_noise,h_f2);
data_x_s_noise_fading2_pwr = mean(abs(x_s_noise_fading2).^2);

% the designing of the backscatter waveform at the user 1
bc_signal_fh = ones(n_ofdm/2,1);
bc_signal_bh = repmat(-1,n_ofdm/2,1);
bc_signal = [bc_signal_fh' bc_signal_bh']';
bc_signal = repmat(bc_signal,n_frame,1);

% backscatter operation of the signal from user1
sym_rem = mod(length(x_s_noise_fading1),n_ofdm);
padding_bc = ones(sym_rem,1);
bc_signal = [bc_signal',padding_bc']';
x_s_noise_fading1_bc = alpha*bc_signal.*x_s_noise_fading1;
data_x_s_noise_fading1_bc_pwr = mean(abs(x_s_noise_fading1_bc).^2);

% the designing of the backscatter waveform at the user 2
bc_signal_fh = ones(n_ofdm/2,1);
bc_signal_bh = repmat(-1,n_ofdm/2,1);
bc_signal = [bc_signal_fh' bc_signal_bh']';
bc_signal = repmat(bc_signal,n_frame,1);

% backscatter operation of the signal from user2
sym_rem = mod(length(x_s_noise_fading2),n_ofdm);
padding_bc = ones(sym_rem,1);
bc_signal = [bc_signal',padding_bc']';
x_s_noise_fading2_bc = alpha*bc_signal.*x_s_noise_fading2;
data_x_s_noise_fading2_bc_pwr = mean(abs(x_s_noise_fading2_bc).^2);

% the channel h12 model
pow_h12 = exp(-(0:n_taps_h12-1));
pow_h12 = h12_cg*pow_h12/norm(pow_h12);

h_h12=(randn(1,1)+1i*randn(1,1)).*sqrt(pow_h12/2);
h_h21 = h_h12;
data_h_h12_pwr = mean(abs(h_h12).^2);

% the signal from user1 passing through channel h21
x_s_noise_fading1_bc_h21 = conv(x_s_noise_fading1_bc,h_h21);
data_x_s_noise_fading1_bc_h21_pwr = mean(abs(x_s_noise_fading1_bc_h21).^2);

% the signal from user2 passing through channel h12
x_s_noise_fading2_bc_h12 = conv(x_s_noise_fading2_bc,h_h12);
data_x_s_noise_fading2_bc_h12_pwr = mean(abs(x_s_noise_fading2_bc_h12).^2);

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

%% calculate the desired information at receivering devices
% user 1
zb1_s_noise = x_s_noise_bd1(n_L1:n_cpe)-x_s_noise_bd1(n_fft+n_L1:n_fft+n_cpe);
for n = 1:1:n_frame-1
    zb1_s_noise =[zb1_s_noise(:)' (x_s_noise_bd1(n*n_ofdm+n_L1:n*n_ofdm+n_cpe)-x_s_noise_bd1(n*n_ofdm+n_fft+n_L1:n*n_ofdm+n_fft+n_cpe))']';
end
%zb1_s_noise = zb1_s_noise/2;

zd1_s_noise = x_s_noise_bd1(n_L1:n_cpe)+x_s_noise_bd1(n_fft+n_L1:n_fft+n_cpe);
for n = 1:1:n_frame-1
    zd1_s_noise =[zd1_s_noise(:)' (x_s_noise_bd1(n*n_ofdm+n_L1:n*n_ofdm+n_cpe)+x_s_noise_bd1(n*n_ofdm+n_fft+n_L1:n*n_ofdm+n_fft+n_cpe))']';
end
%zd1_s_noise = zd1_s_noise/2;

% user 2
zb2_s_noise = x_s_noise_bd2(n_L1:n_cpe)-x_s_noise_bd2(n_fft+n_L1:n_fft+n_cpe);
for n = 1:1:n_frame-1
    zb2_s_noise =[zb2_s_noise(:)' (x_s_noise_bd2(n*n_ofdm+n_L1:n*n_ofdm+n_cpe)-x_s_noise_bd2(n*n_ofdm+n_fft+n_L1:n*n_ofdm+n_fft+n_cpe))']';
end
%zb2_s_noise = zb2_s_noise/2;

zd2_s_noise = x_s_noise_bd2(n_L1:n_cpe)+x_s_noise_bd2(n_fft+n_L1:n_fft+n_cpe);
for n = 1:1:n_frame-1
    zd2_s_noise =[zd2_s_noise(:)' (x_s_noise_bd2(n*n_ofdm+n_L1:n*n_ofdm+n_cpe)+x_s_noise_bd2(n*n_ofdm+n_fft+n_L1:n*n_ofdm+n_fft+n_cpe))']';
end
%zd2_s_noise = zd2_s_noise/2;

%% virtual link information between two BDs
%user1
v12_s_noise = zb1_s_noise.*zd1_s_noise;

v12_s_noise_ch = zb1_s_noise.*zd1_s_noise.*conj(x_s(n_L1:n_cpe)).*conj(x_s(n_L1:n_cpe))./(x_s(n_L1:n_cpe).*conj(x_s(n_L1:n_cpe))).^2;
v12_s_noise_ch_pwr = mean(abs(v12_s_noise_ch).^2);

zb1_s_noise_pwr = mean(zb1_s_noise.*conj(zb1_s_noise));
zd1_s_noise_pwr = mean(zd1_s_noise.*conj(zd1_s_noise));
v12_s_noise_pwr = mean(v12_s_noise.*conj(v12_s_noise));

%user2
v21_s_noise = zb2_s_noise.*zd2_s_noise;
zb2_s_noise_pwr = mean(abs(zb2_s_noise).^2);
zd2_s_noise_pwr = mean(abs(zd2_s_noise).^2);
v21_s_noise_pwr = mean(abs(v21_s_noise).^2);

v21_s_noise_ch = zb2_s_noise.*zd2_s_noise.*conj(x_s(n_L1:n_cpe)).*conj(x_s(n_L1:n_cpe))./(x_s(n_L1:n_cpe).*conj(x_s(n_L1:n_cpe))).^2;
v21_s_noise_ch_pwr = mean(abs(v21_s_noise_ch).^2);