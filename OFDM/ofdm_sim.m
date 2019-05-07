%% simulation of OFDM communication 
% consists of all precedures of the OFDM 
% 
clear;
clc;

%% simulation parameters
% modulation methods: BPSK, QPSK,16QAM, 32QAM,64QAM
mod_method = 'QPSK';

% IFFT/FFT size
n_fft = 128;

% size of cyclic prefix extension
n_cpe = 32;

% target SNR (dB)
snr =20;

% number of channel taps (1 = no channel)
n_taps = 8;

% channel estimaition method: none, LS
ch_est_method = 'LS';

% option to save plot to file
save_file = 0;

% calculate modulation order from modulation method
mod_methods = {'BPSK', 'QPSK','8PSK','16QAM', '32QAM','64QAM'};
mod_order = find(ismember(mod_methods, mod_method));

%% input data to binary stream
im  = imread('Picture1.bmp');
im_bin = dec2bin(im(:))';
im_bin = im_bin(:);

% zero_ints = zeros(128,1)+255;
% zero_bits = dec2bin(zero_ints(:));
% im_bin = zero_bits(:);

% rand_ints = randi(2,25600,1)-1;
% rand_bits = dec2bin(rand_ints(:));
% im_bin = rand_bits(:);
%% binary stream to symbols 
% Parse binary stream into mod_order bit symbols
% pads input signal to appropriate length
sym_rem = mod(mod_order-mod(length(im_bin),mod_order),mod_order);
padding = repmat('0',sym_rem,1);
im_bin_padded = [im_bin;padding];
cons_data = reshape(im_bin_padded, mod_order,length(im_bin_padded)/mod_order)';
cons_sym_id = bin2dec(cons_data);

%% symbol modulation
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

%% use IFFT to move to time domain
% pad input signal to appropriate length
fft_rem = mod(n_fft-mod(length(X),n_fft),n_fft);
X_padded = [X;zeros(fft_rem,1)];
X_blocks = reshape(X_padded, n_fft,length(X_padded)/n_fft);
x = ifft(X_blocks);

% add cyclic prefix extension and shift from parallel to serial
x_cpe = [x(end-n_cpe+1:end,:);x];
x_s = x_cpe(:);

%% add awgn
% calculation data power
data_pwr = mean(abs(x_s).^2);

% calculate the cp data power for n frame OFDM
x_s_cp = x_s(9:32);
for n = 1:1:9
    x_s_cp =[x_s_cp(:)' x_s(160*n+9:160*n+32)']';
end
data_cpe_pwr = mean(abs(x_s_cp).^2);

% data_pwr_cpe_add=0;
% for i = 1:1:5
%     data_pwr_cpe_add = data_pwr_cpe_add + mean(abs(x_s(160*(i-1)+9:160*(i-1)+32).^2));
% end
% data_pwr_cpe = data_pwr_cpe_add/i;
%data_pwr_cpe = mean(abs(x_s(241:320).^2));

% add noise to channel 
noise_pwr = data_pwr/10^(snr/10);
noise = normrnd(0,sqrt(noise_pwr/2),size(x_s))+normrnd(0,sqrt(noise_pwr/2),size(x_s))*1i;
x_s_noise = x_s+noise;

%x_s_noise = x_s;

% measure SNR
snr_meas = 10*log10(mean(abs(x_s.^2))/mean(abs(noise.^2)));

%% apply fading channel
% there is a big difference with or without the parameter 'same'
% with the 'same' returns only the central part of the convolution
g = exp(-(0:n_taps-1));
g = g/norm(g);
x_s_noise_fading = conv(x_s_noise,g,'same');

%x_s_noise_fading = conv(x_s_noise,g);

% calculation data power after channel 
% data_pwr_cpe_add_ac=0;
% for i = 1:1:5
%     data_pwr_cpe_add_ac = data_pwr_cpe_add_ac + mean(abs(x_s_noise_fading(160*(i-1)+9:160*(i-1)+32).^2));
% end
% data_pwr_cpe_ac = data_pwr_cpe_add_ac/i;
%data_pwr_cpe = mean(abs(x_s(241:320).^2));

data_pwr_ac = mean(abs(x_s_noise_fading.^2));

%% use FFT to move to frequency domain
% remmove cyclic prefix extension and shift from serial to parallel
x_p = reshape(x_s_noise_fading,n_fft+n_cpe,length(x_s_noise_fading)/(n_fft+n_cpe));
x_p_cpr = x_p(n_cpe+1:end,:);

% move to frequency domain
X_hat_blocks = fft(x_p_cpr);

%% estimate channel
if n_taps > 1
    switch(ch_est_method)
        case 'none'
        case 'LS'
            G = X_hat_blocks(:,1)./X_blocks(:,1);
            X_hat_blocks = X_hat_blocks./repmat(G,1,size(X_hat_blocks,2));
    end
end

%% symbol demodulation
% remove fft padding
X_hat = X_hat_blocks(:);
X_hat = X_hat(1:end-fft_rem);

% recover data from modulated symbols
rec_syms = knnsearch ([real(symbol_book) imag(symbol_book)],[real(X_hat),imag(X_hat)])-1;

% parse to binary stream & remove symbol padding 
rec_syms_cons = dec2bin(rec_syms);
rec_im_bin = reshape(rec_syms_cons',numel(rec_syms_cons),1);
rec_im_bin = rec_im_bin(1:end-sym_rem);
ber = sum(abs(rec_im_bin - im_bin))/length(im_bin);

%% recover image
rec_im = reshape(rec_im_bin,8,numel(rec_im_bin)/8);
rec_im = uint8(bin2dec(rec_im'));
rec_im = reshape(rec_im,size(im));

%% generate plots
% transmit constellation
subplot(2,2,1)
plot(X,'x','linewidth',2,'markersize',10);
xlim([-2 2]);
ylim([-2 2]);
xlabel('In Phase');
ylabel('Quadrature');
if n_taps > 1
    title(sprintf('\\bfTransmit Constellation\n\\rm%s Modulation\nMultipath Channel taps: %d', mod_method, n_taps));
else
    title(sprintf('\\bfTransmit Constellation\n\\rm%s Modulation', mod_method));
end 
grid on;

% recovered constellation
subplot(2,2,2)
plot(X_hat(1:500:end),'x','markersize',3);
xlim([-2 2]);
ylim([-2 2]);
xlabel('In Phase');
ylabel('Quadrature');
if n_taps > 1
    title(sprintf('\\bfReceived Constellation\n\\rmMeasured SNR: %0.2f dB\nChannel estimation: %s',snr_meas,ch_est_method));
else
    title(sprintf('\\bfReceived Constellation\n\\rmMeasured SNR: %0.2f dB',snr_meas));
end 
grid on;

% original image
subplot(2,2,3);
imshow(im);
title('\bfTransmit Image');

% recovered image
subplot(2,2,4);
imshow(rec_im);
title(sprintf('\\bfRecovered Image\n\\rmBER: %02g',ber));

%position figure
%set(gcf,'Position',[680,287,698,691]);

%save figure
if save_file
    print(sprintf('Plots/%s_%.offft_%.0fcpe_%0fdB_%.0ftaps_%s',mod_method,n_fft,n_cpe,snr,n_taps,ch_est_method),'_dmata');
end

