function ofdm_data = ofdm_module(data,mod_method, n_fft, n_cp)
% modulate the data with OFDM
% inputs:
%       data: input raw data with bits, length = n_fft*mod_order*n_frame
%       mod_method: symbol modulation
%       n_fft: IFFT/FFT size
%       n_cp: size of cyclic prefix extension
%       n_frame: the acount of the OFDM frame
%output:
%       ofdm_data: modulated OFDM data

%n_ofdm = n_fft+n_cp;
% calculate modulation order from modulation method
mod_methods = {'BPSK', 'QPSK','8PSK','16QAM', '32QAM','64QAM'};
mod_order = find(ismember(mod_methods, mod_method));

% input data to binary stream
rand_bits = dec2bin(data(:));
im_bin = rand_bits(:);

% binary stream to symbols 
% Parse binary stream into mod_order bit symbols
% pads input signal to appropriate length
sym_rem = mod(mod_order-mod(length(im_bin),mod_order),mod_order);
padding = repmat('0',sym_rem,1);
im_bin_padded = [im_bin;padding];
cons_data = reshape(im_bin_padded, mod_order,length(im_bin_padded)/mod_order)';
cons_sym_id = bin2dec(cons_data);

% symbol modulation codebook
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
% pad input signal to appropriate length
fft_rem = mod(n_fft-mod(length(X),n_fft),n_fft);
X_padded = [X;zeros(fft_rem,1)];
X_blocks = reshape(X_padded, n_fft,length(X_padded)/n_fft);
x = ifft(X_blocks);

% add cyclic prefix extension and shift from parallel to serial
x_cp = [x(end-n_cp+1:end,:);x];
ofdm_data = x_cp(:);

