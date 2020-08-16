function ofdm_data = ofdm_module(data,mod_method, n_fft, n_cp, cp_flag)
% modulate the data with OFDM
% inputs:
%       data: input raw data with bits, length = n_fft*mod_order*n_frame
%       mod_method: symbol modulation
%       n_fft: IFFT/FFT size
%       n_cp: size of cyclic prefix extension
%       n_frame: the acount of the OFDM frame
%       cp_flag: 1--cp, 0--guard
%output:
%       ofdm_data: modulated OFDM data

if nargin < 4
    cp_flag = 1;
end

%n_ofdm = n_fft+n_cp;
% calculate modulation order from modulation method
mod_methods = {'BPSK', 'QPSK','8PSK','16QAM', '32QAM','64QAM'};
mod_order = find(ismember(mod_methods, mod_method));

% input data to binary stream
rand_bits = dec2bin(data(:));

% binary stream to symbols 
% Parse binary stream into mod_order bit symbols
% pads input signal to appropriate length
sym_rem = mod(mod_order-mod(length(rand_bits),mod_order),mod_order);
padding = repmat('0',sym_rem,1);
rand_bits_padded = [rand_bits;padding];
cons_data = reshape(rand_bits_padded, mod_order,length(rand_bits_padded)/mod_order)';
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
%size_X = size(X)

% use IFFT to move to time domain
% pad input signal to appropriate lengthmod(length(X),n_fft)
fft_rem = mod(n_fft-mod(length(X),n_fft),n_fft);
X_padded = [X;zeros(fft_rem,1)];
X_blocks = reshape(X_padded, n_fft,length(X_padded)/n_fft);
x = ifft(X_blocks);
%data_x_pwr = mean(abs(x).^2)
% add cyclic prefix extension and shift from parallel to serial
guard_zp = zeros(n_cp,1);
guard_op = ones(n_cp,1)*(1+1i)/sqrt(2);
guard_pre = [0;(1+1i)/sqrt(2)];
guard_pre = repmat(guard_pre,n_cp/2,1);

if cp_flag == 1
    x_cp = [x(end-n_cp+1:end,:);x];
elseif cp_flag == 2
    x_cp = [guard_zp;x(1:n_fft-n_cp,:);guard_zp];
elseif cp_flag == 3
    x_cp = [guard_op;x(1:n_fft-n_cp,:);guard_op];
elseif cp_flag == 4
    x_cp = [guard_pre;x(1:n_fft-n_cp,:);guard_pre];
end

ofdm_data = x_cp(:);

ofdm_data = ofdm_data/sqrt(mean(abs(ofdm_data).^2));


