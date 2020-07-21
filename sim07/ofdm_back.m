function back_data = ofdm_back(data, alpha, n_ofdm, n_frame)
% inputs:
%       data: the ofdm data before backscatter 
%       alpha: backscatter coefficient
%       n_ofdm: the number of IFFT of a ofdm frame
%       n_frame: the number of frame of the all data
% outputs:
%       back_data: the ofdm data after backscatter


% the designing of the backscatter waveform at the user 1
bc_signal_fh = ones(n_ofdm/2,1);
bc_signal_bh = ones(n_ofdm/2,1)-2;
bc_signal = [bc_signal_fh; bc_signal_bh];
bc_signal = repmat(bc_signal,n_frame,1);

% backscatter operation of the signal from user1
sym_rem = mod(length(data),n_ofdm);
padding_bc = zeros(sym_rem,1);   %pad zero
bc_signal = [bc_signal;padding_bc];
%size(bc_signal)
%size(data)
back_data = alpha*bc_signal.*data;
