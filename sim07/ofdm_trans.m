function channel_data = ofdm_trans(ofdm_data, h, noise_pw)
% inputs:
%       ofdm_data: modulated data for transmitting
%       h: the Rayleigh channel that data will go through
%       noise_pw: the noise power recieved at receiver
% output: 
%       channel data: the received data at receiver, including data signal and
%       noise signal

channel_data = conv(ofdm_data,h);


noise = normrnd(0,sqrt(noise_pw),size(channel_data))+normrnd(0,sqrt(noise_pw),size(channel_data))*1i;
channel_data = channel_data+noise;