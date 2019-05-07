clc;
clear;
pn = 0.01;
snr = 10:2:30;
pt = 10.^(snr/10)*pn;

num = length(pt);
num_mc = 100000;
v12 = zeros(num, num_mc);
v21 = zeros(num, num_mc);
x_bd1 = zeros(num, num_mc);
x_bd2 = zeros(num, num_mc);

for j = 1:num
    j
    for i = 1: num_mc
        [v12(j,i),v21(j,i),x_bd1(j,i),x_bd2(j,i)] = func_ambc_ofdm_ray_exe(pt(j));
    end
end
corr2(v12(1,:),v21(1,:))

figure(1);
plot([1:num_mc],v12(1,:));
hold on;
plot([1:num_mc],v21(1,:));

mean(v12(10,:))

mean(x_bd1(4,:))
mean(x_bd2(1,:))