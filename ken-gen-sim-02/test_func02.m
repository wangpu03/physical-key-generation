%  set the parameters
%  channels, backscatter coefficient, energy harvesting coefficient and
%  the transmit power
%  the form of the channel model f=g*d^(L/2)--(d-distance,g-CSCG, f-channel)

alpha = 0.4;        %the constant of backscatter coefficient
Pt = 10000;           %transmission power of the RF source
rho = 0.2;          %power split ratio for information receiver, range is [0,1]

theta = 1;
d_f1 = 5.1;
d_f2 = 6;
d_h12 = 3;
sigma_f1 = sqrt(theta*(d_f1)^(-3/2));    %channel f1
f1 = normrnd(0,sigma_f1,[1 10000000])';

sigma_f2 = sqrt(theta*(d_f2)^(-3/2));    %channel f2
f2 = normrnd(0,sigma_f2,[1 10000000])';

sigma_h12 = sqrt(theta*(d_h12)^(-3/2));          %channel h12 or h21
h12 = normrnd(0,sigma_h12,[1 10000000])';

sigma_z1 = sqrt(1);       % noise at the device BD1
z1 = normrnd(0,sigma_z1,[1 10000000])';

sigma_z2 = sqrt(1);       % noise at the device BD2
z2 = normrnd(0,sigma_z2,[1 10000000])';

b12 = alpha*h12.*f2;
b21 = alpha*h12.*f1;


v1 = 4*alpha*(rho^2)*Pt*h12.*f1.*f2;
v2 = 4*alpha*(rho^2)*Pt*h12.*f2.*f1;
v12 = v1+z1;
v21 = v2+z2;
%v12 = 4*alpha*(rho^2)*Pt*h12*f1.*f2+z1;
%v21 = 4*alpha*(rho^2)*Pt*h12*f2.*f1+z2;
h(v12)
h(v21)
h([v12,v21]) 
c_theo = 1/2*log2(1+(var(v1)/(var(z1)+var(z2)+(var(z1)*var(z2)/var(v1)))));
c_real = mi(v12,v21);
% 
% sigma_h_1e = sqrt(5.4);
% h_1e = normrnd(0,sigma_h_1e,[1 1000000])';
% sigma_h_2e = sqrt(5.4);
% h_2e = normrnd(0,sigma_h_2e,[1 1000000])';
% 
% v_e = 4*alpha^2*Pt*h_1e.*h_2e.*f1.*f2;
% c_sk = cmi(v12,v21,v_e);


zeta = 0.8;
e_theo  = 1/2*zeta*(1-rho)^2*(var(f1)+var(f2))*Pt;
e_real = 1/2*zeta*(1-rho)^2*(var(b12+f1)+var(b21+f2))*Pt;