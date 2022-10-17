clear; clc; close all
transversal = 0;

w0 = 1;              %Beam Waist
lambda = 632.8E-9;   %Wavelength
k = 2*pi/(lambda);   %Wavevector magnitude
kt = 20*w0;          %Transverse component of k
zr = k*(w0^2)/2;     %Rayleigh range
gamma = kt*w0/2;

m = 1;
L = 2*w0;

x = linspace(-2*w0, 2*w0, 500);
z = linspace(0,zr/8, 500); 
if transversal
    Z = 0; 
    [X, Y] = meshgrid(x);
else
    X = 0;
    [Z, Y] = meshgrid(z, x);
end
R = sqrt(X.^2 + Y.^2);
TH = angle(X + 1i*Y);
mu = 1 + 1i*Z/zr;  

J = besselj(m, kt*R./mu);
gb = (exp(1i*k*Z)./mu).*exp(-(X.^2+Y.^2)./(mu*w0^2));
W = J.*exp(1i*m*TH);
BG = exp(-1i*(Z*kt^2)./(2*k*mu)).*gb.*W;

figure(1)
pcolor(abs(BG).^2)
shading interp; axis equal; axis off; colormap(hot)