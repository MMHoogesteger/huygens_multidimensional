function c0 = fourierDC(t,f,omega)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dt = t(end)-t(1);
nP = floor(dt/(2*pi/omega));
t1 = t(1);
tP = t1+nP*2*pi/omega;
idL = find(t<tP, 1, 'last' );
 
fP = interp1(t,f,tP);
t = [t(1:idL)];
f = [f(1:idL)];



c0 = omega/2/pi/nP*trapz(f)*(t(2)-t(1));


end

