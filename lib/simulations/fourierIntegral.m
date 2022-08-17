function ck = fourierIntegral(t,f,omega,kmax)
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


k = 0;
c0 = omega/2/pi/nP*trapz(exp(-1i*omega*k*t).*f)*(t(2)-t(1));
ck = zeros(1,kmax);
for k = 1:kmax
    ck(k) = omega/2/pi/nP*trapz(exp(-1i*omega*k*t).*f)*(t(2)-t(1));
end

end

