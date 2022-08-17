N = 3;
syms MT MMT MR MR2 JT phi

theta = sym('theta',[N,1]);
dtheta = sym('dtheta',[N,1]);
psi = sym('psi',[N,1]);
gamma = sym('gamma',[N,1]);
R = sym('R',[N,1]);

Mh = sym(zeros(3+N,3+N));

Mh(1,1) = MT;
Mh(2,2) = MT;
Mh(1,3) = -sum(MR.*sin(theta).*sin(psi+phi) + MMT.*R.*sin(gamma+phi));
Mh(3,1) = Mh(1,3);

Mh(2,3) = sum(MR*sin(theta).*cos(psi+phi) + MMT.*R.*cos(gamma+phi));
Mh(3,2) = Mh(2,3);

Mh(3,3) = sum(MMT.*R.^2 + 2.*MR.*R.*cos(psi-gamma).*sin(theta) + MR2.*sin(theta).^2) + JT;

for met_k = 1:N
    Mh(1,3+met_k) = MR*cos(theta(met_k))*cos(psi(met_k)+phi);
    Mh(3+met_k,1) = Mh(1,3+met_k);

    Mh(2,3+met_k) = MR*cos(theta(met_k))*sin(psi(met_k)+phi);
    Mh(3+met_k,2) = Mh(2,3+met_k);

    Mh(3,3+met_k) = MR.*R(met_k).*sin(psi(met_k)-gamma(met_k)).*cos(theta(met_k));
    Mh(3+met_k,3) = Mh(3,3+met_k);
    Mh(3+met_k,3+met_k) = MR2;
end

MM = Mh;

spy(inv(MM))