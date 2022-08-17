function Hm = hamiltonian_parts_hand(mdl,x, dx, y, dy, phi, dphi, theta, dtheta, mP, k, JP33, kt, MR2, MR, g)

Hx = 0.5.*mP.*dx.^2 + 0.5.*k.*x.^2;
Hy = 0.5.*mP.*dy.^2 + 0.5.*k.*y.^2;
Hphi = 0.5.*JP33.*dphi.^2 + 0.5.*kt.*phi.^2;
Htheta = 0.5*MR2.*dtheta.^2 + MR*g*cos(theta);
Hm = [Hx; Hy; Hphi; Htheta];

Hm = Hm([mdl.DOFs 3+(1:mdl.N)]);


end