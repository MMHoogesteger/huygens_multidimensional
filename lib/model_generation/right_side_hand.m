function RHS = right_side_hand(mdl, Ux, Uy, Uphi, Utheta, MMT, MR, MR2, R, gamma, x, dx, y, dy, phi, dphi, theta, dtheta, psi, c, ct, k, kt, d, g)

RHS = sym(zeros(3+mdl.N,1));

RHS(1) = -MMT*sum(dphi.^2.*R.*cos(gamma+phi))...
    - MR*sum(dtheta.^2.*sin(theta).*cos(psi+phi)...
    + 2.*dtheta.*dphi.*cos(theta).*sin(psi+phi)...
    + dphi.^2.*sin(theta).*cos(psi+phi))...
    + c*dx + k*x - Ux;

RHS(2) = -MMT*sum(dphi.^2.*R.*sin(gamma+phi))...
    + MR*sum(-dtheta.^2.*sin(theta).*sin(psi+phi)...
    + 2.*dtheta.*dphi.*cos(theta).*cos(psi+phi)...
    - dphi.^2.*sin(theta).*sin(psi+phi))...
    + c*dy + k*y - Uy;

RHS(3) = MR*sum(-dtheta.^2.*R.*sin(psi-gamma).*sin(theta)...
    + 2.*dphi.*dtheta.*R.*cos(psi-gamma).*cos(theta))...
    + 2*MR2.*sum(dphi.*dtheta.*sin(theta).*cos(theta))...
    + ct*dphi + kt*phi - Uphi;

for met_k = 1:mdl.N
    theta_k = theta(met_k);
    gamma_k = gamma(met_k);
    psi_k = psi(met_k);
    R_k = R(met_k);
    RHS(3+met_k) = -MR*dphi^2*R_k*cos(psi_k-gamma_k)*cos(theta_k)...
        - MR2*dphi^2*sin(theta_k)*cos(theta_k)...
        + d*dtheta(met_k) - MR*g*sin(theta_k) - Utheta(met_k);
end

RHS = -RHS([mdl.DOFs 3+(1:mdl.N)],:);


end