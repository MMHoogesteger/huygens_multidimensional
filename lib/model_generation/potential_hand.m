function V = potential_hand(x, y, phi, theta, k, kt, g,MR)

V = MR*g*sum(cos(theta))+0.5*k*(x^2+y^2)+0.5*kt*phi^2;


end