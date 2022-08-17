% To get correct figures in Matlab, not important
set(groot, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

% PGFplots line width same as axis
set(groot, 'DefaultLineLineWidth', 0.4);

% Colors for background lines etc.
cGray = [0.3 0.3 0.3];
cBlack = [0 0 0];
cRed = [1 0 0];


set(groot, 'defaultAxesColorOrder', [cBlack; cRed]);