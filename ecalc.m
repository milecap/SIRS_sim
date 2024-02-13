function ep = ecalc(B0,t,e,d)

gamma = 4258; %Hz/G
w0 = 2*pi*gamma*B0;
% ex = e(1)*sin(w0*t) + e(2)*cos(w0*t);
% ey = e(1)*cos(w0*t) - e(2)*sin(w0*t);

ep = Rot('z',w0*t,d)*e;

% ep = [ex, ey, e(3)];