% Test A
% Comaprison to P-control example at:
%   http://ctms.engin.umich.edu/CTMS/index.php?example=Introduction&section=ControlPID

m = 1; b = 10; k = 20; F = 1;

s = tf('s');
P = 1/(s^2 + b*s + k);
C = pid(300, 0, 0);
sys = feedback(C * P, 1);

stepinfo(sys)

t = 0:0.01:2;
step(sys, t)

stepsys = sys * 1/s

syms s t
[stepsys_num, stepsys_den] = tfdata(stepsys);
stepsys_sym = poly2sym(cell2mat(stepsys_num), s) / poly2sym(cell2mat(stepsys_den), s);

ilaplace(stepsys_sym)
