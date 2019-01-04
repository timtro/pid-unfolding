% Test D
% Comparison to PID-control example from:
%   http://ctms.engin.umich.edu/CTMS/index.php?example=Introduction&section=ControlPID

m = 1; b = 10; k = 20; F = 1;

s = tf('s');
P = 1/(s^2 + b*s + k);
C = pid(350, 300, 50);
sys = feedback(C * P, 1);

stepinfo(sys)

tvec = 0:0.01:2;
step(sys, tvec)

stepsys = sys * 1/s

syms s t
[stepsys_num, stepsys_den] = tfdata(stepsys);
stepsys_sym = poly2sym(cell2mat(stepsys_num), s) / poly2sym(cell2mat(stepsys_den), s);

ilaplace(stepsys_sym)
