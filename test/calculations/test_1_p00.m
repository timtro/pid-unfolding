% http://ctms.engin.umich.edu/CTMS/index.php?example=Introduction&section=ControlPID

m = 1; b = 10; k = 20; F = 1;

s = tf('s');
P = 1/(s^2 + b*s + k);
C = pid(300, 0, 0, 0.01);
sys = feedback(C * P, 1);

stepinfo(sys)

t = 0:0.01:2;
step(sys, t)