% test the transfer_function.m
clear;
close all;

Ts = 0.01;
num = [2,3,1];
den = [5,7,5,6];
system_tf = transfer_function(num,den,Ts);

sys = tf(num,den);
figure(1)
[A, B, C, D] = tf2ss(num,den)
step(sys, 10)
% [output_1, time_1] = step(sys, );

sim_time = 0.0;
time = [sim_time];
output = [system_tf.update(0.)];
while sim_time < 10.0
    u = 1;
    y = system_tf.update(u);
    sim_time = sim_time + Ts;
    
    time = cat(1, time, sim_time);
    output = cat(1, output, y);
end

figure(2);
plot(time, output);
title('使用matlab程序画的动态系统离散响应曲线');
grid on;






