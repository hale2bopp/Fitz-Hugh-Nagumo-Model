% 2 cell FitzHugh Nagumo model
%function [T,Y] = twocellFHmodel()

tspan = [0 10];

% system 1 - v1,w1 and s21
% initial conditions
v1_0 =  0;
w1_0 =  0;
s12_0 = 0;
v2_0 =  0;
w2_0 =  0;
s21_0 = 0;

%    global I1 I2 c e gs21 gs12 E1 E2
[T,Y] = ode15s(@osc2cell,tspan,[v1_0 w1_0 s12_0 v2_0 w2_0 s21_0]);
length(T)
% Plot voltage vs time
figure(1)
grid on
plot(T,Y(:,1),'r');
xlabel('Time(s)');
ylabel('Voltage cell 1(Red) and Voltage cell 2 (Blue)');
title('Voltage vs Time-Two Cell FitzHugh Nagumo');

hold on
grid on
plot(T,Y(:,4),'b');
%{
hold on
grid on
plot(T,Y(:,3),'k')
%}

hold off
legend('voltage cell1','voltage cell2','s12')
%end



