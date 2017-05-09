function dydt = osc2cell(t,y)
%   global I1 I2 c e gs21 gs12 E1 E2
%   y(1) : v1, y(2) : w1, y(3) : s12
%   y(4) : v2, y(5) : w2, y(6) : s21
    e = 0.01; c = -0.1;
    E1 = 5; E2 = 5;
    I1 = 0; I2 = 0.12;
    A = 3; B = 3;
    
    gs12 = 0; gs21 = 0.01; 
 
    vth = 0.3; vsl = 0.001;
    
    % a few defined equations
    f = @(v,w) v*(v+c)*(1-v)-w;
    g = @(v,w) v-0.5*w;
    N = @(v) 0.5*( (1+ atan((v-vth)/vsl) ) );
    
    % equations
    dydt = [(f(y(1),y(2))+I1-(gs21*y(6)*(y(1)-E1)))/e   % v1'  = (f(v1,w1) + I1 - gs21*s21*(v1-E1))/e
    g(y(1),y(2))                                        % w1'  = g(v1,w1)    
    A*(N(y(1)))*(1-y(3))-B*y(3)                         % s12' = A*N(v1)*(1-s12) - B*s12
    (f(y(4),y(5))+I2-gs12*y(3)*(y(4)-E2))/e          % v2'  = (f(v2,w2)+I2 - gs12*s12*(v2-E2))/e 
    g(y(4),y(5))                                        % w2'  = g(v2,w2) 
    A*(N(y(4)))*(1-y(6))-B*y(6)                         % s21' = A*N(v2)*(1-s21) - B*s21
    ];       
end