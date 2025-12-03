alpha_1 = 1.1722; 
alpha_2 = -7.9096;
[A, B, C, D] = create_ip_model(alpha_1, alpha_2);


t = 0:0.01:2;
u = zeros(size(t));
x0 = [0 0.01 0 0];

p1 = -10 + 10i;
p2 = -10 - 10i;
p3 = -50;
p4 = -30;
controller_poles = [p1, p2, p3, p4];


K = place(A, B, controller_poles);
L = place(A', C', 3*controller_poles)';

At = [ A-B*K             B*K
       zeros(size(A))    A-L*C ];




sys_cl = ss(A-B*K,B,C,0);

[y,t,x] = lsim(sys_cl,u,t,x0);
plot(t,y)
xlabel('Time (sec)')
ylabel('Ball Position (m)')





function [A, B, C, D] = create_ip_model(alpha_1, alpha_2)
    M_c = 0.915;
    M_r = 0.210;
    I_mot = 3.87e-7; 
    k = 3.7;
    R = 0.00635;
    L = 0.60;
    g = 9.8;

    H_11 = M_c + M_r + I_mot * k^2 / R^2;
    H_12 = -M_r * L * 0.5;
    H_21 = -1;
    H_22 = L * 0.66666666666667;

    H = [H_11 H_12; 
         H_21 H_22];
    C_q = [0.0 0.0;  
           0.0 g];
       
    C_v = [alpha_2 0.0;  
           0.0     0.0];
    
    B_u = [alpha_1; 
           0.0];
    iH = inv(H);
    
    A = [zeros(2,2), eye(2);  
         iH * C_q,   iH * C_v];
    
    B = [0.0; 
         0.0;
        iH * B_u];
    
    % Please update C here with [C_SISO; 0 0 0 0]
    C = [1 -0.4 0 0; 0 0 0 0];
    
    D = [0.0; 0.0];    
end