close all; clear; clc;
%% Problem 1
A = [0 1; 0 0];
B = [0; 1];
C = [1 0];
D = [0];

state = ss(A, B, C, D);

k_desired_sets = [-1+1i, -1-1i; -sqrt(2), -sqrt(2); -5+1i -5-1i];
l_desired_set = [-10, -10];
z = struct();
u = struct();
bode_legend_string = {};


figure;
hold on;
for j = 1:length(k_desired_sets)

    curr_poles = k_desired_sets(j,:);
    K = acker(A, B, curr_poles);
    L = acker(A', C', l_desired_set)';

    At = [A-B*K            B*K;
        zeros(size(A))    A-L*C];

    Bt = [B;
        zeros(size(B))];

    Ct = [C zeros(size(C))];

    new_sys = ss(At, Bt, Ct, D);

    k_r = -1/(Ct*inv(At)*Bt);
    final_tf = k_r*new_sys;
    z.("poles_" + string(j)) = initial(feedback(ss(A,B,C,D)*ss(A-B*K-L*C,L,K,0), 1), [1;0;0;0]);

    u.("poles_" + string(j))= initial(feedback(ss(A-B*K-L*C,L,K,0),ss(A,B,C,D)), [1;0;0;0]);
    bode(ss(A-B*K-L*C,L,K,0)*ss(A,B,C,D));
    [Gm_ss, PM_ss] = margin(final_tf);
    bode_legend_string = horzcat(bode_legend_string, "Poles @ [" + string(curr_poles(1)) + ", " + string(curr_poles(2)) + "], GM = " + string(Gm_ss) + ", PM = " + string(PM_ss));
end
legend(bode_legend_string);
hold off;

%%
struct_fields = string(fieldnames(z));

z_fig = figure();
hold on;
title("Response to Initial Conditions");
xlabel("Time (seconds)");
ylabel("Amplitude");

u_fig = figure();
hold on;
title("Response to Initial Conditions");
xlabel("Time (seconds)");
ylabel("Amplitude");

bode_legend_string = {};

for field = 1:length(struct_fields)
    figure(z_fig);
    plot(z.(struct_fields(field)))

    figure(u_fig);
    plot(u.(struct_fields(field)))


    bode_legend_string = horzcat(bode_legend_string, "Poles @ [" + string(k_desired_sets(field, 1)) + ", " + string(k_desired_sets(field, 2)) + "]");
end

figure(z_fig);
legend(bode_legend_string);

figure(u_fig);
legend(bode_legend_string);

%% Problem 2

%2a
v0 = 5;
[M, C, K0, K2] = matricies();
[a_mat,b_mat,c_mat,d_mat] = generate_state_space_matricies(M, C, K0, K2, v0);
sys = ss(a_mat, b_mat, c_mat, d_mat);
Q = [0 0 0 0;
    0 1 0 0;
    0 0 0 0;
    0 0 0 0];
R = 1;
[K,S,P] = lqr(sys,Q,R);
eig(a_mat - b_mat*K)

% 2b
sigma_d = 0.1;
sigma_n = 10^-3;
mu_d = 0;
mu_n = 0;

sigma_d = sigma_d^2 - mu_d^2;
sigma_n = sigma_n^2 - mu_n^2;

[kalmf,L,P] = kalman(sys, sigma_d, sigma_n, 0);
eig(a_mat - L*c_mat)

% 2c
controller_d_mat = 0;
controller_ss = ss(a_mat-b_mat*K-L*c_mat, L, -K, controller_d_mat);
figure;
rlocus(controller_ss*sys);
figure;
margin(controller_ss*sys);

%% Problem 3

% a
s = tf('s');
t = 0:0.01:150;
u = 0.002*ones(size(t));

P_tf = tf(sys);
C_prime = 5; % from sisotool
P_prime_tf = feedback(P_tf, C_prime);
z = 2.1797;
p = [0 -1];
k = 0.6882;
C_final = zpk(z, p, k);
P_final = P_prime_tf;

td = [0.5 5 50]*10^-3;

step_legend_string = {};
bode_legend_string = {};

step_fig = figure();
hold on;
bode_fig = figure();
hold on;


for i = 1:length(td)
    time_delay = exp(-s*td(i));
    open_loop_tf = P_final*C_final*time_delay;
    closed_loop_tf = feedback(open_loop_tf, 1);

    figure(step_fig);
    lsim(closed_loop_tf, u, t);

    figure(bode_fig);
    bode(open_loop_tf);
    [Gm, Pm] = margin(open_loop_tf);

    step_legend_string = horzcat(step_legend_string, "Time Delay = " + string(td(i)) + " ms");
    bode_legend_string = horzcat(bode_legend_string, "Time Delay = " + string(td(i)) + " ms, GM = " + string(Gm) + ", PM = " + string(Pm));

end

figure(step_fig);
title("Continuous-Time Closed Loop Simulated Step Response with Time Delay");
legend(step_legend_string);

figure(bode_fig);
legend(bode_legend_string);
title("Continuous-Time Closed Loop Bode with Time Delay");

%% b/c
step_fig = figure();
hold on;
bode_fig = figure();
hold on;


step_legend_string = {};
bode_legend_string = {};


Fs = [1e3 100 10];
for i = 1:length(Fs)
    F = Fs(i);
    Ti = 1/F;

    Cd = c2d(C_final, Ti, 'tustin');
    Pd = c2d(P_final, Ti, 'tustin');

    t = 0:Ti:5;
    u = 0.002*ones(size(t));
    discrete_ol = (Cd)*(Pd);
    discrete_cl = feedback(discrete_ol, 1);

    figure(step_fig);
    dlsim(tf(discrete_cl), u, t);

    figure(bode_fig);
    bode(discrete_ol);
    [Gm, Pm] = margin(discrete_ol);

    step_legend_string = horzcat(step_legend_string, "Time Delay = " + string(td(i)) + " ms");
    bode_legend_string = horzcat(bode_legend_string, "Time Delay = " + string(td(i)) + " ms, GM = " + string(Gm) + ", PM = " + string(Pm));
end


figure(step_fig);
title("Discrete-Time Closed Loop Simulated Step Response");
legend(step_legend_string);

figure(bode_fig);
legend(bode_legend_string);
title("Discrete-Time Closed Loop Bode");