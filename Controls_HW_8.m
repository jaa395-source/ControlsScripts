%% Problem 1
% x = [phi delta phi_dot delta_dot];
% xdot = [phi_dot delta_dot phi_ddot delta_ddot];
% xdot = Ax + Bu
% y = delta = Cx + Du = [0 1 0 0]x + D[0]
[M, C1, K0, K2] = matricies();
v0 = 5;
a_mat = [0 0 1 0;...
    0 0 0 1;...
    -inv(M)*(K0 + K2*v0^2) -inv(M)*C1*v0];

b_mat = [0; 0; inv(M)*[0;1]];

c_mat = [0 1 0 0];

d_mat = [0];

%% 1a
P = ss(a_mat, b_mat, c_mat, d_mat);
P_tf = tf(P);
figure;
bode(P)

%% 1b
C_prime = 5; % from sisotool
P_prime_tf = feedback(P_tf, C_prime);
z = 2.1797;
p = [0 -1];
k = 0.6882;
C_final = zpk(z, p, k);
P_final = P_prime_tf;

%% 1c
v0_values = [4, 6];
figure;
hold on;
bode(P_final*C_final);

% Bodes
for i = 1:length(v0_values)

    v0_curr = v0_values(i);
    a_mat = [0 0 1 0;...
        0 0 0 1;...
        -inv(M)*(K0 + K2*v0_curr^2) -inv(M)*C1*v0_curr];

    b_mat = [0; 0; inv(M)*[0;1]];

    c_mat = [0 1 0 0];

    d_mat = [0];

    new_P = ss(a_mat, b_mat, c_mat, d_mat);
    new_P_tf = tf(new_P);
    bode(new_P_tf*C_final);
end

legend("v0 = 5 m/s", "v0 = 4 m/s", "v0 = 6 m/s");
hold off;

% Root Loci
figure;
hold on;
rlocus(P_final*C_final);


for i = 1:length(v0_values)

    v0_curr = v0_values(i);
    a_mat = [0 0 1 0;...
        0 0 0 1;...
        -inv(M)*(K0 + K2*v0_curr^2) -inv(M)*C1*v0_curr];

    b_mat = [0; 0; inv(M)*[0;1]];

    c_mat = [0 1 0 0];

    d_mat = [0];

    new_P = ss(a_mat, b_mat, c_mat, d_mat);
    new_P_tf = tf(new_P);
    rlocus(new_P_tf*C_final);
end

legend("v0 = 5 m/s", "v0 = 4 m/s", "v0 = 6 m/s");
hold off;

% From root loci, v = 6 has poles on the LHP making it still possible to
% just barely have a stable loop but for v = 4 those same poles are in the
% RHP so it is always unstable

%% 1d

v0_values = [4, 5, 6];
step_freq = 0.002;


for i = 1:length(v0_values)

    v0_curr = v0_values(i);
    a_mat = [0 0 1 0;...
        0 0 0 1;...
        -inv(M)*(K0 + K2*v0_curr^2) -inv(M)*C1*v0_curr];

    b_mat = [0; 0; inv(M)*[0;1]];

    c_mat = [0 1 0 0];

    d_mat = [0];

    new_P = ss(a_mat, b_mat, c_mat, d_mat);
    new_P_tf = tf(new_P);
    dc_gain = evalfr(new_P_tf, 0);
    dc_gain = 1/dc_gain;

    figure;
    step(new_P*dc_gain*step_freq)
    legend("v0 = " + v0_curr + " m/s")
end


%% Problem 2
clear i
[M, C1, K0, K2] = matricies();
v0 = 5;
a_mat = [0 0 1 0;...
    0 0 0 1;...
    -inv(M)*(K0 + K2*v0^2) -inv(M)*C1*v0];

b_mat = [0; 0; inv(M)*[0;1]];

c_mat = [0 1 0 0];

d_mat = [0];


% 2a
reachability_mat = [b_mat a_mat*b_mat a_mat^2*b_mat  a_mat^3*b_mat];
if rank(reachability_mat) == size(reachability_mat, 1)
    disp("This system is Reachable!")
end

% 2b
poles = [-2 -10 -1+i  -1-i];
ss_model = ss(a_mat, b_mat, c_mat, d_mat);
gains = place(a_mat, b_mat, poles);

new_a_mat = a_mat-b_mat*gains;
new_ss_model = ss(new_a_mat, b_mat, c_mat, d_mat);

% 2c
step_freq = 0.002;
new_ss_tf = tf(new_ss_model);
dc_gain = evalfr(new_ss_tf, 0);
dc_gain = 1/dc_gain;
new_ss_model = new_ss_model*dc_gain;
step(new_ss_model*step_freq);

%% Problem 3
clear i
[M, C1, K0, K2] = matricies();
v0 = 5;
a_mat = [0 0 1 0;...
    0 0 0 1;...
    -inv(M)*(K0 + K2*v0^2) -inv(M)*C1*v0];

b_mat = [0; 0; inv(M)*[0;1]];

c_mat = [0 1 0 0];

d_mat = [0];

% 3a
observability_mat = [c_mat; c_mat*a_mat; c_mat*a_mat^2;  c_mat*a_mat^3];
if rank(observability_mat) == size(observability_mat, 1)
    disp("This system is Observable!")
end

% 3b
poles = [-4 -20 -2+i -2-i];
ss_model = ss(a_mat, b_mat, c_mat, d_mat);
gains = place(a_mat, c_mat, poles);

new_a_mat = a_mat-c_mat*gains;
new_ss_model = ss(new_a_mat, b_mat, c_mat, d_mat);