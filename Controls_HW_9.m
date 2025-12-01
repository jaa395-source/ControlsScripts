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
legend_string = {};


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
    legend_string = horzcat(legend_string, "Poles @ [" + string(curr_poles(1)) + ", " + string(curr_poles(2)) + "], GM = " + string(Gm_ss) + ", PM = " + string(PM_ss));
end
legend(legend_string);
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

legend_string = {};

for field = 1:length(struct_fields)
    figure(z_fig);
    plot(z.(struct_fields(field)))

    figure(u_fig);
    plot(u.(struct_fields(field)))


    legend_string = horzcat(legend_string, "Poles @ [" + string(k_desired_sets(field, 1)) + ", " + string(k_desired_sets(field, 2)) + "]");
end

figure(z_fig);
legend(legend_string);

figure(u_fig);
legend(legend_string);

%% Problem 2

%2a
v0 = 5;
[M, C, K0, K2] = matricies();
[a_mat,b_mat,c_mat,d_mat] = generate_state_space_matricies(M, C, K0, K2, v0);
Q = [0 0 0 0;
     0 1 0 0;
     0 0 0 0;
     0 0 0 0];
sys = ss(a_mat, b_mat, c_mat, d_mat);
R = 1;
[K,S,P] = lqr(sys,Q,R);
eig(a_mat - b_mat*K)
