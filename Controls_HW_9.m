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
    z.("poles_" + string(j)) = initial(feedback(ss(A-B*K-L*C,L,K,0)*ss(A,B,C,D), 1), [1;0;0;0]);
    u.("poles_" + string(j))= initial(feedback(ss(A-B*K-L*C,L,K,0),ss(A,B,C,D)), [1;0;0;0]);
    bode(ss(A-B*K-L*C,L,K,0)*ss(A,B,C,D));
    [Gm_ss, PM_ss] = margin(ss_tf);
    legend_string = horzcat(legend_string, "Poles @ [" + string(curr_poles(1)) + ", " + string(curr_poles(2)) + "], GM = " + string(Gm_ss) + ", PM = " + string(PM_ss));
end
legend(legend_string);
hold off;

struct_fields = string(fieldnames(z));

z_fig
u_fig = figure();


