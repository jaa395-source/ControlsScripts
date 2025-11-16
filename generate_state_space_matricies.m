function [a_mat,b_mat,c_mat,d_mat] = generate_state_space_matricies(M, C1, K0, K2, v0)
a_mat = [0 0 1 0;...
    0 0 0 1;...
    -inv(M)*(K0 + K2*v0^2) -inv(M)*C1*v0];

b_mat = [0; 0; inv(M)*[0;1]];

c_mat = [0 1 0 0];

d_mat = [0];
end