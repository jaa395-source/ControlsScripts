% MAE 4780/5780 LAB 4 Analysis Template
% Feedback control systems
% Original from Lab 2, Fall 2019 and Lab 4, Fall 2021

clc; clear all;close all;

% fill in your values for alpha from your lab analysis if the ones that are
% close to the ones below, which are taken from the Lab 3 prelab (70% of
% the theoretical alpha1 = 1.5 and alpha2 = -10; otherwise go ahead and use
% these values

alpha_1 = 1.05; 
alpha_2 = -7;
[A, B, C, D] = create_ip_model(alpha_1, alpha_2);

% Fill in your choice of eigenvalues for A-BK and A-LC
all_ev = ???; % need a vector here
all_ev_obs = ???;  % e.g., all_ev_obs = 3*all_ev;

% Generate K matrix
K = place(A, B, all_ev);

% Generate L matrix
L = place(A', C', all_ev_obs)';

% save the controller

save_ss_control(alpha_1, alpha_2, K, L, 'my_ss_control_SIMO'); 

%% below are functions that you don't need to modify
% inverted pendulum model
% function for generating the matrices for the inverted pendulum system
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
    
    C = [1.0, 0.0, 0.0, 0.0; 0.0 1.0 0.0 0.0];
    
    D = [0.0; 0.0];    
end

% function for saving control values
function [] = save_ss_control(alpha_1, alpha_2, K_, L_, save_name)
    K = K_;
    L = L_;
    
    assert(alpha_1 < 2.7, "alpha_1 is too large")
    assert(0.9 < alpha_1, "alpha_1 is too small")
    assert(abs(alpha_2) < 20.0, "the magnitude of alpha_2 is too large") 
    assert(5.0 < abs(alpha_2), "the magnitude of alpha_2 is too small") 
    assert(sign(alpha_2) == -1, "alpha_2 needs to be negative")
    
    [A, B, C, ~] = create_ip_model(alpha_1, alpha_2);
    
    verify_process_matrices()
    verify_controller_observer(A, B, C, K, L)
    
    save(save_name, 'A', 'B', 'C', 'K', 'L', '-v6');
    disp(['controller successfully saved to file: ', save_name, '.mat'])
end

function [] = verify_controller_observer(A, B, C, K, L)
    assert(all(size(L) == [4, 2]), "L must be 4x2")
    assert(all(size(K) == [1, 4]), "K must be 1x4")
end

function [] = verify_process_matrices()
    alpha_verify_1 = 0.5; 
    alpha_verify_2 = 2.0;
    [A, B, C, D] = create_ip_model(alpha_verify_1, alpha_verify_2);

    assert(all(size(A) == [4, 4]), "A must be 4x4")
    assert(all(size(B) == [4, 1]), "B must be 4x1")
    assert(all(size(C) == [2, 4]), "C must be 2x4")
    assert(all(size(D) == [2, 1]), "D must be 2x1")
end 
