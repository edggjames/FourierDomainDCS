function [s1_pos] = CosineFDmodelPositiveOmegaBothR(A,B,F,R1,R2,w)
% Function to determine FD model of G1 for both sides of semi-infinite
% geometry - normalised output. All lengths are in cm, all times are in
% seconds. Returns real component of complex Fourier transform

% Function is dependent on erfz function, which calculates error function 
% in the complex plane, and also CosineFDmodelPositiveOmega function

% Inputs
    % 1) A = 3*mu_a - scalar
    % 2) B = mu_s_p - scalar
    % 3) F = k_0^2*6*alpha_Db - scalar
    % 4) R1 = r_1 in semi-infinite geometry model - scalar
    % 5) R1 = r_1 in semi-infinite geometry model - scalar
    % 6) w = vector of frequency values in rads/second. For positive omega
    %        only
    
% Outputs
    % 1) s1_pos = Normalised (by maximum value) first order spectra of S1

[S1_pos_R1] = CosineFDmodelPositiveOmega(A,B,F,R1,w,false);
[S1_pos_R2] = CosineFDmodelPositiveOmega(A,B,F,R2,w,false);
S1_pos = S1_pos_R1 - S1_pos_R2;
s1_pos = S1_pos./max(S1_pos);

end

