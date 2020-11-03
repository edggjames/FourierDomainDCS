function [S_1,s_1] = CosineFDmodelConvective(A,B,F,R1,R2,w)
% Function to determine FD model of G1 for both sides of semi-infinite
% geometry - normalised output for convective flow model. All lengths are
% in cm, all times are in seconds.

% Function is dependent on besselk function

% Inputs
    % 1) A = 3*mu_a - scalar
    % 2) B = mu_s_p - scalar
    % 3) F = k_0^2*V^2 - scalar
    % 4) R1 = r_1 in semi-infinite geometry model - scalar
    % 5) R2 = r_2 in semi-infinite geometry model - scalar
    % 6) w = vector of frequency values in rads/second. For positive omega
    %        only
    
% Outputs
    % 1) S_1 = un-normalised power spectra
    % 2) s_1 = normalised power spectra

term_1 = sqrt(((A*B)./(B^2*F*R1^2+w.^2))); 
term_2 = besselk(1,sqrt(A*(B*R1^2+w.^2/(B*F))));
term_3 = sqrt(((A*B)./(B^2*F*R2^2+w.^2))); 
term_4 = besselk(1,sqrt(A*(B*R2^2+w.^2/(B*F))));
term_5 = sqrt((B*F)/A)*R1*R2;
term_6 = R2*besselk(1,sqrt(A*B)*R1);
term_7 = R1*besselk(1,sqrt(A*B)*R2);

S_1 = (term_1.*term_2)-(term_3.*term_4);
s_1 = S_1*(term_5/(term_6-term_7));

end

