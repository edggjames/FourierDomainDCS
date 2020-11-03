function [S1_pos] = CosineFDmodelPositiveOmega(A,B,F,R,w,plot_fig)
% Function to implement analytical model of Cosine Fourier transform of
% G_1(\tau,R) into Fourier domain. All lengths are in cm, all times are in
% seconds. Returns real component of complex Fourier transform

% function is dependent on erfz function, which calculates error function 
% in the complex plane

% Inputs
    % 1) A = 3*mu_a - scalar
    % 2) B = mu_s_p - scalar
    % 3) F = k_0^2*6*alpha_Db - scalar
    % 4) R = either r_1 or r_2 in semi-infinite geometry model - scalar
    % 5) w = vector of frequency values in rads/second. For positive omega
    %        only
    % 6) plot_fig = plots real and imaginary parts of FT side by side if
    %               true
    
% Outputs
    % 1) S1 = Un-normalised first order spectra of S1
    
% define repeating values to be used in terms below
const_1 = (1i*B^2*F*R^2)./(2*w);
% vector, length omega, complex
const_2 = (1i*A*w)/(B*F);
% vector, length omega, complex
const_3 = ((1+1i)*B^2*sqrt(F)*R+(2-2i)*sqrt(A*B/F)*w)./(2*sqrt(2)*B*sqrt(w));
% vector, length omega, complex
const_4 = ((1-1i)*B^2*sqrt(F)*R+(2+2i)*sqrt(A*B/F)*w)./(2*sqrt(2)*B*sqrt(w));
% vector, length omega, complex

% define all terms in expression
term_1 = (1./w.^(3/2)).*B.*exp(-const_1/2-const_2).*F^(1/2)*(1/4+1i/4);
% vector, length omega, complex
term_2 = -1i*exp(const_1).*(1-erfz(const_3));
% vector, length omega, complex
term_3 = exp(2*const_2).*(1-erfz(const_4));
% vector, length omega, complex

% combine all terms and extract real component
S1_pos = term_1.*(term_2+term_3);
% vector, length omega, complex

if plot_fig
    % plot real and imaginary components side by side
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)
    loglog(w,real(S1_pos),'r')
    title('Real Component of FD model')
    subplot(1,2,2)
    loglog(w,abs(imag(S1_pos)),'r')
    title('Absolute Imaginary Component of FD model')
end

S1_pos = real(S1_pos);
% vector, length omega, real

end

