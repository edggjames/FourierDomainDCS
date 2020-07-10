function [IRF_neg,IRF_pos,IRF_neg_norm,IRF_pos_norm,IRF_neg_phase,...
    IRF_pos_phase] = ...
    IRF(w_s,tau_e,n,freq_range,const1,const2)
% Function to implement IRF for n-phase detection for certain exposure
% time, frame rate and range of detuning frequency offsets. Based on
% equation 26 from 'Spatiotemporal Heterodyne Detection' Atlan and Gross
% 2007. 

% Inputs
% 1) w_s - camera frame rate (rads/s)
% 2) tau_e - camera exposure time (seconds)
% 3) n - number of images
% 4) freq_range - range of detuning frequencies (Hz)
% 5) const1 - factor to divide tau_e by (should be 2)
% 6) const2 - factor to multiply camera frame rate by (should be 1)

% Outputs
% 1) IRF_neg - negative IRF response (mask 1)
% 2) IRF_pos - positive IRF response (mask 2)
% 3) IRF_neg_norm - IRF_neg normalised by maximum value (mask 1)
% 4) IRF_pos_norm - IRF_pos normalised by maximum value (mask 2)

BPF = sinc(freq_range*tau_e/(const1/2)); 

SUM_neg = 0;
SUM_pos = 0;

w_s = w_s * const2; 

w_detune = freq_range*2*pi;

for k = 1:n
    
    % first term is roots of unity, second is due to frequency of oscillation
    neg = exp(-2*1i*k*pi/n)*exp(-2*1i*k*pi*w_detune/w_s);
    SUM_neg = SUM_neg + neg; 
    
    pos = exp(-2*1i*k*pi/n)*exp(+2*1i*k*pi*w_detune/w_s);
    SUM_pos = SUM_pos + pos; 

end

IRF_neg = BPF .* SUM_neg;
IRF_pos = BPF .* SUM_pos;

% calculate phase of IRF
IRF_neg_phase = atan2(imag(IRF_neg),real(IRF_neg));
IRF_pos_phase = atan2(imag(IRF_pos),real(IRF_pos));

% take square of absolute values
IRF_neg = abs(IRF_neg).^2;
IRF_pos = abs(IRF_pos).^2;

% normalise
IRF_neg_norm = IRF_neg./max(abs(IRF_neg));
IRF_pos_norm = IRF_pos./max(abs(IRF_pos));
end