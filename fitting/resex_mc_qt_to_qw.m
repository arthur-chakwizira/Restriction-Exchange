function [Qw, w] = resex_mc_qt_to_qw(qt, dt)
%computes encoding power spectrum for one waveform
%also returns the frequency grid

if ~isrow(qt); qt = qt'; end
% q = [qt zeros(size(qt)) zeros(size(qt))];
% q = [qt zeros(size(qt)) zeros(size(qt))  zeros(size(qt))  zeros(size(qt)) ];
% q = qt;
len = numel(qt)*30;
qw = (fftshift(fft(qt, len)))*dt;



Qw = (abs(qw)).^2;


Nw = numel(qw);
df = 1/(Nw*dt);
if mod(Nw, 2) == 0; offset = 1; else; offset = 0; end
w = 2*pi*df*(-floor(Nw/2):floor(Nw/2)-offset); %omega
end