function [G,g,f1,f2] = oneOver_n_OctBandFilter(N, n, fs, low_f, high_f)
% [G,g,f1,f2] = oneOver_n_OctBandFilter(N, n, fs, low_f, high_f)
% Creates fractional octave band filters starting from low_f and ending at
% high_f
% 
% USAGE:
% 
% G: frequency domain filters
% g: time domain filters
% f1: frequency band limit, low
% f2: frequency band limit, high
% 
% N: number of fft bins
% n: 1 / n -octave band filter design
% fs: sampling frequency
% low_f: lowest frequency,  first freq. band will be (0, low_f)
% high_f: highest frequency, last freq. band will be (h_f fs/2]
%         where h_f is the closest octave band cut-off freq to high_f
% I, number of filter bands
% 
% EXAMPLES:
% 
% 
% N = 2^16;
% n = 3; % One-third-octave bands, useful for, e.g., room acoustics
% % analysis
% fs = 48000;
% low_f = 62.5; % lowest band central freq 62.5
% high_f = 20000;
% oneOver_n_OctBandFilter(N, n, fs, low_f, high_f)
% 
% N = 2^16;
% n = 12; % One-twelth-octave bands, central frequencies are piano key
% % frequencies, useful for example for the analysis of music
% fs = 48000;
% low_f = 62.5; % lowest band central freq 62.5
% high_f = 20000;
% oneOver_n_OctBandFilter(N, n, fs, low_f, high_f)
% 
% 
% More info on the design of the FIRs is found from:
% [1] "Orthogonal-like fractional-octave-band filters"
% J. Antoni, The Journal of the Acoustical Society of America, 2010


% SDM toolbox : oneOver_n_OctBandFilter
% Sakari Tervo & Jukka Pätynen, Aalto University, 2011-2016
% Sakari.Tervo@aalto.fi and Jukka.Patynen@aalto.fi

if nargin < 2
   fs = 48e3;
   n = 3;
   low_f = 62.5;
end

% n = 3;

w_s = fs*2*pi;

I = 10*n; % This is the number of fractional octave bands that are created

% Mid-frequencies of the freq. band
w_m = zeros(I,1);
% First one is low_f
w_m(1) = low_f*2*pi;%(2^(1/(2*n))); 
for i = 2:I+1
    w_m(i) = 2^(1/n)*w_m(i-1);
end

w_1 = w_m/(2^(1/(2*n)));
% f1 and f2 are the initial low and high cut-off frequencies from
% definitions of the user
% w_2 = w_m*(2^(1/(2*n)));
% f1 = w_1/(2*pi); 
% f2 = w_2/(2*pi);

k = round(N*w_1/w_s); % k_i
km = round(N*w_m/w_s); % k_i
% k2 = round(N*w_2/w_s);

G_1 = zeros(N,I);
G_2 = zeros(N,I);
G = zeros(N,I+1);

for i = 1:I+1;% filter band index
    clear phi
    P = round(N*(w_m(i)-w_1(i))/w_s)+1;
    p = -4*P:1:2*P;

    phi(1,:) = 1/2*(p./P + 1);
    L = 6;

    for l = 2:L
        phi(l,:) = sin(pi/2*phi(l-1,:));
    end
    phi = phi(L,:);
   
    % According to [1]
    %a = (sin((pi/2)*phi));
    %b = (cos((pi/2)*phi));
    
    % Proposed version
    a = (sin((pi/2)*phi)+1)/2;
    b = abs((1-a));

    ind = find(a > 0.5,1,'First');
    ind = ind + 0;
    p = -P:1:2*P;

    a = a(ind+p);
    b = b(ind+p);

    k_vect = (k(i)+0) + p;
    inds = ~ (k_vect < 1 | k_vect > N);
    
    % Increasing part
    G_1(k_vect(inds),i+1) = a(inds);
    % Decreasing part
    G_2(k_vect(inds),i) = b(inds);
end

for i = 1:I+1
    if i == I+1;
        ind1 = find(G_1(:,i)> 0,1,'Last');
         G(1:end,i) = G_1(:,i);
         G(ind1+1:end,i) = 1;
    elseif i == 1
        ind1 = find(G_2(:,i)> 0,1,'First');
         G(1:end,i) = G_2(:,i);
         G(1:ind1,i) = 1;
    else
        inds = G_1(:,i) < G_2(:,i);
        ind1 = find(inds(:)> 0,1,'First');
        inds(1:ind1) = 1;
        
        % ind1 = find(inds(:) < 1,1,'First');
        % inds(1:ind1) = 1;
        ind1 = min(km(i-1),size(G,1));
        G(ind1+1:end,i) = G_2(ind1+1:end,i);
        G(1:ind1,i) = G_1(1:ind1,i);
    end
end

% Find the true cross-over frequencies
kc = zeros(I,1);
for i = 1:I
    ind1 = find(G(:,i+1) > G(:,i),1,'First');;
    if ~isempty(ind1)
        kc(i) = ind1;
    else
        kc(i) = size(G,1);
    end
end
kc = kc-2;

% Update the frequency limits
f1 = [0; kc]/N*fs;
f2 = [kc; fs/2*N/fs]/N*fs;

% Limit the frequency bands according to high_f
inds = find(f2 > high_f);
if ~isempty(inds)
    G(:, inds(1)) = sum(G(:,inds),2);
    G(:,inds(1)+1:end) = [];
end

indsf = f1 < high_f;
if ~isempty(indsf)
    f1 = f1(indsf);
    f2 = f2(indsf);
    kc = kc(indsf(1:length(kc)));
end

% Shift to center delayed time domain version
G = [G(1:end/2+1,:);conj(G(end/2:-1:2,:))];
g = real(fftshift2(ifft(G)));

if nargout < 1
    figure(1)
    semilogx(((0:N-1)/N*fs),20*log10(G))
    hold on
    semilogx(([kc kc]/N*fs),[0 1],'b-')
    hold off
end


