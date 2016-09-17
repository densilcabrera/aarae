function OUT = imagesource_rectangular_room_demo(Lx,Ly,Lz,xs,ys,zs,xr,yr,zr,c,jitter,maxorder,ambiorder)
% image-source model for a rectangular room, by Densil Cabrera written for
% demonstration purposes in 2009, and adapted for AARAE in 2014.
%
% creates an impulse response of a rectangular room that has dimensions of
% Lx, Ly and Lz. the coordinate system's origin is in a corner of the room,
% and all coordinates within the room are positive (xs,ys,zs) are the
% coordinates of the source (xr,yr,zr) are the coordinates of the receiver
% If source and receiver are co-located, then the direct sound is not
% included in the impulse response. c is the speed of sound jitter
% introduces a random offset to each reflection time (and consequently
% amplitude) (try a value of 0.1 m). Jitter is not applied to the direct
% sound. maxorder is the maximum order of the image source calculations
% ambiorder is the ambisonics order (0 - 7). A value of 0 yields a single
% (omnidirectional) channel. A value of 7 yields 64 channels.
%
% Edit the absorption coefficients directly in the code (filtHzalpha). The
% absorption coefficient values are applied to all surfaces.
%
% The function calculates the time and amplitude of the highest order
% reflections first, and this is filtered by the wall reflection filter,
% before the next highest order reflections are added to it, and the wave
% filtered again. This continues until the direct sound is reached (which
% is not filtered by a wall reflection).
%
% a dissipation (air absorption) filter has not yet been implemented.
%
% This function has been adapted to use Nicolas Epain's HOAToolbox within
% the AARAE project for higher order Ambisonics encoding of the impulse
% response. Currently this makes the function VERY SLOW if high order
% reflections and high order ambisonics are generated - hopefully a
% speed-up will be added in a future revision.

if nargin < 13, ambiorder = 0; end % ambisonics order of 0 yields 1 omnidirectional channel
if nargin < 12, maxorder = 50; end % maximum reflection order calculated
if nargin < 11, jitter = 0; end % jitter standard deviation in metres (random offset to timing of reflections)
if nargin < 10, c = 344; end % speed of sound

if nargin < 9
    % default receiver position (room dimensions should be larger)
    xr = 1;
    yr = 1;
    zr = 1;
end

if nargin < 6
    % default source position (in the corner of the room)
    xs = 0;
    ys = 0;
    zs = 0;
end

if nargin < 3
    % input via dialog box
    % default room dimensions (Uni of Sydney reverberant room)
    Lx = 6.35;
    Ly = 5.1;
    Lz = 4;
    
    fs = 44100; % default sampling rate of generated wave
    
    % dialog box for settings
    prompt = {'Length, width and height of the room or box (m)', ...
        'Coordinates of the source (m)', ...
        'Coordinates of the receiver (m)', ...
        'Speed of sound (m/s)', ... 
        'Spatial jitter of reflections (m)', ...
        'Maximum order of reflections', ...
        'Ambisonics order', ...
        'Sampling rate (Hz)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {[num2str(Lx),', ',num2str(Ly),', ',num2str(Lz)], ...
        [num2str(xs),', ',num2str(ys),', ',num2str(zs)], ...
        [num2str(xr),', ',num2str(yr),', ',num2str(zr)], ...
        num2str(c), num2str(jitter), num2str(maxorder), ...
        num2str(ambiorder), num2str(fs)};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if isempty(answer)
        OUT = [];
        return
    else
        L = str2num(answer{1,1});
        Lx = L(1); Ly = L(2); Lz = L(3);
        s = str2num(answer{2,1});
        xs = s(1); ys = s(2); zs = s(3);
        r = str2num(answer{3,1});
        xr = s(1); yr = s(2); zr = s(3);
        c = str2num(answer{4,1});
        jitter = str2num(answer{5,1});
        maxorder = str2num(answer{6,1});
        ambiorder = str2num(answer{7,1});
        fs = str2num(answer{8,1});
    end
    
end



% number of channels
nchans = (ambiorder + 1)^2;
if nchans > 64, nchans = 64; end
hoaFmt = GenerateHoaFmt('res2d',ambiorder,'res3d',ambiorder);

% generate filter to simulate a single wall reflection. Edit the octave
% band absorption coefficents in the second column below
filtHzalpha = ...
    [0      1;
    31.5	0.01;
    63      0.01;
    125     0.01;
    250     0.01;
    500     0.02;
    1000	0.03;
    2000	0.05;
    4000	0.07;
    8000	0.09;
    16000	0.6
    fs/2    1];

% filterorder is one less than the number of 'a' coefficients
% (feedback taps, which are in the denominator) and 'b' coefficients (in the
% numerator). By increasing the filter order you can increase the accuracy
% of the filter's magnitude response, but this also affects the time and
% phase response of the filter.
filterorder = 6;

% Absorption coefficients (alpha) are converted to reflection coefficients
% (1-alpha), and then the square root is taken to convert to a pressure
% coefficient.
[b,a]=yulewalk(filterorder,filtHzalpha(:,1)./(0.5*fs),(1-filtHzalpha(:,2)).^0.5);

% pre-allocate zeros to output wave
out =  zeros(ceil(fs / c * max([Lx,Ly,Lz]) * (maxorder + 2)) , nchans);

% generally it is preferable not to use for loops in Matlab if there is a
% vectorized alternative availble. Vectorized code tends to use more memory
% but is faster to run. The below uses nested for loops anyway...
for o = 1:maxorder + 1
    order = maxorder + 1 - o; % reflection order
    for nx = -order:order
        x = Lx*(nx+mod(nx,2))+xs*(-2*mod(nx,2)+1); % x coordinate
        for ny = -(order-abs(nx)):order-abs(nx)
            y = Ly*(ny+mod(ny,2))+ys*(-2*mod(ny,2)+1); % y coordinate
            for nz = [-(order-(abs(nx)+abs(ny))),order-(abs(nx)+abs(ny))]
                z = Lz*(nz+mod(nz,2))+zs*(-2*mod(nz,2)+1); % z coordinate
                [theta, phi, r] = cart2sph(x-xr,y-yr,z-zr); % angle & distance of image-source to receiver
                if order ~= 0
                    r = r + jitter*randn; % distance with jitter (but don't jitter the direct sound)
                end % if order ~= 0
                k = round(r * fs / c)+1; % sample number
                % the following avoids jitter overshoot, and also avoids
                % infinite amplitude at r = 0
                if k>1 && k <=length(out)
                    if nchans == 1
                        out(k) = out(k) + 1/r;
                    else
                        % This is the slowest part of the loop (calling the
                        % HOA Toolbox)
                    out(k,:) = out(k,:)+ SphericalHarmonicMatrix(hoaFmt,theta,phi)'/r;
                    end
                end % if k>=1
            end % for nz =
        end % for ny =
    end % for nx =
    if o <= maxorder
        out = filter(b,a,out); % filter IR for each order
    end %if o <= maxorder
end


OUT.audio = out;
OUT.fs = fs;



OUT.chanID = makechanID(nchans,1);

OUT.funcallback.name = 'imagesource_rectangular_room_demo.m';
OUT.funcallback.inarg = {Lx,Ly,Lz,xs,ys,zs,xr,yr,zr,c,jitter,maxorder,ambiorder};



% AIR ABSORPTION FILTER
% This could be implemented by filtering the IR into a number of bands, and
% applying a decreasing gain funtion over time for each (the slope of which
% depends on the frequency of the band), and then recombining the bands.
% However, from a pragmatic perspective, it is more efficient to factor
% air absorption into the room surface absorption coefficients.

% BASIC ACOUSTICAL PARAMETERS (BROADBAND)
% Schroeder reverse integration
decay = 10*log10(flipud(cumsum(flipud(out(:,1).^2))));
G = decay(1)+20; % strength factor
decay = decay - max(decay);
Tstart = find(decay <= -5, 1, 'first'); % -5 dB
T20end = find(decay <= -25, 1, 'first'); % -25 dB
T30end = find(decay <= -35, 1, 'first'); % -35 dB
p = polyfit((Tstart:T20end)', decay(Tstart:T20end),1); %linear regression
T20 = 3*((p(2)-25)/p(1)-(p(2)-5)/p(1))/fs; % reverberation time, T20
q = polyfit((Tstart:T30end)', decay(Tstart:T30end),1); %linear regression
T30 = 2*((q(2)-35)/q(1)-(q(2)-5)/q(1))/fs; % reverberation time, T30
IRstart = find(decay < 0, 1, 'first'); % direct sound
C50 = 10*log10(1-10^(decay(IRstart+0.05*fs))/10)-decay(IRstart+0.05*fs); % clarity index
%disp(['G ',num2str(G,3),' dB    T20 ',num2str(T20,3),' s    T30 ',num2str(T30,3),' s    C50 ',num2str(C50,3),' dB'])


% visualisation and auralization
figure
subplot(3,1,1)
plot(((1:length(out(:,1)))-1)./fs,out(:,1),'r')
xlabel('Time (s)')
ylabel('Amplitude')

subplot(3,1,2)
plot(((1:length(out(:,1)))-1)./fs,decay,'g')
hold on
plot(((Tstart:T20end)-1)./fs,((Tstart:T20end)-1).*p(1)+p(2),'b')
plot(((Tstart:T30end)-1)./fs,((Tstart:T30end)-1).*q(1)+q(2),'r')
ylim([-100 0])
xlabel('Time (s)')
ylabel('Reverse-Integrated Level (dB)')
hold off

fftlen = fs*10; % 0.1 Hz resolution
magspectrum = 10*log10(abs(fft(out(:,1),fftlen)).^2);
plotcomponents = 1:floor(fftlen/2);
subplot(3,1,3)
plot((plotcomponents-1)./10,magspectrum(plotcomponents),'k')
xlim ([0 200])
xlabel('Frequency (Hz)')
ylabel('Level (dB)')

msgbox({['G     ',num2str(G,3),' dB'];['T20 ',num2str(T20,3),' s'];['T30 ',num2str(T30,3),' s'];['C50 ',num2str(C50,3),' dB']},'Result')

sound(out(:,1)./max(abs(out(:,1))),fs) % play normalized sound wave (ch1 only)


