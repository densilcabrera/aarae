function [chord2,cLH,ChrdChTimes,chroTSpec] = PitchProfile(Salience,timePoints,order)
% PITCHPROFILE find the pitch profile and apply chord recognition
% Wrapper around FindChord

[r,c] = size(Salience);
Salience = [zeros(r,3) Salience];
Salience(:,120)= zeros(r,1); 

for i = 1:r
  chromata = reshape(Salience(i,:)',12,10);
  chroTSpec(i,:) = mean(chromata,2)';
  chord{i} = FindChord(chroTSpec(i,:));
end

[r,c] = size(chroTSpec);

% Calculate Change Likelihood
cLH = sum(diff(chroTSpec.^2,order).^2,2); 

% Find peaks
cLHdiff = [diff(cLH); 0];  
cLHpeaks = find(cLHdiff(1:end-1) > 0 & cLHdiff(2:end) < 0); 


if isempty(cLHpeaks) % What if there are no peaks;
    % Do nothing
elseif cLHpeaks(1) == 1 % What if the peak is first.
  cLHpeaks = cLHpeaks(2:end);
end


% Which peaks are greater than 30% of the total height?
candidates = (cLHdiff(cLHpeaks-1) ./ max(cLHdiff)) > 0.3;



% Indexes of peaks
cLHpeakIndexes = cLHpeaks(candidates)+1;
ChrdChTimes = [0; timePoints(cLHpeakIndexes); timePoints(end)];
cLHpeakIndexes = [1; cLHpeakIndexes; r];

% So we have worked out at which time the chord is changing. 
% Next is to take a mean from ChrdChTime(x) to ChrdChTime(x+1) and then
% apply the FindChord to that new frame.

try
  for i =1:length(cLHpeakIndexes)-1
    SalFrInd(i) = {[cLHpeakIndexes(i):cLHpeakIndexes(i+1)]'};
  end


  for i = 1:length(SalFrInd)
    SalienceFrame = Salience(SalFrInd{i},:); % Get the chunk between the chord change times
    [r,c] = size(SalienceFrame);
    for j = 1:r
      chromata2(j,:)= mean(reshape(SalienceFrame(j,:),12,10)');
    end
    chroTSpec2(i,:) = mean(chromata2); % Mean of the reshaped chromas
    chord2{i} = FindChord(chroTSpec2(i,:)); % Run FindChord on the chromas
  end
catch
  chord2 = {'No chord'};
  % No chords changing.
end
ChrdChTimes = ChrdChTimes(1:end-1);
end

function [chord,R2]= FindChord(chromata)
% finds the octave-spaced chord chroma profile that most closely matches the chroma profile%
% var	chroma, intervals: string;
% 	root: integer;
% 	rootX: integer;
% 	i: integer;
% 	r, SUMx, SUMx2: longreal; %correlation variables%
% 	rX: longreal;

if sum(chromata) == 0
  chord = 'none';
  return
end

% ChrdPr - one profile per row, each chroma in the 12 columns
ChrdPr= [...
2.28, 0,    0.05, 0,    0,    0.35, 0,    0,    0.19, 0,    0.03, 0;
1.44, 1.37, 0.03, 0.03, 0,    0.22, 0.21, 0,    0.13, 0.12, 0.03, 0.02;
1.43, 0,    1.53, 0,    0.03, 0.18, 0,    0.18, 0.11, 0,    0.21, 0;
1.38, 0.02, 0.03, 1.42, 0,    0.4,  0,    0,    0.66, 0,    0.02, 0.12;
1.68, 0,    0.09, 0,    1.09, 0.16, 0.02, 0,    0.09, 0.17, 0.02, 0;
0.95, 0.09, 0.02, 0.02, 0,    1.82, 0,    0.02, 0.08, 0,    0.25, 0;
1.46, 0,    0.31, 0,    0.03, 0.23, 1.51, 0,    0.29, 0,    0.02, 0.23;
0.45, 0.82, 0.01, 0.05, 0,    1.62, 0.06, 0.02, 0.04, 0.04, 0.24, 0.01;
1.48, 0.04, 0.06, 0.01, 0.48, 1.1,  0.01, 0.01, 0.1,  0.07, 0.15, 0;
0.79, 0.08, 0.72, 0.01, 0.01, 1.45, 0,    0.18, 0.06, 0,    0.46, 0;
0.77, 0.11, 0.02, 0.68, 0,    1.46, 0,    0.01, 0.3,  0,    0.17, 0.05;
0.82, 0,    1.34, 0,    0.06, 0.1,  0.95, 0.1,  0.16, 0,    0.12, 0.15;
1.37, 0,    0.26, 0,    0.81, 0.14, 0.86, 0,    0.17, 0.1,  0.02, 0.1;
1.35, 0,    0.68, 0.08, 0.01, 0.17, 0,    1.39, 0.05, 0.02, 0.1,  0;
1.08, 0.02, 0.23, 0.93, 0.02, 0.31, 1.11, 0,    0.71, 0,    0.02, 0.49;
1.34, 0.01, 0.02, 1.09, 0,    0.32, 0,    0.79, 0.32, 0.02, 0.01, 0.06;
1.62, 0,    0.05, 0.05, 0.55, 0.17, 0.01, 0.58, 0.06, 0.16, 0.01, 0;
1.22, 0.12, 0.07, 0,    1.24, 0.12, 0.07, 0,    1.26, 0.12, 0.06, 0;
0.65, 0.77, 0.04, 0.01, 1.17, 0.05, 0.21, 0,    1.12, 0.24, 0.05, 0.01;
0.33, 1.16, 0.01, 0.04, 0.06, 1.08, 0.1,  0.02, 0.94, 0.03, 0.27, 0.01;
1.11, 0.04, 0.57, 0.1,  0.01, 1.05, 0,    1.11, 0.05, 0.01, 0.31, 0;
0.96, 0.12, 0.63, 0,    1.09, 0.08, 0.05, 0.06, 1.13, 0.08, 0.18, 0;
0.67, 0.09, 1.11, 0,    0.2,  0.09, 0.73, 0.08, 1.11, 0,    0.19, 0.09;
0.51, 0.3,  0.01, 0.47, 0.05, 0.92, 0.01, 0.01, 1.28, 0,    0.19, 0.04;
0.62, 0.31, 0.61, 0.01, 0.13, 1.1,  0.01, 0.15, 1.01, 0,    0.54, 0;
0.51, 0.12, 0.09, 0.47, 0.08, 0.15, 0.49, 0,    1.38, 0,    0.04, 0.2;
0.83, 0.02, 0.61, 0.86, 0.02, 0.61, 0.88, 0.02, 0.61, 0.87, 0.02, 0.62];

interval{1} = '___';
interval{2} = '__M';
interval{3} = '__m';
interval{4} = 'm__';
interval{5} = 'M__';
interval{6} = '_P_';
interval{7} = '_d_';
interval{8} = 'M_M';
interval{9} = '_PM';
interval{10} = 'm_m';
interval{11} = '_Pm';
interval{12} = 'M_m';
interval{13} = '_dm';
interval{14} = 'SP_';
interval{15} = 'md_';
interval{16} = 'mP_';
interval{17} = 'MP_';
interval{18} = 'MA_';
interval{19} = 'mPM';
interval{20} = 'MPM';
interval{21} = 'SPm';
interval{22} = 'MAm';
interval{23} = 'Mdm';
interval{24} = 'mPm';
interval{25} = 'mdm';
interval{26} = 'MPm';
interval{27} = 'mdd';

modOffset = [0 1 2 0 0 5 6 1 5 2 5 2 6 7 0 0 0 8 1 5 7 2 2 5 2 8 6];
chromaNames = {'A ','B flat ','B ','C ','D flat ','D ','E flat ','E ','F ','F sharp ','G ','A flat '};

r       = 0;
SUMx    = 0;
SUMx2   = 0;
SUMx  = sum(chromata);
SUMx2 = sum(chromata.^2);

for i = 1:27
[rX, rootX] = ChordProfile(ChrdPr(i,:), SUMx, SUMx2, chromata);
if r < rX
  r = rX;
  if modOffset(i) == 0 
    root = rootX;
  else
    root = mod((rootX + modOffset(i)),12);
  end
  intervals = interval{i};
end
root = root + 1;
if root == 13
root = 1;
end
chroma = chromaNames{root};

chord = [chroma  intervals];
R2 = r^2;
end
end


% [rX, rootX] = ChordProfile(2.28, 0, 0.05, 0, 0, 0.35, 0, 0, 0.19, 0, 0.03, 0, SUMx, SUMx2, chromata)%0 - unison%
% if r < rX
%   r = rX;
%   root = rootX;
%   intervals = '___';
% end
% [rX, rootX] = ChordProfile(1.44, 1.37, 0.03, 0.03, 0, 0.22, 0.21, 0, 0.13, 0.12, 0.03, 0.02, SUMx, SUMx2, chromata)%01 - major 7th%
% if r < rX
%   r = rX;
%   root = mod((rootX + 1),12);
%   intervals = '__M';
% end
% [rX, rootX] = ChordProfile(1.43, 0, 1.53, 0, 0.03, 0.18, 0, 0.18, 0.11, 0, 0.21, 0, SUMx, SUMx2, chromata)%02 - minor 7th%
% if r < rX
%   r = rX;
%   root = mod((rootX + 2),12);
%   intervals = '__m';
% end
% [rX, rootX] = ChordProfile(1.38, 0.02, 0.03, 1.42, 0, 0.4, 0, 0, 0.66, 0, 0.02, 0.12, SUMx, SUMx2, chromata)%03 - minor 3rd%
% if r < rX
%   r = rX;
%   root = rootX;
%   intervals = 'm__';
% end
% [rX, rootX] = ChordProfile(1.68, 0, 0.09, 0, 1.09, 0.16, 0.02, 0, 0.09, 0.17, 0.02, 0, SUMx, SUMx2, chromata)%04 - major 3rd%
% if r < rX
%   r = rX;
%   root = rootX;
%   intervals = 'M__';
% end
% [rX, rootX] = ChordProfile(0.95, 0.09, 0.02, 0.02, 0, 1.82, 0, 0.02, 0.08, 0, 0.25, 0, SUMx, SUMx2, chromata)%05 - perfect 5th%
% if r < rX
%   r = rX;
%   root = mod((rootX + 5) , 12);
%   intervals = '_P_';
% end
% [rX, rootX] = ChordProfile(1.46, 0, 0.31, 0, 0.03, 0.23, 1.51, 0, 0.29, 0, 0.02, 0.23, SUMx, SUMx2, chromata)%06 - tritone%
% if r < rX
%   r = rX;
%   root = mod((rootX + 6) , 12);
%   intervals = '_d_';
% end
% [rX, rootX] = ChordProfile(0.45, 0.82, 0.01, 0.05, 0, 1.62, 0.06, 0.02, 0.04, 0.04, 0.24, 0.01, 		SUMx, SUMx2, chromata)%015 - incomplete major 7th%
% if r < rX
%   r = rX;
%   root = mod((rootX + 1) , 12);
%   intervals = 'M_M';
% end
% [rX, rootX] = ChordProfile(1.48, 0.04, 0.06, 0.01, 0.48, 1.1, 0.01, 0.01, 0.1, 0.07, 0.15, 0,		SUMx, SUMx2, chromata)%045 - incomplete major 7th%
% if r < rX
%   r = rX;
%   root = mod((rootX + 5) , 12);
%   intervals = '_PM';
% end
% [rX, rootX] = ChordProfile(0.79, 0.08, 0.72, 0.01, 0.01, 1.45, 0, 0.18, 0.06, 0, 0.46, 0,		SUMx, SUMx2, chromata)%025 - incomplete minor 7th%
% if r < rX
%   r = rX;
%   root = mod((rootX + 2) , 12);
%   intervals = 'm_m';
% end
% [rX, rootX] = ChordProfile(0.77, 0.11, 0.02, 0.68, 0, 1.46, 0, 0.01, 0.3, 0, 0.17, 0.05, SUMx, SUMx2, chromata)%035 - incomplete minant 7th%
% if r < rX
%   r = rX;
%   root = mod((rootX + 5) , 12);
%   intervals = '_Pm';
% end
% [rX, rootX] = ChordProfile(0.82, 0, 1.34, 0, 0.06, 0.1, 0.95, 0.1, 0.16, 0, 0.12, 0.15,		SUMx, SUMx2, chromata)%026 - incomplete minant 7th%
% if r < rX
%   r = rX;
%   root = mod((rootX + 2),12);
%   intervals = 'M_m';
% end
% [rX, rootX] = ChordProfile(1.37, 0, 0.26, 0, 0.81, 0.14, 0.86, 0, 0.17, 0.1, 0.02, 0.1,		SUMx, SUMx2, chromata)%046 - incomplete half diminished 7th%
% if r < rX
%   r = rX;
%   root = mod((rootX + 6),12);
%   intervals = '_dm';
% end
% [rX, rootX] = ChordProfile(1.35, 0, 0.68, 0.08, 0.01, 0.17, 0, 1.39, 0.05, 0.02, 0.1, 0,		SUMx, SUMx2, chromata)%027 - suspension%
% if r < rX
%   r = rX;
%   root = mod((rootX + 7),12);
%   intervals = 'SP_';
% end
% [rX, rootX] = ChordProfile(1.08, 0.02, 0.23, 0.93, 0.02, 0.31, 1.11, 0, 0.71, 0, 0.02, 0.49,		SUMx, SUMx2, chromata)%036 - diminished%
% if r < rX
%   r = rX;
%   root = rootX;
%   intervals = 'md_';
% end
% [rX, rootX] = ChordProfile(1.34, 0.01, 0.02, 1.09, 0, 0.32, 0, 0.79, 0.32, 0.02, 0.01, 0.06,		SUMx, SUMx2, chromata)%037 - minor%
% if r < rX
%   r = rX;
%   root = rootX;
%   intervals = 'mP_';
% end
% [rX, rootX] = ChordProfile(1.62, 0, 0.05, 0.05, 0.55, 0.17, 0.01, 0.58, 0.06, 0.16, 0.01, 0,		SUMx, SUMx2, chromata)%047 - major%
% if r < rX
%   r = rX;
%   root = rootX;
%   intervals = 'MP_';
% end
% [rX, rootX] = ChordProfile(1.22, 0.12, 0.07, 0, 1.24, 0.12, 0.07, 0, 1.26, 0.12, 0.06, 0,		SUMx, SUMx2, chromata)%048 - augmented%
% if r < rX
%   r = rX;
%   root = mod((rootX + 8) , 12);
%   intervals = 'MA_';
% end
% [rX, rootX] = ChordProfile(0.65, 0.77, 0.04, 0.01, 1.17, 0.05, 0.21, 0, 1.12, 0.24, 0.05, 0.01,		SUMx, SUMx2, chromata)%0148 - minor chord with major 7th%
% if r < rX
%   r = rX;
%   root = mod((rootX + 1) , 12);
%   intervals = 'mPM';
% end
% [rX, rootX] = ChordProfile(0.33, 1.16, 0.01, 0.04, 0.06, 1.08, 0.1, 0.02, 0.94, 0.03, 0.27, 0.01,		SUMx, SUMx2, chromata)%0158 - major 7th%
% if r < rX
%   r = rX;
%   root = mod((rootX + 5) , 12);
%   intervals = 'MPM';
% end
% [rX, rootX] = ChordProfile(1.11, 0.04, 0.57, 0.1, 0.01, 1.05, 0, 1.11, 0.05, 0.01, 0.31, 0,		SUMx, SUMx2, chromata)%0257 - suspended minor 7th%
% if r < rX
%   r = rX;
%   root = mod((rootX + 7) , 12);
%   intervals = 'SPm';
% end
% [rX, rootX] = ChordProfile(0.96, 0.12, 0.63, 0, 1.09, 0.08, 0.05, 0.06, 1.13, 0.08, 0.18, 0,		SUMx, SUMx2, chromata)%0248 - augmented 7th%
% if r < rX
%   r = rX;
%   root = mod((rootX + 2) , 12);
%   intervals = 'MAm';
% end
% [rX, rootX] = ChordProfile(0.67, 0.09, 1.11, 0, 0.2, 0.09, 0.73, 0.08, 1.11, 0, 0.19, 0.09,		SUMx, SUMx2, chromata)%0268 - french 6th%
% if r < rX
%   r = rX;
%   root = mod((rootX + 2) , 12);
%   intervals = 'Mdm';
% end
% [rX, rootX] = ChordProfile(0.51, 0.3, 0.01, 0.47, 0.05, 0.92, 0.01, 0.01, 1.28, 0, 0.19, 0.04,		SUMx, SUMx2, chromata)%0358 - minor 7th%
% if r < rX
%   r = rX;
%   root = mod((rootX + 5) , 12);
%   intervals = 'mPm';
% end
% [rX, rootX] = ChordProfile(0.62, 0.31, 0.61, 0.01, 0.13, 1.1, 0.01, 0.15, 1.01, 0, 0.54, 0,		SUMx, SUMx2, chromata)%0258 - half diminished 7th%
% if r < rX
%   r = rX;
%   root = mod((rootX + 2) , 12);
%   intervals = 'mdm';
% end
% [rX, rootX] = ChordProfile(0.51, 0.12, 0.09, 0.47, 0.08, 0.15, 0.49, 0, 1.38, 0, 0.04, 0.2,		SUMx, SUMx2, chromata)%0368 - minor 7th%
% if r < rX
%   r = rX;
%   root = mod((rootX + 8) , 12);
%   intervals = 'MPm';
% end
% [rX, rootX] = ChordProfile(0.83, 0.02, 0.61, 0.86, 0.02, 0.61, 0.88, 0.02, 0.61, 0.87, 0.02, 0.62,		SUMx, SUMx2, chromata)%0369 - diminished 7th%
% if r < rX
%   r = rX;
%   root = mod((rootX + 6) , 12);
%   intervals = 'mdd';
% end


function [rX, rootX] =  ChordProfile(chord, SUMx, SUMx2, chromata)
% ChordProfile is called by FindChord for correlation calculations%
% type	localarray = array (0..11)  real;
% var	i, j: integer;
% 	chordsum: longreal;
% 	chordsum2: longreal; %sum  squares%
% 	SUMxy, SUMy, SUMy2, SP, SSx, SSy: longreal; %correlation variables%
% 	chord: localarray;
rX = 0;
% chord(1) = w0;
% chord(2) = w1;
% chord(3) = w2;
% chord(4) = w3;
% chord(5) = w4;
% chord(6) = w5;
% chord(7) = w6;
% chord(8) = w7;
% chord(9) = w8;
% chord(10) = w9;
% chord(11) = w10;
% chord(12) = w11;

SUMy = sum(chord);
SUMy2 = sum(chord.^2);
SSy = SUMy2 - (SUMy^2) / 12;

for i = 1:12
  SUMxy = 0;
  for j = 1:12
    SUMxy = SUMxy + chord(j)  * chromata(mod((i + j),  12) + 1);
    SP    = SUMxy - SUMx      * SUMy / 12;
    SSx   = SUMx2 - (SUMx^2) / 12;
    if (rX < (SP / sqrt(SSx * SSy)))
      rX     = SP / sqrt(SSx * SSy);
      rootX  = i;
    end
  end % for i = 1:12  %
end % ChordProfile%
end

function [key, R2, Mm,tonic] = FindKey(chromata)
%finds the key chroma profile that most closely matches the chroma profile%

% type	localarray = array (0..11)  real;
% var	majorprofile: localarray;
% 	minorprofile: localarray;
% 	chroma: string;
% 	majsum, minsum: real;
% 	majsum2, minsum2: real; %sum  squares%
% 	maj: boolean;
% 	i, j: integer;
% 	r, SUMx, SUMxy, SUMx2, SUMy2, SP, SSx, SSy: longreal; %correlation variables%

majorprofile (1) = 11.5;
majorprofile (2) = 0.2;
majorprofile (3) = 3.7;
majorprofile (4) = 0.5;
majorprofile (5) = 4.5;
majorprofile (6) = 3.9;
majorprofile (7) = 0.1;
majorprofile (8) = 8.3;
majorprofile (9) = 0.5;
majorprofile (10) = 3.7;
majorprofile (11) = 0.7;
majorprofile (12) = 1.7;
minorprofile (1) = 9.6;
minorprofile (2) = 1.1;
minorprofile (3) = 2.8;
minorprofile (4) = 7.1;
minorprofile (5) = 0.8;
minorprofile (6) = 4.3;
minorprofile (7) = 0.0;
minorprofile (8) = 9.4;
minorprofile (9) = 5.6;
minorprofile (10) = 0.3;
minorprofile (11) = 1.3;
minorprofile (12) = 2.1;
r = 0;
majsum = 39.3;
majsum2 = 267.91;
minsum = 44.4;
minsum2 = 296.66;
SUMx = 0;
SUMx2 = 0;
for i = 1:12
  SUMx = SUMx + chromata (i);
  SUMx2 = SUMx2 + sqr(chromata (i));
end
SSy = majsum2 - sqr(majsum) / 12;
for i = 1:12
  SUMxy = 0;
  for j = 1:12
    SUMxy = SUMxy + majorprofile (j) * chromata (mod((i+j),12)+1);
    SP = SUMxy - SUMx * majsum / 12;
    SSx = SUMx2 - sqr(SUMx) / 12;
    if (r < (SP / sqrt(SSx * SSy)))
      r = SP / sqrt(SSx * SSy);
      maj = true;
      tonic = i;
    end
  end
  SSy = minsum2 - sqr(minsum) / 12;
  for i = 1:12
    SUMxy = 0;
    for j = 1:12
      SUMxy = SUMxy + minorprofile (j) * chromata (mod((i+j),12)+1);
      SP = SUMxy - SUMx * minsum / 12;
      SSx = SUMx2 - sqr(SUMx) / 12;
      if (r < (SP / sqrt(SSx * SSy)))
        r = SP / sqrt(SSx * SSy);
        maj = false;
        tonic = i;
      end
    end
    switch tonic
      case 1 
        chroma = 'A ';
      case 2 
        chroma = 'B flat ';
      case 3 
        chroma = 'B ';
      case 4 
        chroma = 'C ';
      case 5 
        chroma = 'D flat ';
      case 6 
        chroma = 'D ';
      case 7 
        chroma = 'E flat ';
      case 8 
        chroma = 'E ';
      case 9 
        chroma = 'F ';
      case 10 
        chroma = 'F sharp ';
      case 11 
        chroma = 'G ';
      case 12 
        chroma = 'A flat ';
    end
    if maj
      key = chroma + ' major';
      Mm = 0;
    else
      key = chroma + ' minor';
      Mm = -3;
    end
    R2 = r^2;
  end %FindKey%
end
end


function TallyChromaAndPitch
%TallyChromaAndPitch chroma and pitch patterns
for i = 1:12
  ChromaSalience(i) = 0;
end
for i = 0:87
  PitchSalience(i) = 0;
end
for i= 1:PitchCount
  FindChromaName(CompoundPitch(i));
  ChromaSalience (round(Pitch)) = ChromaSalience (round(Pitch)) + Salience(i) * (100 - Cents)/100;
  if (Note <= 87) and (Note >= 0)
    PitchSalience(Note) = PitchSalience(Note) + Salience(i) * (100 - Cents)/100;
  end
  if Cents ~= 0 
    if strcmp(PlusMinus,'-')
      ChromaSalience(mod(round(Pitch) + 11),12) = 	ChromaSalience(mod(round(Pitch) + 11), 12) + Salience(i) * Cents/100;
      if Note > 0
        PitchSalience(Note - 1) = PitchSalience(Note - 1)	+ Salience(i) * Cents/100;
      end
    elseif strcmp(PlusMinus,'+')
      ChromaSalience(mod(round(Pitch) + 1), 12) =		ChromaSalience(mod(round(Pitch) + 1), 12)+ Salience(i) * Cents/100;
      if Note < 87
        PitchSalience(Note + 1) = PitchSalience(Note + 1)	+ Salience(i) * Cents/100;
      end
    end

  end 
end 
end 

function Chroma = FindChromaName(f)
% FINDCHROMANAME computes the chroma and detuning of the frequency 'f'
Pitch = 12 * ((log2(f) - log2(ATune) + 50) - trunc(log2(f) - log2(ATune) + 50));
Cents = round(100 * (Pitch - round(Pitch)));
PlusMinus = '+';
if Cents < 0
  PlusMinus = '-';
end
Cents = abs(Cents);
Note = round(12 * log2(f/ATune) + 48);
if round(Pitch) == 12
  Pitch = 0;
end
switch round(Pitch)
  case {0, 12}
    Chroma = 'A ';
  case 1
    Chroma = 'A#';
  case 2
    Chroma = 'B ';
  case 3
    Chroma = 'C ';
  case 4
    Chroma = 'C#';
  case 5
    Chroma = 'D ';
  case 6
    Chroma = 'D#';
  case 7
    Chroma = 'E ';
  case 8
    Chroma = 'F ';
  case 9
    Chroma = 'F#';
  case 10
    Chroma = 'G ';
  case 11 
    Chroma = 'G#';
end %case round(Pitch) of}
end %procedure FindChromaName (f: longreal);}

    
    

function y = sqr(x)
y = x^2;
end