function M = HSVplotcolours1(chans, bands, Vlim)
% define plot colours in HSV colour-space
%
% Value (V) is used for for channel. The default value range is [0.4, 1],
% but this can be changed using the third input argument Vlim.
% 
% Hue (H) is used for band, with equally distributed hues around the hue 
% circle.
%
% Saturation (S) is maximum (1)

% Saturation = 1
S = ones(bands,chans);

% Value
if nargin == 3
    minV = Vlim(1);
    maxV = Vlim(2);
else
    maxV = 1;
    minV = 0.4;
end
if chans == 1
    V = minV;
else
    V = minV:(maxV-minV)/(chans-1):maxV;
end
V = repmat(V,[bands,1]);

% Hue
if bands == 1
    H = 0;
else
    H = 0:1/(bands):(1-1/bands);
end
H = repmat(H',[1,chans]);

% convert HSV to RGB
M = hsv2rgb(cat(3,H,S,V));

end % eof