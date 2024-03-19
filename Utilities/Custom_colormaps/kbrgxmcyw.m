function cmap = kbrgxmcyw(n,test)
% This function returns an n-value colormap, named 'kbrgxmcyw'. The
% colormap has the following attributes:
%
% * lightness increases almost monotonically from black to white as per Lab
%   color scheme (which is a perception model).
% * a large number of hues is distingushable, potentially allowing more
%   detailed reading of the data than typical colormaps.
% * the center of the colormap is mid-grey, allowing the colormap to
%   represent bipolar data.
%
% Data for the colormap was generated from:
%   https://davidjohnstone.net/lch-lab-colour-gradient-picker
% To generate it, the background color was specified as white, and 250
% colors were interpolated in Lab space from 26 stops.
% Settings can be recreated via:
% https://davidjohnstone.net/lch-lab-colour-gradient-picker#000000,001129,00026d,320066,560057,690038,7f0000,693800,515100,456225,307038,347a56,777777,9e7876,bd7696,bd86bc,b797d1,aaa8e3,94bbf4,77d3d2,87e1ae,a0f099,d0f4a3,fef881,fff3b9,ffffff
%
% The process of generating the colormap was to first set the HSV hue at 30
% degree intervals (but with grey dividing the upper and lower halves), and
% then set the Lab lightness at approximately even increments over the
% whole colormap range (note that the lightness adjustment shifts the hues
% a little). The name of the colormap refers to the order of hues behind
% this approach. The rationale for the order is that the colours black,
% blue, red, green, magenta, cyan, yellow, and white are inherently ordered
% from dark to light at their maximum HSV values in an RGB cube.
%
% The optional first input argument (n) determines the number of rows
% (colours) in the colormap. It is an integer in the range 1 to 10000 (but
% see below for negative integers). An odd-number of n has mid-grey in the
% centre, which can be useful when n is small. If the argument is omitted
% then a 250-value colormap is returned. To flip the colormap use a
% negative value of n.
%
% To quickly test the colormap using Matlab's peaks function, use the
% second input argument (true).
%
% Example function calls:
%  kbrgxmcyw; % returns a 250-value colormap
%  kbrgxmcyw(17); % returns a 17-value colormap
%  kbrgxmcyw(17,1); % returns a 17-value colormap and demostrates it
%  kbrgxmcyw(-128,1); % returns a flipped 128-value colormap and demostrates it
%
%
% Densil Cabrera 2023


cmaphex = {'#000000', '#000206', '#00040c', '#000511', '#000715', '#000919', '#000b1c', '#000d1f', '#000e22', '#001026', '#001129', '#001030', '#001036', '#000f3d', '#000e43', '#000e4a', '#000c51', '#000a58', '#00085f', '#000566', '#00026d', '#00016c', '#00016c', '#02006b', '#0f006a', '#180069', '#1f0069', '#250068', '#2a0067', '#2e0067', '#330066', '#370064', '#3b0063', '#3f0061', '#430060', '#47005e', '#4a005d', '#4d005b', '#50005a', '#540058', '#560056', '#590053', '#5b0050', '#5d004d', '#5f004a', '#610047', '#630044', '#650041', '#66003e', '#68003a', '#6a0037', '#6c0032', '#6f002d', '#710028', '#730023', '#76001e', '#780019', '#7a0013', '#7c000c', '#7e0004', '#7f0000', '#7c0c00', '#7a1500', '#781c00', '#762200', '#742700', '#722b00', '#6f2f00', '#6d3200', '#6b3600', '#683900', '#663c00', '#643e00', '#624100', '#604400', '#5d4600', '#5b4900', '#584b00', '#564d00', '#535000', '#515200', '#505300', '#4f5500', '#4e5701', '#4d5807', '#4c5a0e', '#4b5c14', '#495d19', '#485f1e', '#466122', '#446326', '#436428', '#41652a', '#3f672b', '#3e682d', '#3b6a2f', '#396b31', '#376c33', '#346e35', '#326f37', '#307039', '#31713c', '#31723f', '#327342', '#327445', '#337548', '#33764b', '#33774e', '#347851', '#347954', '#387a57', '#407a5b', '#487a5e', '#4f7a61', '#567965', '#5c7968', '#63796b', '#69786f', '#6e7872', '#747775', '#797777', '#7d7777', '#817777', '#857877', '#897877', '#8d7876', '#917876', '#957876', '#987876', '#9c7876', '#a07878', '#a3787b', '#a6787e', '#a97881', '#ac7784', '#af7788', '#b2778b', '#b5778e', '#b97691', '#bc7695', '#bd7798', '#bd789c', '#bd7aa0', '#bd7ca3', '#bd7da7', '#bd7fab', '#bd81af', '#bd82b3', '#bd84b7', '#bd85ba', '#bd87bd', '#bc89bf', '#bc8ac1', '#bb8cc4', '#bb8ec6', '#ba90c8', '#b991ca', '#b993cc', '#b895ce', '#b796d0', '#b698d2', '#b59ad4', '#b49cd6', '#b39dd8', '#b19fd9', '#b0a1db', '#afa2dd', '#ada4df', '#aca6e1', '#aaa7e2', '#a9a9e4', '#a7abe6', '#a5ade8', '#a3afe9', '#a1b1eb', '#9eb3ed', '#9cb5ee', '#9ab7f0', '#97b9f2', '#95baf4', '#92bdf2', '#90bfee', '#8dc2eb', '#8bc4e7', '#88c7e4', '#85c9e1', '#82cbdd', '#7fceda', '#7bd0d6', '#78d2d3', '#79d4cf', '#7bd5cc', '#7dd7c8', '#7ed8c5', '#80dac1', '#81dbbd', '#83dcba', '#84deb6', '#86dfb2', '#87e1af', '#89e2ac', '#8ce4aa', '#8fe5a8', '#91e7a6', '#94e8a4', '#96eaa2', '#99eba0', '#9bed9e', '#9dee9c', '#a0f099', '#a4f09a', '#a9f19b', '#aef19c', '#b3f29d', '#b8f29e', '#bdf29f', '#c2f3a0', '#c6f3a1', '#cbf4a2', '#cff4a3', '#d4f4a0', '#d9f59d', '#def59a', '#e3f696', '#e8f693', '#ecf68f', '#f1f78c', '#f5f788', '#f9f885', '#fef881', '#fff886', '#fff78c', '#fff792', '#fff698', '#fff69d', '#fff5a3', '#fff5a8', '#fff4ae', '#fff4b3', '#fff3b9', '#fff4c0', '#fff5c7', '#fff7ce', '#fff8d5', '#fff9dc', '#fffae3', '#fffbea', '#fffdf1', '#fffef8', '#ffffff'};
[R,G,B] = deal(zeros(250,1));
for i = 1:250
    R(i) = hex2dec(cmaphex{i}(2:3))./255;
    G(i) = hex2dec(cmaphex{i}(4:5))./255;
    B(i) = hex2dec(cmaphex{i}(6:7))./255;
end
cmap = [R,G,B];
if nargin < 2
    test = false;
end
if nargin == 0
    return
end
if n == 250 || abs(round(n)) == 0
    if test
        testsurf(cmap)
    end
    return
end
if n<0
    doflip = true;
else
    doflip = false;
end
n = abs(round(n));
if n > 10000
    n = 10000; 
    disp('kbrgxmcyw colormap has a maximum of 10000 values (tweak the code if you need more)')
end
% upsample by interpolation from 250 to 10000 values
cmap = [interp(R,40),interp(G,40),interp(B,40)];
cmap(cmap>1) = 1;
cmap(cmap<0) = 0;
len = length(cmap);
% decimate to n values, maintaining first and last values 
% (except when n is 1)
if n > 2
    cmap = cmap(round([1:(len/(n-1)):((n-1.5)*len/(n-1)),len]),:);
elseif n == 2
    cmap = cmap([1 len],:);
elseif n == 1
    cmap = cmap(round(len/2),:);
end
if doflip, cmap = flip(cmap); end
if test, testsurf(cmap); end
end

function testsurf(cmap)
figure1 = figure('Name',['Colormap test, n = ', num2str(size(cmap,1))], 'Color','w');
axes1 = axes('Parent',figure1);
surf(peaks(80),'Parent',axes1,'EdgeAlpha',0.25);
colormap(cmap)
set(axes1,'CLim',[-8 8]);
colorbar(axes1);
end