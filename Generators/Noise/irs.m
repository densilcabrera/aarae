function  OUT = irs(n,cycles,fs)
% This function is used to generate an inverse repeated repeated sequence
% (IRS), made from a maximum length sequence signal (mls), which can be
% used to measure impulse responses.
%
% The IRS signal consists of an MLS signal repeated once, with all even
% samples multiplied by -1. The n input argument (or first field of the
% dialog box) determines the bit length of the MLS signal (note that the
% IRS signal is double the length of the MLS signal used). The number of
% samples in one cycle of the IRS signal is 2*(2^n-1).
%
%
% This function calls code by M.R.P. Thomas  - please refer to
% the following folder AARAE/Generators/Noise/mls, which contains his code,
% documentation and license.
%
% The function outputs the IRS sequence (in the audio field), together with
% the time-reversed IRS signal (as audio2). The signal is time reversed for
% compatability with the general use of audio2 as an inverse filter.
% However, you should not normally use AARAE's '*' button (which convolves
% audio with audio2) to obtain the impulse response, although it will
% probably still work to an extent, because it does linear convolution
% rather than the required circular convolution (or cross-correlation with
% the non-reversed signal). Instead use the processor CircXcorrforIR to
% derive the impulse response, which is in AARAE's Cross & auto functions
% folder (in Processors).
%
% code by Densil Cabrera
% Version 1 (1 August 2014)


if nargin == 0
    param = inputdlg({'Bit length of MLS [2-24]';...
                       'Number of cycles [2 or more]';...
                       'Sampling frequency [Hz]'},...
                       'IRS input parameters',1,{'16';'2';'48000'});
    param = str2num(char(param));
    if length(param) < 3, param = []; end
    if ~isempty(param)
        n = round(param(1));
        cycles = round(param(2));
        if cycles < 2, cycles = 2; end
        fs = param(3);
    end
else
    param = [];
end
if ~isempty(param) || nargin ~= 0
    
    if n >=2 && n <=24
        %inefficient, but avoids any need to edit Thomas' code
        [~, mls] = GenerateMLSSequence(2, round(n), 0);
        irs = [mls'; mls'];
        irs(2:2:end) = -irs(2:2:end);
        
    else
        OUT = [];
        return
    end

    OUT.audio = [repmat(irs,[cycles,1]);zeros(size(irs))];
    OUT.audio2 = flipud(irs);
    OUT.fs = fs;
    OUT.tag = ['IRS' num2str(n)];
    OUT.properties.n = n;
    OUT.properties.combinehalves = 1; % used to set dialog box default in CircXcorrforIR
    OUT.properties.cycles = cycles;
    OUT.funcallback.name = 'irs.m';
    OUT.funcallback.inarg = {n,cycles,fs};
else
    OUT = [];
end
end % End of function