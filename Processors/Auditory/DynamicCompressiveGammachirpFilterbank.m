function out = DynamicCompressiveGammachirpFilterbank(in,fs,compressive,NumCh,Hicutoff,Locutoff,OutMidCrct)
% This function interfaces AARAE with Toshio Irino's GCFBv208 filters
%
% Note that this uses mex files, which may not be compatible with all
% systems.
%
% Please read the license for this toolbox:
%
% Toshio Irino, Wakayama university (National university corporation, Japan)
% Copyright (c) 2006, 2007
%
% Permission to use, copy, modify, and distribute this software without 
% fee is hereby granted FOR RESEARCH/EDUCATION PURPOSES only, provided 
% that this copyright notice appears in all copies and in all supporting 
% documentation, and that the software is not redistributed for any 
% fee (except for a nominal shipping charge). 
%
% For any other uses of this software, in original or modified form, 
% including but not limited to consulting, production or distribution
% in whole or in part, specific prior permission must be obtained 
% from the author.
% Signal processing methods and algorithms implemented by this
% software may be claimed by patents owned by ATR or others.
%
% The author makes no representation about the suitability of this 
% software for any purpose.  It is provided "as is" without warranty 
% of any kind, either expressed or implied.  
% Beware of the bugs.

% Interface file by Densil Cabrera 
% Version 0 (5 November 2013)
if exist('OutMidCrct','var'), GCparam.OutMidCrct = OutMidCrct; end
if nargin < 7, GCparam.OutMidCrct = 1; end
if nargin < 6, Locutoff = 100; end
if nargin < 5, Hicutoff = 6000; end
if exist('NumCh','var'), GCparam.NumCh = NumCh; end
if nargin < 4, GCparam.NumCh = 75; end
if nargin < 3
% Dialog box for settings 
% Prompt for 3 named parameters:
    prompt = {'Passive or Compressive GC filter [0 | 1]', ...
              'Number of bands', ...
              'Highest band frequency (Hz)', ...
              'Lowest band frequency (Hz)', ...
              'Outer and middle ear correction [0 | 1]'};
% Title of the dialog box
    dlg_title = 'Settings'; 
    num_lines = 1;
% Default values
    def = {'1','75','6000','100','1'}; 
 
    answer = inputdlg(prompt,dlg_title,num_lines,def);
     
 % Set function variables from dialog box user-input values
     if ~isempty(answer)
         compressive = str2num(answer{1,1});
         GCparam.NumCh = str2num(answer{2,1});
         Hicutoff = str2num(answer{3,1});
         Locutoff = str2num(answer{4,1});
         GCparam.OutMidCrct = str2num(answer{5,1});
     else
         out = [];
         return
     end
end
if isstruct(in)
% get audio signal from AARAE field
in = choose_from_higher_dimensions(in,2,1); 
    audio = in.audio;
    GCparam.fs = in.fs;
elseif ~isempty(answer)
    SndIn = in;
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        GCparam.fs = str2num(char(fs));
    else
        GCparam.fs = fs;
    end
end
if ~isempty(GCparam.fs) && ~isempty(compressive) && ~isempty(GCparam.NumCh) && ~isempty(Hicutoff) && ~isempty(Locutoff) && ~isempty(GCparam.OutMidCrct)
    % set gammachirp parameters
    GCparam.FRange = [Locutoff Hicutoff];

    GCparam = GCFBv208_SetParam(GCparam);
    
    [cGCout,pGCout]=deal(zeros(GCparam.NumCh,length(audio),size(audio,2)));
    for ch = 1:size(audio,2)
        SndIn = audio(:,ch,1,1,1,1)';
%     SndIn = mean(SndIn,3); % mixdown the third dimension, if it exists - not needed anymore
%     SndIn = SndIn(:,1); % select the first channel in the case of multichannel
%     SndIn = SndIn'; % transpose (1 row of audio signal is required for GCFBv208)

    

    [cGCout(:,:,ch), pGCout(:,:,ch), GCparam, GCresp] = GCFBv208(SndIn,GCparam);
    
    
    end
    
    
    if ~isstruct(in)
        if compressive
            out = permute(cGCout,[2,3,1]);
        else
            out = permute(pGCout,[2,3,1]);
        end
    else
        if compressive
            out.audio = permute(cGCout,[2,3,1]);
        else
            out.audio = permute(pGCout,[2,3,1]);
        end
        out.funcallback.name = 'DynamicCompressiveGammachirpFilterbank.m';
        out.funcallback.inarg = {GCparam.fs,compressive,GCparam.NumCh,Hicutoff,Locutoff,GCparam.OutMidCrct};
    end
end