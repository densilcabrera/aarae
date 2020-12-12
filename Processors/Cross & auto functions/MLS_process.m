function OUT = MLS_process(IN, offset, DCCoupling, d2stack, cycles, n)
% This function is used to analyse signals that were recorded using the MLS
% generator within AARAE. The output is an impulse response.
%
% This function calls code by M.R.P. Thomas  - please refer to
% the following folder AARAE/Generators/Noise/mls, which contains his code,
% documentation and license.
%
% SETTINGS
%   Offset: Positive offset causes negative time shift (default = 0).
%
%   DC recovery: 1 to recover DC value, 0 to not do so (e.g. for
%   loudspeaker measurements. See Thomas' documentation for more
%   information.
%
%   IR stack in dimension 2 (if available): in AARAE, dimenson 2 is used
%   for channels, and if it is singleton, then multiple IRs can be stacked
%   in dimension 2 instead of in dimension 4. If d2stack == 1, then this
%   will be done; if d2stack == -1 then no stacking is done, otherwise IRs
%   are stacked in dimension 4.
%   
% code by Densil Cabrera
% version 1 (31 July 2014)

if isstruct(IN)
    audio = IN.audio;
    %fs = IN.fs;
    %     if isfield(IN,'audio2')
    %         audio2 = IN.audio2;
    %     end
    if isfield(IN,'properties')
        if isfield(IN.properties,'cycles') && isfield(IN.properties,'n')
            cycles = IN.properties.cycles;
            n = IN.properties.n;
        else
            disp('Required properties fields not found')
        end
    else
        disp('Required properties fields not found')
    end
    OUT = IN;
else
    audio = IN;
end


% Dialog box parameters
if nargin ==1
    param = inputdlg({'Offset in samples';... %
        'DC recovery [0 | 1]';
        'IR stack in dimension 2 if available [0 | 1]'},...% use -1 for no stacking
        'MLS process settings',...
        [1 60],...
        {'0';'1';'1'});
    param = str2num(char(param));
    if length(param) < 3, param = []; end
    if ~isempty(param)
        offset = param(1);
        DCCoupling = param(2);
        d2stack = param(3);
    else
        OUT=[];
        return
    end
else
    param = [];
end

[~, chans, bands, dim4, dim5, dim6] = size(audio);


% average phase-complementary pairs of signals (if they exist)
if isfield(IN,'properties')
    if isfield(IN.properties,'complementarysignals')
        if IN.properties.complementarysignals && isfield(IN.properties,'complementarysignalsoffset')
            if isfield(IN.properties,'startflag')
                startflag = IN.properties.startflag;
            else
                startflag = 1;
            end
            for d=1:length(startflag)
                startindex1 = startflag(d);
                endindex1 = startindex1 + IN.properties.complementarysignalsoffset-1;
                startindex2 = startflag(d) + IN.properties.complementarysignalsoffset;
                endindex2 = startindex2 + IN.properties.complementarysignalsoffset-1;
                audio(startindex1:endindex1,:,:,:,:,:) = ...
                mean(cat(7,audio(startindex1:endindex1,:,:,:,:,:),...
                -audio(startindex2:endindex2,:,:,:,:,:)),7);
            end
        end
    end
end

% Stack IRs in dimension 4 (or 2) if AARAE's multi-cycle mode was used
if d2stack ~= -1
    if isfield(IN,'properties')
        if isfield(IN.properties,'startflag') && dim4==1
            startflag = IN.properties.startflag;
            dim4 = length(startflag);
            audiotemp = zeros((cycles+1)*(2^n-1),chans,bands,dim4,dim5,dim6);
            for d=1:dim4
                audiotemp(:,:,:,d,:,:) = ...
                    audio(startflag(d):startflag(d)+(cycles+1)*(2^n-1)-1,:,:,1,:,:);
            end
        end
    end
    
    
    if exist('audiotemp','var')
        audio = audiotemp;
    end
    
    if d2stack == 1 && chans == 1 && dim4 > 1
        audio = permute(audio,[1,4,3,2,5,6]);
        chans = dim4;
        dim4 = 1;
    end
end  
    
impalign = 0; % not used by generator
ir = zeros(2^n,chans,bands,dim4,dim5,dim6);
for d6 = 1:dim6
    for d5 = 1:dim5
        for d4=1:dim4
            for b = 1:bands
                ir(:,:,b,d4,d5,d6) = AnalyseMLSSequence(audio(:,:,b,d4,d5,d6),offset,cycles,n,DCCoupling,impalign);
            end
        end
    end
end


    
    


if isstruct(IN)
    OUT = IN;
    %OUT = rmfield(OUT,'audio2');
    OUT.audio = ir;
    OUT.funcallback.name = 'MLS_process.m';
    OUT.funcallback.inarg = {offset, DCCoupling, d2stack, cycles, n};
else
    OUT = ir;
end