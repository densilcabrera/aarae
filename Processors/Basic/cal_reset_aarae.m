function out = cal_reset_aarae(in,outcal,incal)
% This function resets the .cal field of an AARAE audio object to a
% user-specified value. Gain is applied to the underlying audio waveform to
% compensate for the .cal field change. For multichannel audio, this
% function will reset the .cal to be the same for all channels.
%
% This function can also be used to interpret the cal value of audio so as
% to apply gain (hence it could be called by another AARAE function to
% simplify the code used to interpret the cal value). The third input
% argument (incal) is used to specify the original cal value for cases
% where the input is not in the AARAE structure format (it is not used if 
% the input is an AARAE structure).
%
% Code by Densil Cabrera
% version 1.01 (2 June 2014)

if isstruct(in)
    %audio = in.audio;
    if isfield(in,'cal')
        cal = in.cal;
        [len,chans,bands,dim4,dim5,dim6] = size(in.audio);
        if length(cal) ~= chans
            if length(cal) == 1
                cal = repmat(cal,[1,chans]);
            else
                warndlg('Cal size does not match audio channels','AARAE info','modal');
                out = [];
                return
            end
        end
    else
        h = warndlg('Cal field does not exist','AARAE info','modal');
        uiwait(h)
        out = [];
        return
    end
    
else
    [len,chans,bands,dim4,dim5,dim6] = size(in);
    if ~exist ('incal','var')
        cal = zeros(1,chans); % default 0 dB cal for input
    else
        cal = incal;
        if length(incal) ~= chans
            if length(incal) == 1
                cal = repmat(incal,[1,chans]);
            else
                h=warndlg('incal size does not match audio channels','AARAE info','modal');
                uiwait(h)
                out = [];
                return
            end
        end
    end
end




if nargin == 1
    prompt = {'Desired cal value (dB)'};
    dlg_title = 'Reset calibration';
    num_lines = 1;
    def = {'0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if ~isempty(answer)
        outcal = str2num(char(answer{1}));
    else
        out = [];
        return
    end
end


if ~isempty(outcal)
    if isstruct(in)
        out = in;
        gain = 10.^((cal - outcal)./20);
        out.cal = repmat(outcal,[1,chans]);
        out.audio = out.audio .* repmat(gain,[len,1,bands,dim4,dim5,dim6]);
        % to do: consider how to operate on audio2
        out.funcallback.name = 'cal_reset_aarae.m';
        out.funcallback.inarg = {outcal};
    else
        gain = 10.^((cal - outcal)./20);
        out = in .* repmat(gain,[len,1,bands,dim4,dim5,dim6]);
    end
else
    out = [];
end

end