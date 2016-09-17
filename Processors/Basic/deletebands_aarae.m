function out = deletebands_aarae(in,S)
% This function allows you to select the bands that you wish to retain,
% and delete the remaining bands

if isfield(in,'bandID')
    param = in.bandID;
else
    param = [];
end

if ~isempty(param) 
    if nargin < 2
        [S,ok] = listdlg('Name','Band selection',...
            'PromptString','Delete unselected bands',...
            'ListString',num2str(param'),...
            'ListSize', [160 320]);
    else
        ok = 1;
    end

    if ok == 1 && ~isempty(S)
        try
            out.audio = in.audio(:,:,S,:,:,:);
            if isfield(in,'bandID')
                out.bandID = in.bandID(S);
            else
                out.bandID = S;
            end
            out.funcallback.name = 'deletebands_aarae.m';
            out.funcallback.inarg = {S};
        catch
            out = [];
        end
    else
        out = [];
    end
else
    out = [];
    h = warndlg('No bands available','AARAE info','modal');
    uiwait(h)
end