function OUT = AARAE_workflow_processor(IN,filename)
% This function provides a simple method for running user-written
% workflow functions, using code adapted from AARAE's log file (or code in 
% a similar style).
%
% Workflow code should be put in AARAE's Workflows directory. Some examples
% of workflow functions are already in that directory.
%
% One purpose of using a workflow might be to avoid (or reduce instances
% of) dialog boxes, which can sometimes be an inefficient way of entering
% function parameters. Hence the parameters can be written in (and read
% from) the workflow file code instead. This can greatly speed up
% processing if you are doing the same thing repeatedly.
%
% Another purpose of using a workflow might be to concatenate a sequence of
% processors to achieve a particular goal. Again, this can reduce the
% amount of tedious interaction with AARAE's user interface if this
% workflow needs to be run many times.
%
% Note that not all dialog boxes can be avoided by using a workflow. Only
% dialog boxes that correspond to function call inputs can be avoided.
%
% Code by Densil Cabrera

if ~isdir([cd '/Workflows']), mkdir([cd '/Workflows']); end

ok = false;
while ~ok
    if nargin ==1
        %cd('/Workflows'); % default directory for workflows
        [filename, pathname] = uigetfile({'*.m'}, ...
            'Select the workflow function that you wish to run.',...
            [cd '/Workflows']);
        
        if isequal(filename,0) || isequal(pathname,0)
            % user pressed 'cancel'
            OUT = [];
            return
        end
    end
    
    % remove suffix if present
    C = strsplit(filename,'.');
    filename = C{1};
    try
        functionhandle = str2func(filename);
        if isempty(IN) && nargin(functionhandle) == 1
            IN = choose_audio;
            if isempty(IN)
                OUT = [];
                return
            end
%             h = warndlg('The selected workflow function requires audio input. Please try again.','AARAE info');
%             pause(1)
%             delete(h)
        else
            ok = true;
        end
    catch
        warndlg('Unable to interpret the selected function.','AARAE info','modal');
    end
    
end

 
try
    
    if nargin(functionhandle) == 1
        OUT = functionhandle(IN);
    else
        OUT = functionhandle();
    end
    OUT.funcallback.name = 'AARAE_workflow_processor.m';
    OUT.funcallback.inarg = {filename};

    % classify as 'results' if the output does not have audio, or does have
    % tables
    

    if isfield(OUT,'audio')
        if isempty(OUT.audio)
            % empty audio field
            OUT = rmfield(OUT,'audio');
            datatype = 4;
            OUT.datatype = datatypefield(datatype);
        end
    else 
        % no audio field
        datatype = 4;
        OUT.datatype = datatypefield(datatype);
    end
    % the following lines are commented out to preserve the audio field if
    % there are tables 26 Oct 2019
%     if isfield(OUT,'tables')
%         % tables are present
%         if isfield(OUT,'audio')
%             % remove the audio field so that the tables are available in
%             % the GUI
%             OUT = rmfield(OUT,'audio');
%         end
%         datatype = 4;
%         OUT.datatype = datatypefield(datatype);
%     end
    if ~isfield(OUT,'datatype')
        datatype = 1;
        OUT.datatype = datatypefield(datatype);
    end


catch err
    h=warndlg('AARAE workflow abandoned either because of an error in the workflow function or because the selected input data was inappropriate for the workflow function. Please refer to the Matlab error report in the Command Window.','AARAE info','modal');
    msgString = getReport(err);
    disp('AARAE workflow abandoned either because of an error in the workflow function or because the selected input data was inappropriate for the workflow function (see below).');
    disp(msgString); % displays the error message without creating an error.
    uiwait(h)
    OUT = [];
    return
end




%**************************************************************************
% Copyright (c) 2015, Densil Cabrera
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%  * Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%  * Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%  * Neither the name of the University of Sydney nor the names of its contributors
%    may be used to endorse or promote products derived from this software
%    without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%**************************************************************************