function [OUT,varargout] = squarewave(dur, fs, f, duty)
% Square-wave generator

    if nargin == 0 % If the function is called within the AARAE environment it
        % won't have any input arguments, this is when the inputdlg
        % function becomes useful.

        param = inputdlg({'Duration [s]';... % These are the input box titles in the
            'Sampling frequency [samples/s]';...
            'Fundamental frequency [Hz]';...
            'Duty cycle [0-100]'},...% inputdlg window.
            'Window title',... % This is the dialog window title.
            [1 60],... % You can define the number of rows per
            ...        % input box and the number of character
            ...        % spaces that each box can display at once
            ...        % per row.
            {'1';'48000';'20';'50'}); % And the preset answers for your dialog.

        param = str2num(char(param)); % Since inputs are usually numbers it's a
        % good idea to turn strings into numbers.

        if length(param) < 4, param = []; end % You should check that the user
        % has input all the required
        % fields.
        if ~isempty(param) % If they have, you can then assign the dialog's
            % inputs to your function's input parameters.
            dur = param(1);
            fs = param(2);
            f = param(3);
            duty = param(4);
        end
    else
        param = [];
    end
    
    % To make your function work as standalone you can check that the user has
    % either entered some parameters as inputs or that the inputs have been
    % acquired through the input dialog window.
    if ~isempty(param) || nargin ~= 0
        % If there are some input parameters to work with then your code can be
        % executed! You may copy and paste some code you've written or write
        % something new in the lines below, such as:
        
        t = linspace(0,dur,fs*dur);
        audio = square(2*pi*f.*t',duty);
        
        % And once you have your result, you should set it up in an output form
        % that AARAE can understand.
        
        if nargin == 0
            OUT.audio = audio; % You NEED to provide the audio you generated.
            OUT.fs = fs;       % You NEED to provide the sampling frequency of your audio.
            OUT.tag = ['Square' num2str(dur)];      % You may assign it a name to be identified in AARAE.
            OUT.funcallback.name = 'squarewave.m';
            OUT.funcallback.inarg = {dur, fs, f, duty};
        end
        
        % You may choose to increase the functionality of your code by allowing
        % it to operate outside the AARAE environment you may want to output
        % independent variables instead of a structure...
        if nargin ~= 0
            OUT = audio;
            varargout{1} = fs;
            %varargout{2} = ?;
        end
    else
        % AARAE requires that in case that the user doesn't input enough
        % arguments to generate audio to output an empty variable.
        OUT = [];
    end
    
end % End of function

%**************************************************************************
% Copyright (c) 2014, Daniel Jimenez & Densil Cabrera
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
%  * Neither the name of the <ORGANISATION> nor the names of its contributors
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