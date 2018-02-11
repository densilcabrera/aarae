function [OUT, varargout] = generator_template(input_1, input_2)
% This function can be used as a template for adapting your audio
% generating functions to work in the AARAE environment.
%
% AARAE generators require that the output of this function is given in the
% form of a structure type variable. You may use 'OUT' as the name of your
% output structure. You can design your function to take as many input
% arguments as you require, but usually, since AARAE doesn't allow to input
% parameters as if you were executing the function from MATLAB's Command
% Window, these parameters are requested upon the functions' call in AARAE
% through an input dialog window (See MATLAB help inputdlg). Nonetheless,
% it's useful to declare the input arguments if you'd like your function to
% be useful as a standalone function in MATLAB.
%
% You can also use these first few lines to write a brief description of
% what your function does. This will be displayed as a tooltip when the
% user hoovers the mouse over the box where your function is displayed in
% AARAE
%
% The next few lines show an example on how you may use MATLAB's built-in
% inputdlg function to allow the user to type in the input arguments your
% function requires to work.

if nargin == 0 % If the function is called within the AARAE environment it
    % won't have any input arguments, this is when the inputdlg
    % function becomes useful.
    
    param = inputdlg({'Input 1';... % These are the input box titles in the
        'Input 2'},...% inputdlg window.
        'Window title',... % This is the dialog window title.
        [1 60],... % You can define the number of rows per
        ...        % input box and the number of character
        ...        % spaces that each box can display at once
        ...        % per row.
        {'1';'48000'}); % And the default answers for your dialog.
    
    param = str2num(char(param)); % Since inputs are usually numbers it's a
    % good idea to turn strings into numbers.
    
    if length(param) < 2, param = []; end % You should check that the user
    % has input all the required fields.
    if ~isempty(param) % If they have, you can then assign the dialog's
        % inputs to your function's input parameters.
        input_1 = param(1);
        input_2 = param(2);
    else
        % get out of here if the user presses 'cancel'
        OUT = [];
        return
    end
else
    param = [];
end


% Normally generators do not have audio inputs. If you wish to make a
% generator that has audio input, first think about whether it would be
% better classified as a processor. If it is really best as a generator,
% then audio can be input using AARAE's choose_audio function.
% If you wish to do the same with this function outside
% the AARAE environment (as a stand-alone), then it might be easiest to
% have the input audio as an input argument to the function.
if false % change to true if you wish to enable the following
    
    % Use a menu & dialog box to select a wav file or audio within AARAE
    selection = choose_audio; % call AARAE's choose_audio function
    if ~isempty(selection)
        audio = selection.audio; % audio data
        fs = selection.fs; % sampling rate
        [len, chans, bands] = size(audio); % input audio dimensions
    end
end
    
    % To make your function work as standalone you can check that the user has
    % either entered some parameters as inputs or that the inputs have been
    % acquired through the input dialog window.
    if ~isempty(param) || nargin ~= 0
        % If there are some input parameters to work with then your code can be
        % executed! You may copy and paste some code you've written or write
        % something new in the lines below, such as:
        
        t = linspace(0,input_1,input_2*input_1);
        audio = cos(2*pi*1000.*t');
        fs = input_2;
        
        % And once you have your result, you should set it up in an output form
        % that AARAE can understand.
        
        
        OUT.audio = audio; % You NEED to provide the audio you generated.
        %OUT.audio2 = ?;     You may provide additional audio derived from your function.
        OUT.fs = fs;       % You NEED to provide the sampling frequency of your audio.
        %OUT.tag = tag;      You may assign it a name to be identified in AARAE.
        %OUT.properties.prop1 = prop1;
        %OUT.properties.prop2 = prop2; You may provide additional info
        %OUT.properties.prop3 = prop3; about your generated signal in
        %                              a structure type field called
        %                              .properties
        OUT.funcallback.name = 'generator_template.m'; % Provide AARAE
        % with the name of your function
        OUT.funcallback.inarg = {input_1,input_2}; % assign all of the
        % input parameters that could be used to call the function
        % without dialog box to the output field param (as a cell
        % array).

        

    else
        % AARAE requires that in case that the user doesn't input enough
        % arguments to generate audio to output an empty variable.
        OUT = [];
    end
    
end % End of function

%**************************************************************************
% Copyright (c) <YEAR>, <OWNER>
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