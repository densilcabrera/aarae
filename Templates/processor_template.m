function [OUT varargout] = processor_template(IN, input_1, input_2)
% This function can be used as a template for adapting your audio
% processing functions to work in the AARAE environment.
%
% AARAE processors take the audio information stored in the AARAE tree
% display and process the input to produce an output. Unlike generator and
% calculator functions in AARAE, processor functions require an input in
% the form of a structure type variable (IN) and will output a structure
% type variable with the processed audio (OUT). The input structure will
% ALWAYS have at least the fields .audio, .fs and .datatype that you can
% use to process the audio. Processors may as well include additional
% fields to the output structure such as .bandID.
%
% You can also use these first few lines to write a brief description of
% what your function does. This will be displayed as a tooltip when the
% user hoovers the mouse over the box where your function is displayed in
% AARAE
%
% The next few lines show an example on how you may use MATLAB's built-in
% inputdlg function to allow the user to type in the input arguments your
% function requires to work.

if nargin ==1 % If the function is called within the AARAE environment it
    % will have at least one input parameter which is the audio
    % data structure that your function will process, therefore
    % you can check that the user has input all input parameters
    % if you want your function to work as standalone outside the
    % AARAE environment. You can use this input dialog to request
    % for the additional parameters that you require for your
    % function to work if they're not part of the AARAE structure
    
    param = inputdlg({'Input 1';... % These are the input box titles in the
        'Input 2'},...% inputdlg window.
        'Window title',... % This is the dialog window title.
        [1 30],... % You can define the number of rows per
        ...        % input box and the number of character
        ...        % spaces that each box can display at once
        ...        % per row.
        {'2';'3'}); % And the preset answers for your dialog.
    
    param = str2num(char(param)); % Since inputs are usually numbers it's a
    % good idea to turn strings into numbers.
    
    if length(param) < 2, param = []; end % You should check that the user
    % has input all the required
    % fields.
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
if isstruct(IN) % You should check that the function is being called within
    % the AARAE environment, if so, you can extract the
    % information you need to run your processor.
    
    % Ideally you should be able to process up to 6 dimensions of audio,
    % but if your processor cannot, then you can limit the number of
    % dimensions with the following function call:
    maxdim = 6;  % change '6' to the maximum number dimensions
    IN = choose_from_higher_dimensions(IN,maxdim,1);
    
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
    
    
    % Sometimes you may wish to get another audio input
    % The following uses an in-built AARAE function 'choose_audio' to do
    % this for you. If you wish to do the same with this function outside
    % the AARAE environment (as a stand-alone), then it might be easiest to
    % have the additional audio as an input argument to the function.
    if false % change to true if you wish to enable the following
        
        % Use a menu & dialog box to select a wav file or audio within AARAE 
        selection = choose_audio; % call AARAE's choose_audio function
        if ~isempty(selection)
            auxiliary_audio = selection.audio; % additional audio data
            fs2 = selection.fs; % sampling rate
            
            if ~(fs2 == fs)
                % match sampling rates if desired
                auxiliary_audio = resample(auxiliary_audio,fs,fs2);
            end
            [len2, chans2, bands2] = size(auxiliary_audio); % new wave dimensions
        end
    end
    
    
elseif ~isempty(param) || nargin > 1
    % If for example you want to enable your function to
    % run as a standalone MATLAB function, you can use
    % the IN input argument as an array of type double
    % and assign it to the audio your function is going
    % to process.
    audio = IN;
    fs = input_1;
end

% To make your function work as standalone you can check that the user has
% either entered at least an audio variable and it's sampling frequency.
if ~isempty(audio) && ~isempty(fs)
    % If the requirements are met, your code can be executed!
    % You may copy and paste some code you've written or write something
    % new in the lines below, such as:
    
    % It is often important to know the dimensions of your audio
    % Ideally your processor should be able to work with audio of up to 6
    % dimensions - where dimension 1 is time, dimension 2 is channels,
    % dimension 3 is bands, dimension 4 is used for multicycle test signal
    % analysis, dimension 5 is used for measurement with sequential
    % multichannel output from AARAE.
    [len,chans,bands,dim4,dim5,dim6] = size(audio);
    
    % Sometimes you might want to mixdown multiband audio
    if bands > 1
        audio = sum(audio,3);
        disp('Multiband audio has been mixed into a single band')
    end
    
    
    % Here is an example of a very simple processor operation
    audio = flipdim(audio,1);
    
    % And once you have your result, you should set it up in an output form
    % that AARAE can understand.
    
    if isstruct(IN)
        OUT = IN; % You can replicate the input structure for your output
        OUT.audio = audio; % And modify the fields you processed
        % Or simply output the fields you consider necessary after
        % processing the input audio data, AARAE will figure out what has
        % changed and complete the structure. But remember, it HAS TO BE a
        % structure if you're returning more than one field:
        %   OUT = audio;  if you just want to return the audio,
        %   or,
        %   OUT.audio = audio; if you want to return two fields.
        %   OUT.fs = fs;
        %OUT.properties.prop1 = prop1;
        %OUT.properties.prop2 = prop2; You may provide additional info
        %OUT.properties.prop3 = prop3; about your generated signal in
        %                              a structure type field called
        %                              .properties
        OUT.funcallback.name = 'processor_template.m'; % Provide AARAE
        % with the name of your function 
        OUT.funcallback.inarg = {input_1,input_2}; % assign all of the 
        % input parameters that could be used to call the function 
        % without dialog box to the output field param (as a cell
        % array) in order to allow batch processing.
    else
        % You may increase the functionality of your code by allowing the
        % output to be used as standalone and returning individual
        % arguments instead of a structure.
        OUT = audio;
    end
    varargout{1} = fs;
    % The processed audio data will be automatically displayed in AARAE's main
    % window as long as your output contains audio stored either as a single
    % variable: OUT = audio;, or it's stored in a structure along with any other
    % parameters: OUT.audio = audio;
else
    % AARAE requires that in case that the user doesn't input enough
    % arguments to generate audio to output an empty variable.
    OUT = [];
end

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