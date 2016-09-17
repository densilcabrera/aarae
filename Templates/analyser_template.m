function OUT = directivitydemo_ASM2016(IN,Bandfilter,Autocrop,Starttime,Endtime)
% This function 
%


if nargin ==1 

    param = inputdlg({'Broadband (0), Octave band (1) or 1/3-octave band (3)?';... % These are the input box titles in the
                      'Autocrop start: No [0]; As ensemble [1]; Each channel individually [2]?';...
                      'Start time (s)';...
                      'End time (s)'},...
                      'Directivity demo',... 
                      [1 60],... 
                      {'1';'1';'0';'0.01'}); % Default values

    param = str2num(char(param)); 

    if length(param) < 4, param = []; end 
    if ~isempty(param) 
        Bandfilter = param(1);
        Autocrop = param(2);
        Starttime = param(3);
        Endtime = param(4);
    else
        % get out of here if the user presses 'cancel'
        OUT = [];
        return
    end
else
    param = [];
end


% *************************************************************************
if isstruct(IN) 
                
    
    IN = choose_from_higher_dimensions(IN,3,1); 

    
    
    if isfield(IN,'cal')
        IN = cal_reset_aarae(IN,0); % apply gain to set cal offset to 0 dB
    end
    
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
    
    
    % chanID is a cell array of strings describing each channel
    if isfield(IN,'chanID') % Get the channel ID if it exists
        chanID = IN.chanID;
    else
        % or make a chanID using AARAE's utility function
        chanID = makechanID(size(audio,2),0);
    end
    % ChanID has several formats (format 0 is used in the above call of
    % makechanID). Format 0 is simply a listing of channel numbers. Other
    % formats identify spherical harmonics (order and degree) or provide
    % spatial coordinates (Cartesian or polar) corresponding each channel.
    % Refer to the comment in makechanID for specific information about
    % this.
    % If you wish to interpret a pre-defined format of chanID, use AARAE's
    % utility function readchanID
    
    if size(audio,3)> 1
        Bandfilter = -1;
        % bandID is a vector, usually listing the centre frequencies of the
        % bands
        if isfield(IN,'bandID') % Get the band ID if it exists
            bandID = IN.bandID;
        else
            % asssign ordinal band numbers if bandID does not exist (as an
            % example of how to deal with missing data)
            bandID = 1:size(audio,3);
        end
    end
    
    % The name of the audio input could be useful in generating figures
    % (for example, in the title of a figure).
    if isfield(IN,'name') % Get the AARAE name if it exists
        name = IN.name; % this is a string
    else
        name = ''; % empty string - can be concatenated without any problems
    end
    

    % *********************************************************************
    
    

else
    %input must be an aarae audio structure
    OUT = [];
    return
end
% *************************************************************************





% Check that the required data exists for analysis to run
if ~isempty(audio) && ~isempty(fs)
    
    if Autocrop == 1
            audio = autocropstart_aarae(audio,-20,2);
    elseif Autocrop == 2
            audio = autocropstart_aarae(audio,-20,2);
    end
    
    % the audio's length, number of channels
    [len,chans] = size(audio);
    

    

    
    
    
    % *** EXAMPLES OF CALLING OTHER AARAE FUNCTIONS ***
    % AARAE has many functions in it, some of which may be particularly
    % useful when writing a processor or analyser. Here are a few examples
    % of potentially useful AARAE functions:
    
    % cal_reset_aarae can be used to simplify the process of applying the
    % calibration offset to your audio. In the following case, gain is
    % applied to the audio waveform such that the cal value is changed to 0
    % dB. The cal_reset_aarae function is in the Processors>Basic folder.
    audio = cal_reset_aarae(audio,0,cal);
    
    % autocropstart_aarae is especially useful for impulse response
    % analysis, when data prior to the impulse response start needs to be
    % cropped (deleted). In the following example, the start of the audio
    % wave that is before the first sample -20 dB below the peak is
    % cropped (as suggested in ISO3382-1). The autocropstart_aarae function
    % is in the Processors>Basic folder.

    
    % Aweight, along with similar functions in the Processors>Filters 
    % folder, provides an easy way to implement a weighting filter.
    audio = Aweight(audio,fs);
    
    % octbandfilter_viaFFT is a versatile octave band filterbank. In the
    % example below only the first three input arguments are used, but more
    % can be used if there are particular requirements for the octave band
    % filters' performance. There are other similar functions in the
    % Processors>Filterbanks folder.
    [audio,bandID] = octbandfilter_viaFFT(audio,fs,...
        [125,250,500,1000,2000,4000]);
    
    
    
    % *** CREATING A TABLE OR MULTIPLE TABLES ***
    % You may include tables to display your results using AARAE's
    % disptables.m function - this is just an easy way to display the
    % built-in uitable function in MATLAB. It has several advantages,
    % including:
    %   * automatically sizing the table(s);
    %   * allowing multiple tables to be concatenated;
    %   * allowing concatenated tables to be copied to the clipboard
    %     by clicking on the grey space between the tables in the figure;
    %   * if its output is used as described below, returning data to
    %     the AARAE GUI in a format that can be browsed using a bar-plots.
    %   * and values from tables created this way are written into the log
    %     file for the AARAE session, which provides another way of
    %     accessing the results together with a complete record of the
    %     processing that led to the results.
    
    fig1 = figure('Name',['My results table for ' name]); % use the name if it is not empty
    table1 = uitable('Data',[duration maximum minimum],...
                'ColumnName',{'Duration','Maximum','Minimum'},...
                'RowName',{'Results'});
    disptables(fig1,table1); % see below for a better alternative!
    % It is usually a very good idea to return your table(s) data to AARAE
    % as a function output field. Doing this creates a table leaf in the
    % Results section of the GUI, which can be explored as bar plots within
    % the GUI, and also writes the table contents to the log file (as a
    % comma delimited table that can easily be interpreted by a
    % spreadsheet). To do this, simply use the output of the
    % disptables funtion as follows:
    %       [~,table] = disptables(fig1,table1);
    % And include your table in the output data structure
    %       OUT.tables = table;
    %
    %
    % If you have multiple tables to combine in the figure, you can
    % concatenate them:
    %       disptables(fig1,[table3 table2 table1]);
    % (Note that the concatenation is in descending order - i.e. the first
    % listed table in the square brackets will be displayed at the bottom,
    % and the last listed table will be displayed at the top.)
    %
    % You may export these tables to be displayed as bar plots as if you
    % were doing it for a single table:
    %       [~,tables] = disptables(fig1,[table3 table2 table1]);
    % In the output section, we will use:      OUT.tables = tables;
    %
    % The disptables function will take care of allocating each table to a
    % different barplot, there is no need to generate more than one .tables
    % field to display all your tables.
    
    % *** WRITING TO THE LOG FILE ***
    % There is usually no reason to write to the log file, because a lot of
    % information is automatically written to it (depending on the outputs,
    % including tables). However, if you wish to write to the log file,
    % this is easily done using AARAE's logtext utility:
    logtext('%% This is a test of writing to the log file.\n');
    % For consistency, it is best to start all non-executable lines with %%
    % (which becomes a single % when written to the log file). Hence
    % non-executable lines will be written as Matlab comments.
    
    % *** MAKING PLOTS ***
    % You may also include figures to display your results as plots.
    t = linspace(0,duration,len);
    figure('Name', ['Test figure of ' name]); % use the name of the input audio if it is not empty
    plot(t,audio);
    % All figures created by your function are stored in the AARAE
    % environment under the results box. 
    
    
    
    % *** CREATING A RESULTS LEAF (FOR BIG NON-AUDIO DATA) ***
    % You may output data to be plotted in a variety of charts, including
    % lines, mesh, surf, imagesc, loglog, and others depending on the
    % number of dimensions of your data using the doresultleaf.m function:
    % E.g.:
    %
    doresultleaf(myresultvariable,'Type [units]',{'Time','channels','Frequency'},...
                 'Time',      t,                's',           true,...
                 'channels',  chanID,           'categorical', [],...
                 'Frequency', num2cell(bandfc), 'Hz',          false,...
                 'name','my_results');
    %
    % Input arguments:
    % #1: Your data variable. It can be multidimensional, make sure you
    %     specify what each dimension is.
    % #2: What is your data variable representing? is it level? is it
    %     reverb time? make sure you label it appropriately and assign
    %     units to it, this second argument is a single string.
    % #3: This is a cell array where each cell contains the name of each
    %     dimension.
    %
    % #4: Name of your first dimension. (String)
    % #5: Matrix that matches your first dimension, in this case time.
    % #6: Units of your first dimension. (String)
    % #7: Can this dimension be used as a category? (true, false, [])
    %
    % Replicate arguments 4 to 7 for as many dimensions as your data
    % variable has.
    %
    % The second last input argument is the string 'name', this helps the
    % function identify that the last input argument is the name that will
    % be displayed in AARAEs categorical tree under the results leaf.
    
    
    
    % And once you have your result, you should set it up in an output form
    % that AARAE can understand.
    if isstruct(IN)
        % The output of this function, running within AARAE, is a
        % structure. As a minimum, the structure requires the function 
        % callback. However, in the unusual event that you wish to
        % replicate the input structre before adding to it, it is easiest
        % to do this first:
        
        % *** OUTPUTTING AUDIO ***
        % Most analysers do not output audio. If you wish to output audio,
        % first consider whether your analyser should be a processor
        % instead. If it should be an analyser, then audio can be output as
        % follows:
        % OUT = IN; % You can replicate the input structure for your output
        % OUT.audio = audio; % And modify the fields you processed
        % However, for an analyser, you might not wish to output audio (in
        % which case the two lines above might not be wanted.
        
        % *** FUNCTION CALLBACKS ***
        % The following outputs are required so that AARAE can repeat the
        % analysis without user interaction (e.g. for batch processing),
        % and for writing the log file (and hence for generating
        % workflows).
        OUT.funcallback.name = 'analyser_template.m'; % Provide AARAE
        % with the name of your function 
        OUT.funcallback.inarg = {input_1,input_2,input_3,input_4}; 
        % assign all of the input parameters that could be used to call the
        % function without dialog box to the output field param (as a cell
        % array) in order to allow batch analysing. Do not include the
        % first (audio) input here. If there are no input arguments (apart
        % from the audio input), then use: 
        % OUT.funcallback.inarg = {};
        
        
        % If you generated tables using AARAE's disptables (as described 
        % above, you can attach the tables to the output structure
        %       OUT.tables = tables;
        
        
        % *** OUTPUTTING NEW FIELDS ***
        % Or simply output the fields you consider necessary after
        % processing the input audio data, AARAE will figure out what has
        % changed and complete the structure. But remember, it HAS TO BE a
        % structure if you're returning more than one field:
        
        OUT.properties.duration = duration;
        OUT.properties.maximum = maximum;
        OUT.properties.minimum = minimum;
        % The advantages of providing outputs as subfields of properties
        % are that AARAE has a button that opens a window to display
        % properties, and that properties are also written to the log file
        % of the AARAE session. Outputing fields without making them
        % subfields of properties is possible, but this makes the output
        % data harder to access.
        % Note that the above outputs might be considered to be redundant
        % if OUT.tables is used, as described above. Generally the
        % properties fields output is suitable for small data only (e.g.
        % single values). Tables is best for small and medium data, while a
        % results leaf is is best for big data.
        
        
        % AARAE will only display the output in the main GUI if it includes
        % tables, audio or other fields. It will not display if only
        % function callbacks are returned. Result leaves are not created by
        % the functions main output, but instead are created by calling
        % AARAE's doresultleaf as described above, and these will be
        % displayed in AARAE's main GUI regardless of the function's main 
        % outputs.
    else
        % You may increase the functionality of your code by allowing the
        % output to be used as standalone and returning individual
        % arguments instead of a structure.
        OUT = audio;
    end
    % some additional outputs (not used by AARAE)
    varargout{1} = duration;
    varargout{2} = maximum;
    varargout{3} = minimum;
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