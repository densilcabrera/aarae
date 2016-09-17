function OUT = Tonality_thirdoctmethod(IN,fs,A_weighted,Criterion)

% This function implements the one third octave band analysis methods
% described in: ISO 1996-2 and AS/NZS 2107, which could analysis the input 
% audio waveform in terms of tonality.
% The difference between these two standards is :
% In the ISO 1996-2 method, the input signal will be analyzed directly, 
% however, in AS/NZS 2107 method, the input signal will be A-weighted first.
% Moreover, there are two different criteria for tonality analysis, which
% are "both adjacent bands" and "average level difference". 
% User could select which method they want.
% code by Jianyang Xun & Densil Cabrera
% Version 1.00 May 2015
% *************************************************************************
% The next few lines is to call MATLAB's built-in inputdlg function to 
% allow the user to select which method they want to use for tonality 
% analysis.
if nargin ==1 % If the function is called within the AARAE environment it
              % will have at least one input parameter which is the audio
              % data structure.
    param = inputdlg({'Apply A-weighted [0|1]       ( 0 for ISO 1996-2 ,1 for AS/NZS 2107)';... % In ISO 1996-2 there is
        % no need to apply A-weighting filter to signal need to be analysis, 
        % however, in AS/NZS 2107 standard, it is necessary to apply 
        % A-weighting filter first. 
                      'Peak Criterion on: average level difference [0],both adjacent bands [1]'},...% inputdlg window.
                      'One Third Octave Band Tonality Analysis Settings',... % This is the dialog window title.
                      [1 30],... % Define the number of rows per
                      ...        % input box and the number of character
                      ...        % spaces that each box can display at once
                      ...        % per row. 
                      {'1';'1'}); % And the preset answers of the dialog.

    param = str2num(char(param)); % Since inputs are usually numbers it's a
                                  % good idea to turn strings into numbers.
                             

    if length(param) < 1, param = []; end % Check that the user 
                                          % has input all the required
                                          % fields.
    if ~isempty(param) % If they have input something,  then assign the dialog's
                       % inputs to this function's input parameters.
        A_weighted = param(1);
  
        Criterion = param(2);
        
    else
        % get out of here if the user presses 'cancel'
        OUT = [];
        return
    end
else
    param = [];
end

% *************************************************************************
if isstruct(IN) % Because AARAE is a Matlab structure, hence we need to
    % test whether this function is being called within AARAE environment.
    %This step is being used to test whether it is struct.
  
    %The following utility function call
    % reduces the number of dimensions in the audio field of the input
    % structure if they are greater than 3. For further information 
    % see the comments in "choose_from_higher_dimensions"
    %in AARAE's utilities directory.
    IN = choose_from_higher_dimensions(IN,3,1);
    
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
    
    
  
    if isfield(IN,'cal') %If the user select Calibrate, it will be applied
        %on the input signal 
        audio = cal_reset_aarae(audio, 0, IN.cal);
    end
    
    % chanID is a cell array of strings describing each channel
    if isfield(IN,'chanID') % Get the channel ID if it exists
        chanID = IN.chanID;
    else
        % or make a chanID using AARAE's utility function
        chanID = makechanID(size(audio,2),0);
    end
    
    
elseif ~isempty(param) || nargin > 1
    % If for example you want to enable your function to
    % run as a standalone MATLAB function, you can use
    % the IN input argument as an array of type double
    % and assign it to the audio your function is going
    % to process.
    audio = IN;
end
% *************************************************************************

% Check that the required data exists for analysis to run
if ~isempty(audio) && ~isempty(fs) 

    bands = size(audio,3);
    
    if bands > 1
        audio = sum(audio,3); % mixdown bands if multiband
    end
    
    if A_weighted == 1
    audio = Aweight(audio,fs);
    % Aweight, use AARAE's A-weighting filter: in Processors/Filters
    end %If user select to Apply A-weighting filter on input signal
    
    % Apply 1/3-octave band filterbank
    frequencies = [25,31.5,40,50,63,80,100,125,160,200,250,315,400,...
      500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000];
    
    audio = thirdoctbandfilter(audio,fs,frequencies,1);
   % call the function of "thirdoctbandfilter" in
   % AARAE-processor-filterbanks
    L = 10*log10(mean(mean(audio.^2,1),2));
    %transfer the filtered audio into dB domain
    L = permute(L,[1,3,2]); % rearrange the dimension
    
    
    
    % ****************** 
    if Criterion ==1 
    % if the user select the both adjacent band criterion
    % According to AS/NZS 2107:2000 and IS0 1996-2.
    % If a frequency component is a tonal component, 
    % it should meet the requirement that the time-average sound pressure 
    % level in the frequency band is higher than the time average sound     
    % pressure levels of both adjacent one-third octave bands by 
    % a certain level difference.

    % 15 dB in the low-frequency one-third-octave bands (25 Hz to 125 Hz);
    % 8 dB in middle-frequency bands (160 Hz to 400 Hz);
    % 5 dB in high-frequency bands (500 Hz to 10 000 Hz).
  
    Low_band = L(1:8); % (25 Hz to 125 Hz)
    Mid_band = L(9:13); %(160 Hz to 400 Hz)
    High_band = L(14:27);%(500 Hz to 10 000 Hz)
    % Hence these three lines above is to divide input signal to three
    % freuqency bands
    [~,locs_low] = findpeaks(Low_band,'Threshold',15);
    [~,locs_mid] = findpeaks(Mid_band,'Threshold',8);
    [~,locs_high] = findpeaks(High_band,'Threshold',5);
    % Matlab built-in fucntion "findpeaks" is used in this part, and the
    % threshold is the level diferences according to the standard.
 % ****************** ****************** ******************
    % Because the Matlab built-in function "Findpeaks" can not detect
    % whether the first and the last components are peaks. 
    % Hence this section will navigate whether the first value and the last 
    % value are tonal components.


    First_low = L(1)-L(2);
    if First_low>=15    % if the first vlaue in low band is higher than
        % the second value (more than 15 dB) 
        Final_locs_low(1) = 1;
        % the first component is a tonal component
        for ii=2:size(locs_low,2)+1
            Final_locs_low(ii) = locs_low(ii-1);
            % and the locations found by "findpeaks" need to be changed
        end
    else
        Final_locs_low=locs_low;
        % If not, it won't change anything
    end
    
    %if the last is higher than the previous one (more than 15 dB)
    
    Last_low = L(8)-L(7);
    if Last_low>=15
        Final_locs_low(size(Final_locs_low,2)+1)=8;
        % the last component is a tonal component
    else
        Final_locs_low=Final_locs_low;
    end
    %
    if ~isempty(Final_locs_low) %if there exists peaks
    Final_pks_low = L(Final_locs_low);
    Final_lowpeak_frequency = frequencies (Final_locs_low);
    else
        Final_pks_low = [];
        Final_lowpeak_frequency = [];
    end
        
    

    
    % same method to detect mid frequency band
    
    First_mid = L(9)-L(10);
    if First_mid>=8
        Final_locs_mid(1) = 9;
        
        for ii=2:size(locs_mid,2)+1
            Final_locs_mid(ii) = locs_mid(ii-1)+8;
        end
    else
        Final_locs_mid=locs_mid+8;
    end
    
    %if the last is higher than the previous one
    
    Last_mid = L(13)-L(12);
    if Last_mid>=8
        Final_locs_mid(size(Final_locs_mid,2)+1)=13;
    else
        Final_locs_mid=Final_locs_mid;
    end
    %
    if ~isempty(Final_locs_mid)
    Final_pks_mid = L(Final_locs_mid);
    Final_midpeak_frequency = frequencies (Final_locs_mid);
    else
        Final_pks_mid = [];
    Final_midpeak_frequency = [];
    end
    %
 
      
    % same method to detect high frequency band
    
    First_high = L(14)-L(15);
    if First_high>=5
        Final_locs_high(1) = 14;
        
        for ii=2:size(locs_high,2)+1
            Final_locs_high(ii) = locs_high(ii-1)+13;
        end
    else
        Final_locs_high = locs_high+13;
    end
    
    %if the last is higher than the previous one
    
    Last_high = L(27)-L(26);
    if Last_high>=5
        Final_locs_high(size(Final_locs_high,2)+1)=27;
    else
        Final_locs_high=Final_locs_high;
    end
    %
    if ~isempty(Final_locs_high)
        Final_pks_high = L(Final_locs_high);
        Final_highpeak_frequency = frequencies (Final_locs_high);
    else
        Final_pks_high = [];
        Final_highpeak_frequency = [];
    end

    % combine the peak components together
    Final_pks = [Final_pks_low Final_pks_mid Final_pks_high];
    Final_frequency = [Final_lowpeak_frequency Final_midpeak_frequency Final_highpeak_frequency];
    Final_locs = [Final_locs_low Final_locs_mid Final_locs_high];
    
    else % if select the average level difference criterion
    Low_band = L(1:8);
    Mid_band = L(9:13);
    High_band = L(14:27);
    locs_low = [];
    locs_mid = [];
    locs_high = [];% set the initial peak location to be empty
    
    j=1; % initial location index is 1 
    for i=2:7 % for low freuqency band, except 25 Hz and 125 Hz (the boundary£©

      if [L(i)-L(i-1)+L(i)-L(i+1)]/2>=15 % for average method, if the componenet 
        % is higher than the both adjacent average level for a certain
        % value, it shoud be a tonal component.
        locs_low(j) = [i]; % if it is a tonal component, it will be set to the peak location index 
           j= j+1; 
      end
    end
    
    % same method above to detect the first value
    First_low = L(1)-L(2);
    if First_low>=15
        Final_locs_low(1) = 1;
        
        for ii=2:size(locs_low,2)+1
            Final_locs_low(ii) = locs_low(ii-1);
        end
    else
        Final_locs_low=locs_low;
    end
    
    %if the last is higher than the previous one
    
    Last_low = L(8)-L(7);
    if Last_low>=15
        Final_locs_low(size(Final_locs_low,2)+1)=8;
    else
        Final_locs_low=Final_locs_low;
    end
    %
    if ~isempty(Final_locs_low)
    Final_pks_low = L(Final_locs_low);
    Final_lowpeak_frequency = frequencies (Final_locs_low);
    else
        Final_pks_low = [];
        Final_lowpeak_frequency = [];
    end

    
    % for mid frequency band 

    jj=1;
    for i=10:12
    if [L(i)-L(i-1)+L(i)-L(i+1)]/2>=8
        locs_mid (jj)= [i];
        jj = jj+1;
    end
     end   
  
     First_mid = L(9)-L(10);
    if First_mid>=8
        Final_locs_mid(1) = 9;
       
        for ii=2:size(locs_mid,2)+1
            Final_locs_mid(ii) = locs_mid(ii-1)
            %+8;
        end
          else
        Final_locs_mid=locs_mid;
            end
    
    %if the last is higher than the previous one
    
    Last_mid = L(13)-L(12);
    if Last_mid>=8
        Final_locs_mid(size(Final_locs_mid,2)+1)=13;
    else
        Final_locs_mid=Final_locs_mid;
    end
    %
    if ~isempty(Final_locs_mid)
    Final_pks_mid = L(Final_locs_mid);
    Final_midpeak_frequency = frequencies (Final_locs_mid);
    else
        Final_pks_mid = [];
    Final_midpeak_frequency = [];
    end

  
  % for high frequency band 
  
  jjj = 1
  for i=15:26
    if [L(i)-L(i-1)+L(i)-L(i+1)]/2>=5
        locs_high(jjj) = [i];
        jjj = jjj+1;
    end
  end    

    % high
    
    First_high = L(14)-L(15);
    if First_high>=5
        Final_locs_high(1) = 14;
        
        for ii=2:size(locs_high,2)+1
            Final_locs_high(ii) = locs_high(ii-1)
            %+13;
        end
    else
        Final_locs_high = locs_high
        %+13;
    end
    
    %if the last is higher than the previous one
    
    Last_high = L(27)-L(26);
    if Last_high>=5
        Final_locs_high(size(Final_locs_high,2)+1)=27;
    else
        Final_locs_high=Final_locs_high;
    end
    %
    if ~isempty(Final_locs_high)
        Final_pks_high = L(Final_locs_high);
        Final_highpeak_frequency = frequencies (Final_locs_high);
    else
        Final_pks_high = [];
        Final_highpeak_frequency = [];
    end
    % 
    % High_peak = [Final_pks_high;Final_highpeak_frequency];
    
    % combine all the peaks together 
    Final_pks = [Final_pks_low Final_pks_mid Final_pks_high];
    Final_frequency = [Final_lowpeak_frequency Final_midpeak_frequency Final_highpeak_frequency];
    Final_locs = [Final_locs_low Final_locs_mid Final_locs_high];
    end
    % ****************** 
    
    %plot the results
    
           figure('name', '1/3-Octave Band Tonality Analyisis')
           
           ymax = 10*ceil(max(L+5)/10);
           ymin = 10*floor(min(L)/10);

            width = 0.5;
      
            
            
            bar(1:length(frequencies),L,width,'FaceColor',[0,0.7,0],...
                'EdgeColor',[0,0,0],'DisplayName', 'Non-tonal components','BaseValue',ymin);
            % plot the non-tonal components in green 
            hold on
            
            peakbars = nan(1,length(frequencies));

            for k = 1:length(frequencies)
                if ~isempty(find(Final_locs==k, 1))
                    peakbars(k) = L(k);
                end
            end
             % plot the tonal components in red
            bar(1:length(frequencies),peakbars,width,'stacked','FaceColor',[1,0,0], ...
                'EdgeColor',[0,0,0],'DisplayName', 'Tonal components','BaseValue',ymin);
            
            
            
            %hold on

            % x-axis
            set(gca,'XTick',1:length(frequencies),'XTickLabel',num2cell(frequencies))

            xlabel('1/3-Octave Band Centre Frequency (Hz)')


            % y-axis
            
            
            ylabel('Level (dB)')
            ylim([ymin ymax])


          legend('show','Location','EastOutside');
                  
            

            for k = 1:length(frequencies)
                if ~isempty(find(Final_locs==k, 1))
                text(k-0.25,ymax-(ymax-ymin)*0.025, ...
                    num2str(round(L(k)*10)/10),'Color',[1,0,0])
                else
                    text(k-0.25,ymax-(ymax-ymin)*0.025, ...
                    num2str(round(L(k)*10)/10),'Color',[0,0.7,0])
                end
            end
    
    
            
    
    
    
    if isstruct(IN)
        
 
        OUT.funcallback.name = 'Tonality_thirdoctmethod.m'; % Provide AARAE
        % with the name of your function
        OUT.funcallback.inarg = {fs,A_weighted,Criterion};

    else
        % You may increase the functionality of your code by allowing the
        % output to be used as standalone and returning individual
        % arguments instead of a structure.
        %OUT = audio;
        OUT = [];
    end
    
else
    % AARAE requires that in case that the user doesn't input enough
    % arguments to generate audio to output an empty variable.
    OUT = [];
end

%**************************************************************************
% Copyright (c) <2015>, <Jianyang Xun & Densil Cabrera>
% All rights reserved
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