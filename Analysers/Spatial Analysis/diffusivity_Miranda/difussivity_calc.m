function diffusivity = difussivity_calc(IRs, fs, hoaFmt, diffusivity_setup)

if nargin < 4
    diffusivity_setup = get_diffusivity_setup;
end

if diffusivity_setup.calculation == 0;
    diffusivity = [];
else
%% Initialize variable for calculations

analysis_length_insamps = round((diffusivity_setup.analysis_length./1000).*fs);
num_analys_windows = ceil(length(IRs) / (analysis_length_insamps*(1-(diffusivity_setup.overlap/100)))); 
overlap = diffusivity_setup.overlap;
IRslength = size(IRs,1);

%% Wideband calc

[hpF lpF] = band_pass_filter(diffusivity_setup.wideband_edgebands(1), diffusivity_setup.wideband_edgebands(2), fs);

wideband_sigs = filter(lpF,filter(hpF, IRs));

if strcmp(diffusivity_setup.method,'Gover') == 1;
   for I = 0:num_analys_windows-1
      analysisstart = round((I .* analysis_length_insamps + 1)-(analysis_length_insamps.*I.*(overlap/100))); 
      analysisend = round(analysisstart + analysis_length_insamps);
      if analysisend > IRslength; 
          analysisend = IRslength;
      end
      diffusivity.analysis_times(I+1,1) = analysisstart./fs;
      diffusivity.diff_wideband(I+1,1) = GoverDiffuseness(wideband_sigs(analysisstart:analysisend,:),hoaFmt);
   end
elseif strcmp(diffusivity_setup.method,'HOA') == 1;
   for I = 0:num_analys_windows-1
      analysisstart = round((I .* analysis_length_insamps + 1)-(analysis_length_insamps.*I.*(overlap/100))); 
      analysisend = round(analysisstart + analysis_length_insamps);
      if analysisend > IRslength; 
          analysisend = IRslength;
      end
      diffusivity.analysis_times(I+1,1) = analysisstart./fs;
      diffusivity.diff_wideband(I+1,1) = HoaDiffuseness(wideband_sigs(analysisstart:analysisend,:),hoaFmt);
   end    
end

%% Band calculation

number_of_bands = length(diffusivity_setup.center_bands);

octavebandsigs = zeros(size(IRs,1),size(IRs,2),number_of_bands);

octavebandsigs = octbandfilter(IRs,fs,diffusivity_setup.center_bands,1);

diffusivity.bands = diffusivity_setup.center_bands;

if strcmp(diffusivity_setup.method,'Gover') == 1;
    for J = 1:number_of_bands;
       for I = 0:num_analys_windows-1
          analysisstart = round((I .* analysis_length_insamps + 1)-(analysis_length_insamps.*I.*(overlap/100))); 
          analysisend = round(analysisstart + analysis_length_insamps);
          if analysisend > IRslength; 
              analysisend = IRslength;
          end
          diffusivity.diff_octband(I+1,J) = GoverDiffuseness(octavebandsigs(analysisstart:analysisend,:,J),hoaFmt);
       end
    end
elseif strcmp(diffusivity_setup.method,'HOA') == 1;
    for J = 1:number_of_bands;
       for I = 0:num_analys_windows-1
          analysisstart = round((I .* analysis_length_insamps + 1)-(analysis_length_insamps.*I.*(overlap/100))); 
          analysisend = round(analysisstart + analysis_length_insamps);
          if analysisend > IRslength; 
              analysisend = IRslength;
          end
          diffusivity.diff_octband(I+1,J) = HoaDiffuseness(octavebandsigs(analysisstart:analysisend,:,J),hoaFmt);
       end
    end
end
end
