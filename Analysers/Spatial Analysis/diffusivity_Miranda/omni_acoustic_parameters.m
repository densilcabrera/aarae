function omniparameters = omni_acoustic_parameters(IRs, fs, omni_setup)

if nargin < 3;
    omni_setup = get_omni_setup;
end

if omni_setup.calculation == 0;
    omniparameters = [];
else
%% Initialize variables for calculations

analysis_length_insamps = round((omni_setup.analysis_length./1000).*fs);
num_analys_windows = ceil(length(IRs) / (analysis_length_insamps*(1-(omni_setup.overlap/100)))); 
overlap = omni_setup.overlap;
IRslength = size(IRs,1);
omni_signal = IRs(:,1);

%% Signal filtering

[hpF lpF] = band_pass_filter(omni_setup.wideband_edgebands(1), omni_setup.wideband_edgebands(2), fs);
wideband_sigs = filter(lpF,filter(hpF, omni_signal));

number_of_bands = length(omni_setup.center_bands);
octavebandsigs = zeros(size(omni_signal,1),number_of_bands);

octavebandsigs = octbandfilter(IRs,fs,omni_setup.center_bands,1);

omniparameters.bands = omni_setup.center_bands;

%% Magnitude calc

for I = 0:num_analys_windows-1
  analysisstart = round((I .* analysis_length_insamps + 1)-(analysis_length_insamps.*I.*(overlap/100))); 
  
  analysisend = round(analysisstart + analysis_length_insamps);
  if analysisend > IRslength; 

      analysisend = IRslength;
  end
  
  omniparameters.analysis_times(I+1,1) = analysisstart./fs;
  
  omniparameters.magnitude(I+1,1) = sqrt(sum(wideband_sigs(analysisstart:analysisend,:).^2));
  omniparameters.magnitudedB(I+1,1) = mag2db(omniparameters.magnitude(I+1,1)); 
  
end

for I = 0:num_analys_windows-1
    for J = 1:4;
      analysisstart = round((I .* analysis_length_insamps + 1)-(analysis_length_insamps.*I.*(overlap/100))); 

      analysisend = round(analysisstart + analysis_length_insamps);
      if analysisend > IRslength; 

          analysisend = IRslength;
      end


      omniparameters.magnitude_oct_band(I+1,J) = sqrt(sum(octavebandsigs(analysisstart:analysisend,J).^2));
      omniparameters.magnitudedB_oct_band(I+1,J) = mag2db(omniparameters.magnitude_oct_band(I+1,J)); 
    end
end

if omni_setup.normalize == 1;
    omniparameters.magnitudedB = omniparameters.magnitudedB-max(omniparameters.magnitudedB);
    omniparameters.magnitudedB_oct_band = omniparameters.magnitudedB_oct_band-max(max(omniparameters.magnitudedB_oct_band));
end

%% ST Parameters

ten_ms = (10./1000).*fs;
twenty_ms = (20./1000).*fs;
one_h_ms = (100./1000).*fs;
one_sec = fs;

for I = 1:4;
    [~, firstarrival] = max(octavebandsigs(:,I));
    direct_sound(I) = sum(octavebandsigs(firstarrival:firstarrival+ten_ms,I).^2);
    early_part(I) = sum(octavebandsigs(firstarrival+twenty_ms:firstarrival+one_h_ms,I).^2);
    if firstarrival+one_sec > IRslength;
        disp('File length smaller than 1 sec for STlate. Trimming to max file length')
        late_part(I) = sum(octavebandsigs(firstarrival+one_h_ms:IRslength,I).^2);
    else
        late_part(I) = sum(octavebandsigs(firstarrival+one_h_ms:firstarrival+one_sec,I).^2);
    end
    omniparameters.STearly_band(I) = 10.*log10(early_part(I)./ direct_sound(I));
    omniparameters.STlate_band(I) = 10.*log10(late_part(I)./ direct_sound(I));
end
    

omniparameters.STearly = mean(omniparameters.STearly_band);
omniparameters.STlate = mean(omniparameters.STlate_band);

%% Reverberation time

nmbSmp = size(octavebandsigs,1) ;

omniparameters.SchCur = zeros(nmbSmp,number_of_bands) ;
    
for I = 1 : number_of_bands;
    cur = 10*log10(flipud(cumsum(flipud(octavebandsigs(:,I).^2)))) ;
    omniparameters.SchCur(:,I) = cur - max(cur) ;
end

    
for I = 1 : number_of_bands;

    fst = find(omniparameters.SchCur(:,I)<=-05,1,'first') ;
    lst = find(omniparameters.SchCur(:,I)<=-35,1,'first') ;
    lst25 = find(omniparameters.SchCur(:,I)<=-25,1,'first') ;
    coe = LinearRegression((fst:lst)'/fs,omniparameters.SchCur(fst:lst,I)) ;
    coe25 = LinearRegression((fst:lst25)'/fs,omniparameters.SchCur(fst:lst25,I)) ;

    omniparameters.T30(I) = -60/coe(1) ;
    omniparameters.T20(I) = -60/coe25(1) ;

end
    





end

% Linear regression subroutine
function coe = LinearRegression(tme,eng)
    
    % Number of points to fit
    nmb = length(tme) ;

    % Matrix we need to invert
    Mat = [tme(:) ones(nmb,1)] ;
    
    % Linear coefficients
    coe = pinv(Mat) * eng ;


    
