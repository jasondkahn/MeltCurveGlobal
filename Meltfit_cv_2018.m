function Fitvals = Meltfit_cv_2018(filestem,varargin)
% Usage example: >> Meltfit_2018('DN2') or >> Meltfit_2018('DN2',6)
% varargin is Tmin, Tmax, Scale_factor, Debug flag
%  or >> Meltfit_cv_2018('DN2',6,80,100,2)
% Scale factor = 1000 means the low conc melt is weighted 1000x more
% Leave it blank to let the program evaluate the scale factor
% variable args are Tmin,Tmax,scale factor,debug
% The file filestem.xls(x) must contain columns as follows:
% Col 1:  temperature array for high conc dsDNA melt, in Centigrade
% Col 2:  absorbance array for high conc melt
% Col 3:  temperature array for low conc dsDNA melt, in C
% Col 4:  absorbance array for low conc melt
% Col 5:  temperature array for ssDNA #1 reference, in C
% Col 6:  absorbance array for ssDNA #1
% Col 7:  temperature array for ssDNA #2 reference, in C
% Col 8:  absorbance array for ssDNA #2
% Col 9:  Row 1: epsilon(Ref temp) for ssDNA #1  (Ref temp default is 20 °C)
%         Row 2: epsilon(Ref temp) for ssDNA #2
% Col 10: Row 1: Melt volume in uL, should be the same for all (is this needed?)
%         Row 2: Volume of the higher conc stock (typically 40 micomolar stock) 
%                that was used for ss1 -> hi conc melt  (sample in columns
%                1 and 2)
%         Row 3: Volume of ss1 stock -> lo conc melt (sample in columns 3 and 4)
%         Row 4: Volume of ss1 stock -> ss1 ref sample (sample in columns 5 and 6)
%         Rows 5-7 are the corresponding ss2 -> hi, ss2 -> lo, ss2 -> ss2 ref.
%              Omit the latter three if they are the same as for ss1.
% Col 11: Optional temperature array for baseline absorbance
% Col 12: Optional absorbance array for baseline 
%
% There can be as many lines of text (apparently there _must_ be a couple of
% lines??) at the top as you like, but no numbers at all!
% There can't be anything to the right of the columns containing numbers,
% and no text among the numbers in columns 1-12. The printout will emerge
% as a png.
% Input is tricky if there is any taint of a number above the data. That
% will be the only thing read, everything else will be garbage. Absymally
% stupid. One solution is to put everything but the data on a separate
% sheet, but make sure that the first sheet has the data.
%
% Set default parameters and constants.
%
pathlength = 1; % (cm) Change by hand here if needed
RefT = 20; % reference temp for cocncentration calculations
Ratio_tolerance = 0.05 ;  % tolerance for expected vs. observed Hi T absorbances.
Tspace = 3;  % Range of temps over which to average absorbances near reference T
input_Tmax = 90; % default max temp
debug = 0; % default do not show debugging
%
% Change default values for Tmin and Tmax if user provides them
% Can provide Tmin only or both Tmin and Tmax, not Tmax only 
% Usually the bad data is at the very lowest temp.
%
numrequired = 0;
numpossible = 4;
numvarargs = length(varargin);
if numvarargs < numrequired
    error('Meltfit_2018:NotEnoughInputs', ...
        'Meltfit_2018 requires at least %d inputs',numrequired);
end
if numvarargs > numpossible
    error('Meltfit_2018:TooManyInputs', ...
        'Meltfit_2018 can handle at most %d inputs',numpossible);
end
% now put these defaults into the valuesToUse cell array, 
scale_flag = 0;
if (numvarargs ~= 0)
    optargs(1:numvarargs) = varargin;
    input_Tmin = optargs{1};
    if (~isempty(input_Tmin))
        Tmin = input_Tmin;
        if (Tmin > RefT)
            disp(['Meltfit needs data including ref temp ',num2str(RefT),...
            ' °C, you entered Tmin = ', num2str(Tmin)]);
            return
        end
        temp_Tmin = Tmin;
    else
        temp_Tmin = 0;
    end
    if (numvarargs >= 2)
        input_Tmax = optargs{2};
        if (input_Tmax < temp_Tmin)
            disp('Please enter filestem, Tmin, Tmax > Tmin, num iterations, debug flag');
            return
        elseif (input_Tmax < RefT)
            disp(['Meltfit needs data including ref temp',num2str(RefT),...
                ' °C, you entered Tmax = ', num2str(input_Tmax)]);
        return
        end
    end
%     if (numvarargs >= 3)
%         if (~isempty(optargs{3}))
%             numtries = optargs{3};
%         end
%         if (numtries < 1 || mod(numtries,1) ~= 0)
%             disp(['Number of iterations = ',num2str(numtries),...
%                 ' must be an integer >= 1']);
%             return
%         end
%     end
    if (numvarargs >= 3) 
        % Define a scale factor e.g. to look at one melt alone by setting
        % it very high or low.
        if (~isempty(optargs{3}))
            scale_factor = optargs{3};
            scale_flag = 1;
        end
    end
    if (numvarargs >= 4) 
        if (~isempty(optargs{4}))
            debug = optargs{4};
        end
        if (debug < 0 || mod(debug,1) ~= 0)
            disp(['Debug number = ',num2str(debug),...
                ' must be zero or a positive integer >=1']);
            return
        end
    end
end
%
% Now slurp in the data from the Excel file
% Snippet from https://www.mathworks.com/matlabcentral/answers/...
%    49414-check-if-a-file-exists

if (exist([filestem '.xlsx'],'file'))
  % File exists.  Do stuff....
  [A , ~]= xlsread([filestem '.xlsx']);
elseif (exist([filestem '.xls'],'file'))
  % File exists.  Do stuff....
  [A , ~]= xlsread([filestem '.xls']);
else
  % File does not exist.
  disp(['Failing: file does not exist: ',filestem, '.xl*'])
  return
end
% Bring in Absorbance vs. Temperature arrays
% The isnan stuff strips out NaN's so that the reference array and the melt
% array need not be the same length.
T_hiCmelt = A(~isnan(A(:,1)),1);
if (debug > 1)
    T_hiCmelt(1:5)
end
% Print out to check if there is an input problem
% Or remove trailing semicolons below to print out results
% Better yet just do [A,B] = xlsread('filename') from the command line and
% look at the results
Abs_hiCmelt = A(~isnan(A(:,2)),2);
if (T_hiCmelt(end) < T_hiCmelt(1))
    T_hiCmelt = flipud(T_hiCmelt);
    Abs_hiCmelt = flipud(Abs_hiCmelt);
end
T_loCmelt = A(~isnan(A(:,3)),3);
Abs_loCmelt = A(~isnan(A(:,4)),4);
if (T_loCmelt(end) < T_loCmelt(1))
    T_loCmelt = flipud(T_loCmelt);
    Abs_loCmelt = flipud(Abs_loCmelt);
end
Tlow = min([T_hiCmelt ; T_loCmelt]);
Thigh = max([T_hiCmelt ; T_loCmelt]);
T_ss1ref = A(~isnan(A(:,5)),5);
Abs_ss1ref = A(~isnan(A(:,6)),6);
if (T_ss1ref(end) < T_ss1ref(1))
    T_ss1ref = flipud(T_ss1ref);
    T_ss1ref = flipud(T_ss1ref);
end
T_ss2ref = A(~isnan(A(:,7)),7);
Abs_ss2ref = A(~isnan(A(:,8)),8);
if (T_ss2ref(end) < T_ss2ref(1))
    T_ss2ref = flipud(T_ss2ref);
    Abs_ss2ref = flipud(Abs_ss2ref);
end
% If the absorbance blank columns are not empty subtract off interpolated
% absorbance
if size(A,2) > 10
    T_blank = A(~isnan(A(:,11)),11);
    Abs_blank = A(~isnan(A(:,12)),12);
    if (T_blank(end) < T_blank(1))
        T_blank = flipud(T_blank);
        Abs_blank = flipud(Abs_blank);
    end
    Abs_hiCmelt = Abs_hiCmelt - interp1(T_blank,Abs_blank,T_hiCmelt,'pchip');
    Abs_loCmelt = Abs_loCmelt - interp1(T_blank,Abs_blank,T_loCmelt,'pchip');
    Abs_ss1ref = Abs_ss1ref - interp1(T_blank,Abs_blank,T_ss1ref,'pchip');
    Abs_ss2ref = Abs_ss2ref - interp1(T_blank,Abs_blank,T_ss2ref,'pchip');
    disp('Subtracted background')
else
    disp('Did not find blank columns 11 and 12, using raw absorbances')
end
if ~exist('scale_factor','var')
    scale_factor = sqrt(mean(Abs_hiCmelt)/mean(Abs_loCmelt));
    disp(['Computed scale factor as ',num2str(scale_factor)])
    if(~isreal(scale_factor))
        disp('Scale factor calculated as a complex #. Probably a baseline subtraction issue. Bailing.')
        return
    end
end
% These numbers have to be there. Should check? Just print them out.
Epsil_ss1_RefT = A(1,9);
Epsil_ss2_RefT = A(2,9);
disp(['Extinction coeffs ',num2str(Epsil_ss1_RefT), ' and ',num2str(Epsil_ss2_RefT)])
% A(1,10) is the melt volume Meltvol. Uncomment line below if it becomes
% relevant
%Meltvol = A(1,10);
%~ = A(1,10);
% Turns out meltvol is not needed because we get []'s from reference samples
V_ss1_Hi = A(2,10);
V_ss1_Lo = A(3,10);
V_ss1_ssref1 = A(4,10);
if (~isnan(A(5,10)))
    V_ss2_Hi = A(5,10);
else
    V_ss2_Hi = V_ss1_Hi;
end
if (~isnan(A(6,10)))
    V_ss2_Lo = A(6,10);
else
    V_ss2_Lo = V_ss1_Lo;
end
if (~isnan(A(7,10)))
    V_ss2_ssref2 = A(7,10);
else
    V_ss2_ssref2 = V_ss1_ssref1;
end
%
if (~exist('Tmin','var'))
    % Why should T be an integer??
    Tmin = max ([ min(T_hiCmelt) , min(T_loCmelt) , min(T_ss1ref) , min(T_ss2ref)]);
end
%if (~exist('Tmax','var'))
    %Tmax = floor( min ([ input_Tmax , max(T_hiCmelt) , max(T_loCmelt) , 
    % max(T_ss1ref) , max(T_ss2ref)]));
    Tmax = min([input_Tmax, max(T_hiCmelt),max(T_loCmelt),max(T_ss1ref),max(T_ss2ref)]);
%end
disp(['T range is ',num2str(Tmin),' - ',num2str(Tmax)]);
if (Tmin > Tmax)
    disp('Tmin > Tmax')
    return
end
if (RefT < Tmin || RefT > Tmax)
    disp(['RefT ',num2str(RefT),' is outside the T range ',...
        num2str(Tmin),' - ',num2str(Tmax)]);
    return
end
% Trim down to values between Tmin and Tmax


% Now calculate concentrations and evaluate whether they are reasonable
% Need to pull out values nearest RefT °C from sssref absorbance arrays.
T_ss1_near_RefT = find (T_ss1ref < (RefT + Tspace) & T_ss1ref > (RefT - Tspace));
T_ss2_near_RefT = find (T_ss2ref < (RefT + Tspace) & T_ss2ref > (RefT - Tspace));
% Interpolate from data near RefT; allow extrapolation as long as there is
% data from within Tspace of RefT. i.e. could use 85 C if there are at
% least two data points between (85 - Tspace) and 84.98.
p3 = polyfit(T_ss1ref(T_ss1_near_RefT),Abs_ss1ref(T_ss1_near_RefT),1);
Abs_ss1ref_RefT = polyval(p3,RefT);
p4 = polyfit(T_ss2ref(T_ss2_near_RefT),Abs_ss2ref(T_ss2_near_RefT),1);
Abs_ss2ref_RefT = polyval(p4,RefT);
% Best guess reference ssDNA concentrations directly from absorbance on ref
% samples.
Conc1_ref = Abs_ss1ref_RefT/(Epsil_ss1_RefT*pathlength);
disp(['Best initial guess [ss1] = ',num2str(Conc1_ref)])
Conc2_ref = Abs_ss2ref_RefT/(Epsil_ss2_RefT*pathlength);
disp(['Best initial guess [ss2] = ',num2str(Conc2_ref)])
% To calculate the concentration of ssDNA in the melts takes some
% calculation that depends on initial estimates of Tm, so do this in the
% fitting routine.
%
% Now calculate the predicted Abs in the melt sample at the high end of the T range
% If the true values are close to the predicted values, just scale the
% volumes slightly to bring them into accordance. If not, flag an error.
% Use a high ref temp that is just inside the range of values from the ds
% melt
Tmax_interp_min = 5;  % Minimum range of temp for interpolation of Hi T absorbance
%
% Pull out values between Tmin and Tmax
%T_hiC = find(T_hiCmelt(T_hiCmelt > Tmin & T_hiCmelt < Tmax))
T_hiC = find(T_hiCmelt > Tmin & T_hiCmelt < Tmax);
T_loC = find(T_loCmelt > Tmin & T_loCmelt < Tmax);
% At this point we have all the concentrations, extinction coefficients,
% etc. Now pruned the data to Tmin < T < Tmax
Abs_hiCmelt = Abs_hiCmelt(T_hiC);
Abs_loCmelt = Abs_loCmelt(T_loC);
T_hiCmelt = T_hiCmelt(T_hiC);
T_loCmelt = T_loCmelt(T_loC);
L = length(Abs_hiCmelt);
Abs_der_Hi = zeros(1,L-2);
for I=2:L-1
    Abs_der_Hi(I-1) = (Abs_hiCmelt(I-1)-Abs_hiCmelt(I+1))/...
        (T_hiCmelt(I-1) - T_hiCmelt(I+1));
end
L = length(Abs_loCmelt);
Abs_der_Lo = zeros(1,L-2);
for I=2:L-1
    Abs_der_Lo(I-1) = (Abs_loCmelt(I-1)-Abs_loCmelt(I+1))/...
        (T_loCmelt(I-1) - T_loCmelt(I+1));
end

% To figure out exact concentrations, we need some estimate of how much of
% the high temp data should match the ss reference. We can just guess at
% the concntrations of Watson and Crick based on input volumes for this
% purpose.

C1_tot_Himelt_init = (V_ss1_Hi/V_ss1_ssref1)*(Conc1_ref);
C2_tot_Himelt_init = (V_ss2_Hi/V_ss2_ssref2)*(Conc2_ref);
Wa_tot_hi_est = max(C1_tot_Himelt_init,C2_tot_Himelt_init);
Cr_tot_hi_est = min(C1_tot_Himelt_init,C2_tot_Himelt_init);
C1_tot_Lomelt_init = (V_ss1_Lo/V_ss1_ssref1)*(Conc1_ref);
C2_tot_Lomelt_init = (V_ss2_Lo/V_ss2_ssref2)*(Conc2_ref);
Wa_tot_lo_est = max(C1_tot_Lomelt_init,C2_tot_Lomelt_init);
Cr_tot_lo_est = min(C1_tot_Lomelt_init,C2_tot_Lomelt_init);

R = 1.987;
[maxslope,maxindex] = max(Abs_der_Hi);
guesstm = T_hiCmelt(maxindex);
disp(['Initial guess for Tm is ',num2str(guesstm),' C'])
guessDH0 = -6*R*(maxslope/(max(Abs_hiCmelt)-min(Abs_hiCmelt)))*(guesstm+273.15)^2;
disp(['Initial guess for enthalpy is ',num2str(guessDH0),' cal/mole'])
guessDS0 = (guessDH0/(guesstm+273.15)) - R*log(Wa_tot_hi_est - Cr_tot_hi_est/2);
disp(['Initial guess for entropy is ',num2str(guessDS0),' cal/moleK'])
thermvals = [guessDH0 guessDS0];
fds_Hi = Fpred(thermvals,T_hiCmelt+273.15,Wa_tot_hi_est,Cr_tot_hi_est);
fds_Lo = Fpred(thermvals,T_loCmelt+273.15,Wa_tot_lo_est,Cr_tot_lo_est);
ssmin_Hi = find(fds_Hi < 0.03,1,'first');
ssmin_Lo = find(fds_Lo < 0.03,1,'first');
Tmax_interp_hi_begin = min(length(T_hiCmelt)-Tmax_interp_min,ssmin_Hi);
Tmax_interp_lo_begin = min(length(T_loCmelt)-Tmax_interp_min,ssmin_Lo);
% Pull out the temperatures where DNA is almost all single stranded
%T_HiC_near_Tmax = find (T_hiCmelt > (Tmax - Tmax_interpdist) );
T_HiC_near_Tmax = T_hiCmelt(Tmax_interp_hi_begin:end);
T_LoC_near_Tmax = T_loCmelt(Tmax_interp_lo_begin:end);
% These are indices in the ss1 ref arrays that match the actual temps of
% the melt arrays
try
    T_ss1ref_near_Tmax_Hi = find (T_ss1ref > T_hiCmelt(Tmax_interp_hi_begin));
catch
    if(isempty(Tmax_interp_hi_begin))
        disp('Tmax_interp_hi_begin not found.')
    else
        disp(['Tmax_interp_hi_begin = ',num2str(Tmax_interp_hi_begin)])
    end
    disp('Something is probably wrong with the input data file. See plot')
    figure
    plot(T_hiCmelt,Abs_hiCmelt,'o')
    return
end
T_ss2ref_near_Tmax_Hi = find (T_ss2ref > T_hiCmelt(Tmax_interp_hi_begin));
T_ss1ref_near_Tmax_Lo = find (T_ss1ref > T_hiCmelt(Tmax_interp_lo_begin));
T_ss2ref_near_Tmax_Lo = find (T_ss2ref > T_hiCmelt(Tmax_interp_lo_begin));
%
% Interpolate to a common set of temperatures
try
    Abs_HiCmelt_Tmax = Abs_hiCmelt(Tmax_interp_hi_begin:end);
catch
    if(isempty(Tmax_interp_hi_begin))
        disp('Tmax_interp_hi_begin not found.')
    else
        disp(['Tmax_interp_hi_begin = ',num2str(Tmax_interp_hi_begin)])
    end
    disp('Something is probably wrong with the input data file. See plot')
    figure
    plot(T_hiCmelt,Abs_hiCmelt,'o')
    return
end  
try
Abs_ss1ref_Tmax_Hi = interp1(T_ss1ref(T_ss1ref_near_Tmax_Hi),...
    Abs_ss1ref(T_ss1ref_near_Tmax_Hi),T_HiC_near_Tmax,'pchip');
catch
    if(isempty(T_ss1ref_near_Tmax_Hi))
        disp('T_ss1ref_near_Tmax_Hi not found.')
    else
        disp(['T_ss1ref_near_Tmax_Hi = ',num2str(T_ss1ref_near_Tmax_Hi)])
    end
    disp('Something is probably wrong with the input data file. See plot')
    figure
    plot(T_hiCmelt,Abs_hiCmelt,'o')
    return
end
Abs_ss2ref_Tmax_Hi = interp1(T_ss2ref(T_ss2ref_near_Tmax_Hi),...
    Abs_ss2ref(T_ss2ref_near_Tmax_Hi),T_HiC_near_Tmax,'pchip');
Abs_LoCmelt_Tmax = Abs_loCmelt(Tmax_interp_lo_begin:end);
Abs_ss1ref_Tmax_Lo = interp1(T_ss1ref(T_ss1ref_near_Tmax_Lo),...
    Abs_ss1ref(T_ss1ref_near_Tmax_Lo),T_LoC_near_Tmax,'pchip');
Abs_ss2ref_Tmax_Lo = interp1(T_ss2ref(T_ss2ref_near_Tmax_Lo),...
    Abs_ss2ref(T_ss2ref_near_Tmax_Lo),T_LoC_near_Tmax,'pchip');
% These are vectors covering the high temp range
Theo_Abs_HiCmelt_Tmax = V_ss1_Hi/V_ss1_ssref1*Abs_ss1ref_Tmax_Hi + ...
    V_ss2_Hi/V_ss2_ssref2*Abs_ss2ref_Tmax_Hi;
Theo_Abs_LoCmelt_Tmax = V_ss1_Lo/V_ss1_ssref1*Abs_ss1ref_Tmax_Lo + ...
    V_ss2_Lo/V_ss2_ssref2*Abs_ss2ref_Tmax_Lo;
% Average ratios over the range
Ratio_Hi = mean(Theo_Abs_HiCmelt_Tmax./Abs_HiCmelt_Tmax);
Ratio_Lo = mean(Theo_Abs_LoCmelt_Tmax./Abs_LoCmelt_Tmax);
%
Equal_concs_Hi = mean(Abs_HiCmelt_Tmax./(Abs_ss1ref_Tmax_Hi/Conc1_ref + ...
    Abs_ss2ref_Tmax_Hi/Conc2_ref));
Equal_concs_Lo = mean(Abs_LoCmelt_Tmax./(Abs_ss1ref_Tmax_Lo/Conc1_ref + ...
    Abs_ss2ref_Tmax_Lo/Conc2_ref));
%
% Rescale according to observed values to make fit exact if possible.
% This calculation is not sufficient -- real concentrations of C1 and C2
% are not always proportional to each other.  In other words, the error
% can be different in pipetting into the individual melts.  Below we refine
% estimates based on fitting to a consistent extent of hypochromicity.
disp('Rescaled concentrations:')
C1_tot_Himelt = (V_ss1_Hi/V_ss1_ssref1)*(Conc1_ref/Ratio_Hi);
C2_tot_Himelt = (V_ss2_Hi/V_ss2_ssref2)*(Conc2_ref/Ratio_Hi);
C1_tot_Lomelt = (V_ss1_Lo/V_ss1_ssref1)*(Conc1_ref/Ratio_Lo);
C2_tot_Lomelt = (V_ss2_Lo/V_ss2_ssref2)*(Conc2_ref/Ratio_Lo);
disp(['C1_tot_Himelt = ', num2str(C1_tot_Himelt),'C2_tot_Himelt = ', num2str(C2_tot_Himelt),...
    'C1_tot_Lomelt = ', num2str(C1_tot_Himelt),'C2_tot_Lomelt = ', num2str(C2_tot_Himelt)]);
disp(['Actual/theo HiC abs at Tmax is ',num2str(Ratio_Hi),' and ',...
        'actual/theo LoC abs at Tmax is ',num2str(Ratio_Lo)])
%    
if ((abs(Ratio_Hi - 1) > Ratio_tolerance) || (abs(Ratio_Lo - 1) > Ratio_tolerance))
    % Flag possible pipetting error
    disp(['QC Caution: actual/theo abs ratio differs from ideal by more than ',...
        'ratio tolerance ',num2str(Ratio_tolerance)])
end

% We could just combine it in here rather than splitting out the routine -- we are
% specializing anyway! OTOH I might want to plug in a different fitting
% routine at some point again
% Fitvals is a structure.
% 
% Do all graphing inside the function call -- otherwise management becomes
% very hard.

[Fitvals,~] = AbsfitwE_cv_2018(filestem,...
     C1_tot_Himelt_init,C2_tot_Himelt_init,C1_tot_Lomelt_init,C2_tot_Lomelt_init,...
     C1_tot_Himelt,C2_tot_Himelt,C1_tot_Lomelt,C2_tot_Lomelt,...
     RefT,Conc1_ref,Conc2_ref,T_hiCmelt,Abs_hiCmelt,T_loCmelt,Abs_loCmelt,...
     T_ss1ref,T_ss2ref,Abs_ss1ref,Epsil_ss1_RefT,Abs_ss2ref,Epsil_ss2_RefT,...
     Equal_concs_Hi,Equal_concs_Lo,...
     debug,scale_factor,scale_flag,Tlow,Tmin,Tmax,Thigh);
% function Fitvals = AbsfitwE_cv_2018(filestem,C1_tot_Himelt,C2_tot_Himelt,...
%    C1_tot_Lomelt,C2_tot_Lomelt,RefT,...
%     Conc1_ref,Conc2_ref,T_hiCmelt,Abs_hiCmelt,T_loCmelt,Abs_loCmelt,T_ss1ref,...
%       Abs_ss1ref,T_ss2ref,Abs_ss2ref,...
%     Equal_concs_Hi,Equal_concs_Lo,debug,scale_factor)
% A    




% format short e;
% if Refspec
%     disp('Found Refspec, running AbsfitwE1_2016')
%     Fitvals = AbsfitwE1_2016(filestem,CT,T,Abso,Tref,Aref,ExW,inum);
% elseif Nobaseline
%     disp('Nobaseline, running AbsfitwE3')
%     Fitvals = AbsfitwE3(filestem,CT,T,Abso,Refhypo,ExW,EpsW,EpsC,inum);
% else
%     disp('Default, running AbsfitwE2_2016')
%     Fitvals = AbsfitwE2_2016(filestem,CT,T,Abso,ExW,EpsW,EpsC,inum);
% end    
%Print file name (eps format) is Fitvals.Name
%Fitvals.Nums = [thermvals(1) ci(1,:) thermvals(2) ci(2,:) DG37 DG37max DG37min ...
%        TM TMmax TMmin mds mss Amin Amax resnorm flag];
% Flags and handle to plot graphic is in Fitvals.Outflags
%
% Output 2-state results:
%
%fieldnames(Fitvals)
% Need to change answrite function below to change which numbers are
% output. Also need to change the output of the fitting routine to make
% Fitvals concordant

% answrite(Fitvals.Nums);
% if (Fitvals.Wa_is_ss1_Hi)
%     disp('Wa_Hi is Single Strand 1')
% else
%     disp('Wa_Hi is Single Strand 2')
% end
% if (Fitvals.Wa_is_ss1_Lo)
%     disp('Wa_Lo is Single Strand 1')
% else
%     disp('Wa_Lo is Single Strand 2')
% end
% disp(['Tmin, Tamx = ',num2str(Tmin),' ',num2str(Tmax)])
% disp(['Concentrations iterated ',num2str(Fitvals.Iter),' times'])
% disp(['Graph printout file is ',Fitvals.Name])
% Uncomment four lines below to output matout file
%savefile2 = [filestem '.matout'];
%savevals = [Fitvals.Nums];
%save(savefile2,'savevals','-ASCII','-tabs');
%disp(['2-state save file is ',savefile2])
%
% if Nobaseline
%     savevals'
%     return
% end
%Flag1 = Fitvals.Outflags(1);
%Flag2 = Fitvals.Outflags(2);
%
end
% 
% function answrite(pvals)
% disp(' ')
% disp(['Delta H0  : ',num2str(pvals(1:3))])
% disp(['Delta S0  : ',num2str(pvals(4:6))])
% disp(['Delta G037: ',num2str(pvals(7:9))])
% disp(['TM_Hi     : ',num2str(pvals(10:12))])
% disp(['TM_Lo     : ',num2str(pvals(13:15))])
% disp(['Wa, Cr Hi : ',num2str(pvals(16:17))])
% disp(['Wa, Cr Lo : ',num2str(pvals(18:19))])
% disp(['mds,EpDS  : ',num2str(pvals(20:21))])
% disp(['Resnorm   : ',num2str(pvals(22))])
% disp(' ')
% end
