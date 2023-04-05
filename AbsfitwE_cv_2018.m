function [Fitvals,Pfiles] = AbsfitwE_cv_2018(filestem,C1_tot_Himelt_init,...
    C2_tot_Himelt_init,C1_tot_Lomelt_init,C2_tot_Lomelt_init,...
    C1_tot_Himelt,C2_tot_Himelt,C1_tot_Lomelt,C2_tot_Lomelt,RefT,...
    Conc1_ref,Conc2_ref,T_hiCmelt,Abs_hiCmelt,T_loCmelt,Abs_loCmelt,...
    T_ss1ref,T_ss2ref,Abs_ss1ref,Epsil_ss1_RefT,Abs_ss2ref,Epsil_ss2_RefT,...
    Equal_concs_Hi,Equal_concs_Lo,debug,scale_factor,scale_flag,Tlow,Tmin,Tmax,Thigh)

format ShortE
% All arguments required. This should never be called directly by the user.

% Fits absorbance vs. temperature data to a two-state melting model.
% according to math at http://www.biochem.umd.edu/biochem/kahn/teach_res/hybridization.html
% with significant extensions due to presence of excess Watson strand
% Subtracts off a scaled supplied reference A vs. T curve for the excess ss
% and therefore does not calculate a linear ssDNA baseline because it has the real data
% See the file:
% /Users/jdk/Documents/LNA\ project/Dealing\ with\ different\ si.doc
% Input total strand concentrations, absorbance vs temperature arrays
% 

%Calculate Absor(Tmin), Absor(Tmax), Baseline slopes
% Then fit calculated A vs. T and dA/dT vs. T
%   to delta H and delta S
% Generate graphs of A vs. T, dA/dT vs. T, residual vs. T
% Output report, eventually export best-fit theoretical numbers back to Excel 
% Report should contain sequence, graphs, table of numbers (not done yet)
% Need to add in all the three-state stuff
%
%debug = 1 ; % Set to 1 for debugging, 0 for production (now set in
%calling routine)
% Define the Fitvals structure
% Seems to need dummy names for plot files?
Fitvals = struct('Params',[],'Resnorms',[],'Conc_errors',[]);
Pfiles = {'f1','f2','f3','f4','f5','f6'};
if (debug >= 1)
    debf = figure;  % create figure and handle for updating progress
    debf.Units = 'inches';
    debf.Position = [7 5 12 6];
end
fig_handle_mat = zeros(6,1);
Hi_color = [0 0.5 0];
Lo_color = [0.8 0.5 0];
%one_last = 0;
%flag = [1 1];   % The first is for short baseline problems,
% the second for negative baseline slopes.
% Any number other than 1 indicates a problem.
% Set nineplots to "false" here if we want to shut off automatic 3-state
% fitting
%nineplots = true;
dslim = 8;
R = 1.987;
% Bizarre dependence on Conc_tol suggests something is wrong
% Conc_tol = 0.001;
% switched_Hi = 0;
% switched_Lo = 0;
%numtries = 100;
% Interpolate ssDNA reference spectra to match the temperatures of the melt
% spectra. This takes care of pruning the reference data to the Tmin < T < Tmax range.
% At this point we are also defining the Watson strand as the one in
% excess.
% We need to allow for Watson and Crick to be different strands in the two
% melts!
% Need to allow concentrations of the two strands to vary


param_matrix = zeros(11,18);
resarray = zeros(1,6);
conc_err_matrix = zeros(4,6);
% Limit for flagging concentration problems - ~ avg of 10%
Conc_err_concern = 0.4;
%
L = length(Abs_hiCmelt);
Abs_der_Hi = zeros(1,L-2);
for J=2:L-1
    Abs_der_Hi(J-1)  = (Abs_hiCmelt(J-1)-Abs_hiCmelt(J+1))/(T_hiCmelt(J-1) - T_hiCmelt(J+1));
end
L = length(Abs_loCmelt);
Abs_der_Lo = zeros(1,L-2);
for J=2:L-1
    Abs_der_Lo(J-1)  = (Abs_loCmelt(J-1)-Abs_loCmelt(J+1))/(T_loCmelt(J-1) - T_loCmelt(J+1));
end


[maxslope,maxindex] = max(Abs_der_Hi);
guesstm = T_hiCmelt(maxindex+1);
%Tm_Hi = guesstm;
guessDH0 = -6*R*(maxslope/(max(Abs_hiCmelt)-min(Abs_hiCmelt)))*(guesstm+273.15)^2;
if (debug > 1)
    disp(['AbsfitwE_2018 Initial guess for Tm_Hi is ',num2str(guesstm),' C'])
    disp(['Initial guess for enthalpy is ',num2str(guessDH0),' cal/mole'])
end    

%scale_factor = mean(Abs_hiCmelt)/mean(Abs_loCmelt);
% To look at only the high CT melt scale_factor = 0.001;
% To look only at low CT, scale_factor = 1000;
% This is the beginning of the big fitting loop

%resarray = zeros(numtries,1);
%deltarray = zeros(numtries,2);

% Do all these interpolations only once and then can simply pick among them.
Abs_ss1_ref_interp_Hi = interp1(T_ss1ref,Abs_ss1ref,T_hiCmelt,'pchip');
Abs_ss2_ref_interp_Hi = interp1(T_ss2ref,Abs_ss2ref,T_hiCmelt,'pchip');
Abs_ss1_ref_interp_Lo = interp1(T_ss1ref,Abs_ss1ref,T_loCmelt,'pchip');
Abs_ss2_ref_interp_Lo = interp1(T_ss2ref,Abs_ss2ref,T_loCmelt,'pchip');

% This is the loop where we do the calculations

for iter = 1:6
    % for the first two iterations (iter == 1 or 2) determine Wa and Cr from input
    % concentrations and do not let the concentrations vary
    if (iter == 1)
        if (C1_tot_Himelt_init >= C2_tot_Himelt_init)
    % Strand #1 is the initial Watson strand (in excess)
            Wa_is_ss1_Hi = true;
%    Wa_is_ss2_Hi = false;
            Abs_Wa_ref_interp_Hi = Abs_ss1_ref_interp_Hi;
            Abs_Cr_ref_interp_Hi = Abs_ss2_ref_interp_Hi;
            Wa_tot_Hi = C1_tot_Himelt_init;
            Cr_tot_Hi = C2_tot_Himelt_init;
            Wa_ref_ss_Hi = Conc1_ref;
            Cr_ref_ss_Hi = Conc2_ref;
            if ((debug >= 2) && (iter == 1))
                figure
                hold on
                plot(T_ss1ref,Abs_ss1ref,'o',T_hiCmelt,Abs_Wa_ref_interp_Hi)
                plot(T_ss2ref,Abs_ss2ref,'o',T_hiCmelt,Abs_Cr_ref_interp_Hi)
            end
        else
    % Strand #2 is the Watson strand (in excess)
            Wa_is_ss1_Hi = false;
%    Wa_is_ss2_Hi = true;
            Abs_Wa_ref_interp_Hi = Abs_ss2_ref_interp_Hi;
            Abs_Cr_ref_interp_Hi = Abs_ss1_ref_interp_Hi;
            Wa_tot_Hi = C2_tot_Himelt_init;
            Cr_tot_Hi = C1_tot_Himelt_init;
            Wa_ref_ss_Hi = Conc2_ref;
            Cr_ref_ss_Hi = Conc1_ref;
            if (debug >= 2 && iter == 1)
                figure
                hold on
                plot(T_ss2ref,Abs_ss2ref,'o',T_hiCmelt,Abs_Wa_ref_interp_Hi)
                plot(T_ss1ref,Abs_ss1ref,'o',T_hiCmelt,Abs_Cr_ref_interp_Hi)
            end
        end 
        % Allow W and C to be independent for the two melts even though
        % usually they will be the same.
        if (C1_tot_Lomelt_init >= C2_tot_Lomelt_init)
            Wa_is_ss1_Lo = true;
%    Wa_is_ss2_Lo = false;
            Abs_Wa_ref_interp_Lo = Abs_ss1_ref_interp_Lo;
            Abs_Cr_ref_interp_Lo = Abs_ss2_ref_interp_Lo;
            Wa_tot_Lo = C1_tot_Lomelt_init;
            Cr_tot_Lo = C2_tot_Lomelt_init;
            Wa_ref_ss_Lo = Conc1_ref;
            Cr_ref_ss_Lo = Conc2_ref;
        else
            Wa_is_ss1_Lo = false;
%     Wa_is_ss2_Lo = true;
            Abs_Wa_ref_interp_Lo = Abs_ss2_ref_interp_Lo;
            Abs_Cr_ref_interp_Lo = Abs_ss1_ref_interp_Lo;
            Wa_tot_Lo = C2_tot_Lomelt_init;
            Cr_tot_Lo = C1_tot_Lomelt_init;
            Wa_ref_ss_Lo = Conc2_ref;
            Cr_ref_ss_Lo = Conc1_ref;
        end
        % Save these initial apparent values
         Wa_tot_Hi_init = Wa_tot_Hi;
         Cr_tot_Hi_init = Cr_tot_Hi;
         Wa_tot_Lo_init = Wa_tot_Lo;
         Cr_tot_Lo_init = Cr_tot_Lo;
    end
    
    
    if (iter == 2)
        if (C1_tot_Himelt >= C2_tot_Himelt)
    % Strand #1 is the initial Watson strand (in excess)
            Wa_is_ss1_Hi = true;
%    Wa_is_ss2_Hi = false;
            Abs_Wa_ref_interp_Hi = Abs_ss1_ref_interp_Hi;
            Abs_Cr_ref_interp_Hi = Abs_ss2_ref_interp_Hi;
            Wa_tot_Hi = C1_tot_Himelt;
            Cr_tot_Hi = C2_tot_Himelt;
            Wa_ref_ss_Hi = Conc1_ref;
            Cr_ref_ss_Hi = Conc2_ref;
            if ((debug >= 2) && (iter == 1))
                figure
                hold on
                plot(T_ss1ref,Abs_ss1ref,'o',T_hiCmelt,Abs_Wa_ref_interp_Hi)
                plot(T_ss2ref,Abs_ss2ref,'o',T_hiCmelt,Abs_Cr_ref_interp_Hi)
            end
        else
    % Strand #2 is the Watson strand (in excess)
            Wa_is_ss1_Hi = false;
%    Wa_is_ss2_Hi = true;
            Abs_Wa_ref_interp_Hi = Abs_ss2_ref_interp_Hi;
            Abs_Cr_ref_interp_Hi = Abs_ss1_ref_interp_Hi;
            Wa_tot_Hi = C2_tot_Himelt;
            Cr_tot_Hi = C1_tot_Himelt;
            Wa_ref_ss_Hi = Conc2_ref;
            Cr_ref_ss_Hi = Conc1_ref;
            if (debug >= 2 && iter == 2)
                figure
                hold on
                plot(T_ss2ref,Abs_ss2ref,'o',T_hiCmelt,Abs_Wa_ref_interp_Hi)
                plot(T_ss1ref,Abs_ss1ref,'o',T_hiCmelt,Abs_Cr_ref_interp_Hi)
            end
        end 
        % Allow W and C to be independent for the two melts even though
        % usually they will be the same.
        if (C1_tot_Lomelt >= C2_tot_Lomelt)
            Wa_is_ss1_Lo = true;
%    Wa_is_ss2_Lo = false;
            Abs_Wa_ref_interp_Lo = Abs_ss1_ref_interp_Lo;
            Abs_Cr_ref_interp_Lo = Abs_ss2_ref_interp_Lo;
            Wa_tot_Lo = C1_tot_Lomelt;
            Cr_tot_Lo = C2_tot_Lomelt;
            Wa_ref_ss_Lo = Conc1_ref;
            Cr_ref_ss_Lo = Conc2_ref;
        else
            Wa_is_ss1_Lo = false;
%     Wa_is_ss2_Lo = true;
            Abs_Wa_ref_interp_Lo = Abs_ss2_ref_interp_Lo;
            Abs_Cr_ref_interp_Lo = Abs_ss1_ref_interp_Lo;
            Wa_tot_Lo = C2_tot_Lomelt;
            Cr_tot_Lo = C1_tot_Lomelt;
            Wa_ref_ss_Lo = Conc2_ref;
            Cr_ref_ss_Lo = Conc1_ref;
        end
        % Save these initial apparent values
         Wa_tot_Hi_init = Wa_tot_Hi;
         Cr_tot_Hi_init = Cr_tot_Hi;
         Wa_tot_Lo_init = Wa_tot_Lo;
         Cr_tot_Lo_init = Cr_tot_Lo;
    end
% In later iterations we are optimizing concentrations making all four
% possible assumptions about who is Wa and who is Cr
% we set the initial concentrations equal to the concentration determined in the calling
% routine that would correspond to the correct total Abs, and then and allow Wa to go up
% and Cr to go down
    if (iter > 2)
        Wa_tot_Hi = Equal_concs_Hi;
        Cr_tot_Hi = Equal_concs_Hi;
        Wa_tot_Lo = Equal_concs_Lo;
        Cr_tot_Lo = Equal_concs_Lo;
    end
% For iter == 2 we are assuming that Wa = ss1 for Hi and Lo melts
% we set the initial concentrations equal and allow Wa to go up and Cr to
% go down
    if (iter == 3)  
        Wa_is_ss1_Hi = true;
%    Wa_is_ss2_Hi = false;
        Abs_Wa_ref_interp_Hi = Abs_ss1_ref_interp_Hi;
        Abs_Cr_ref_interp_Hi = Abs_ss2_ref_interp_Hi;
        Wa_ref_ss_Hi = Conc1_ref;
        Cr_ref_ss_Hi = Conc2_ref;
        Wa_is_ss1_Lo = true;
%    Wa_is_ss2_Lo = false;
        Abs_Wa_ref_interp_Lo = Abs_ss1_ref_interp_Lo;
        Abs_Cr_ref_interp_Lo = Abs_ss2_ref_interp_Lo;
        Wa_ref_ss_Lo = Conc1_ref;
        Cr_ref_ss_Lo = Conc2_ref;
    end
% iter == 3, we are optimizing concentrations assuming that Wa = ss1 for Hi melt and
% Wa = ss2 for Lo melt
    if (iter == 4)  
        Wa_is_ss1_Hi = true;
        Abs_Wa_ref_interp_Hi = Abs_ss1_ref_interp_Hi;
        Abs_Cr_ref_interp_Hi = Abs_ss2_ref_interp_Hi;
        Wa_ref_ss_Hi = Conc1_ref;
        Cr_ref_ss_Hi = Conc2_ref;
        Wa_is_ss1_Lo = false;
        Abs_Wa_ref_interp_Lo = Abs_ss2_ref_interp_Lo;
        Abs_Cr_ref_interp_Lo = Abs_ss1_ref_interp_Lo;
        Wa_ref_ss_Lo = Conc2_ref;
        Cr_ref_ss_Lo = Conc1_ref;
    end
% iter == 4 we are optimizing concentrations assuming that Wa = ss2 for Hi melt and
% Wa = ss1 for Lo melt
% we set the initial concentrations equal and allow Wa to go up and Cr to
% go down
    if (iter == 5)  
        Wa_is_ss1_Hi = false;
        Abs_Wa_ref_interp_Hi = Abs_ss2_ref_interp_Hi;
        Abs_Cr_ref_interp_Hi = Abs_ss1_ref_interp_Hi;
        Wa_ref_ss_Hi = Conc2_ref;
        Cr_ref_ss_Hi = Conc1_ref;
        Wa_is_ss1_Lo = true;
        Abs_Wa_ref_interp_Lo = Abs_ss1_ref_interp_Lo;
        Abs_Cr_ref_interp_Lo = Abs_ss2_ref_interp_Lo;
        Wa_ref_ss_Lo = Conc1_ref;
        Cr_ref_ss_Lo = Conc2_ref;
    end
% iter == 5 we are optimizing concentrations assuming that Wa = ss2 for Hi melt and
% Wa = ss1 for Lo melt
% we set the initial concentrations equal and allow Wa to go up and Cr to
% go down
    if (iter == 6)  
        Wa_is_ss1_Hi = false;
        Abs_Wa_ref_interp_Hi = Abs_ss2_ref_interp_Hi;
        Abs_Cr_ref_interp_Hi = Abs_ss1_ref_interp_Hi;
        Wa_ref_ss_Hi = Conc2_ref;
        Cr_ref_ss_Hi = Conc1_ref;
    %
        Wa_is_ss1_Lo = false;
        Abs_Wa_ref_interp_Lo = Abs_ss2_ref_interp_Lo;
        Abs_Cr_ref_interp_Lo = Abs_ss1_ref_interp_Lo;
        Wa_ref_ss_Lo = Conc2_ref;
        Cr_ref_ss_Lo = Conc1_ref;
    end

% Now back to steps done for each iteration
    
% Calculate limits (done in calling function now)
% At high T, Wa*epsilWa + Cr*epsilCr = observed high T absorbance
% If Wa = Cr we have Wa = Cr = observed high T absorbance/(epsilWa + epsilCr)

%Concatenate vectors for global fitting
% At this point we scale the concatenated vectors so the hi conc melt does
% not dominate the fit. We need to make sure not to use the _all vectors
% for any actual regeneration of real fits.
% [Now redone to allow passing a cell array of melts to Apred]
% Watch the semicolon to concatenate column vectors rather than make matrices!
%scale_factor = 1;
%scaling_vector = [ones(size(Abs_Wa_ref_interp_Hi)) ; ...
%    scale_factor*ones(size(Abs_Wa_ref_interp_Lo))];
%TK = [T_hiCmelt ; T_loCmelt] + 273.15; 
%Abs_expt_all = [Abs_hiCmelt ; Abs_loCmelt];

% Need to recalculate these each time. Maybe we could save a few milliseconds
% by precalculating and reassigning.
    Epsil_Wa_vec_Hi = Abs_Wa_ref_interp_Hi/Wa_ref_ss_Hi;
    Epsil_Wa_vec_Lo = Abs_Wa_ref_interp_Lo./Wa_ref_ss_Lo;  
    Epsil_Cr_vec_Hi = Abs_Cr_ref_interp_Hi./Cr_ref_ss_Hi;
    Epsil_Cr_vec_Lo = Abs_Cr_ref_interp_Lo./Cr_ref_ss_Lo;

%Epsil_Wa_vec_all = [Epsil_Wa_vec_Hi ; Epsil_Wa_vec_Lo];
%Epsil_Cr_vec_all = [Epsil_Cr_vec_Hi ; Epsil_Cr_vec_Lo];
%
% These vectors of concentrations are unnecessary if we are passing a cell
% array?
%
%Wa_tot_vec_Hi = Wa_tot_Hi*ones(size(Abs_Wa_ref_interp_Hi));
%Wa_tot_vec_Lo = Wa_tot_Lo*ones(size(Abs_Wa_ref_interp_Lo));
%Wa_tot_vec_all = [Wa_tot_vec_Hi ; Wa_tot_vec_Lo];
%Cr_tot_vec_Hi = Cr_tot_Hi*ones(size(Abs_Cr_ref_interp_Hi));
%Cr_tot_vec_Lo = Cr_tot_Lo*ones(size(Abs_Cr_ref_interp_Lo));
%Cr_tot_vec_all = [Cr_tot_vec_Hi ; Cr_tot_vec_Lo];
% OR... build cell arrays instead
    meltca = { 1, scale_factor ; ...
        Wa_tot_Hi , Wa_tot_Lo ; ...
        Cr_tot_Hi , Cr_tot_Lo ; ...
        T_hiCmelt+273.15 , T_loCmelt+273.15 ; ...
        Abs_hiCmelt , Abs_loCmelt ; ...
        Epsil_Wa_vec_Hi , Epsil_Wa_vec_Lo ;...
        Epsil_Cr_vec_Hi , Epsil_Cr_vec_Lo};
    
% Establish initial estimates in first iteration  and independently in the second iteration  
    if (iter == 1 || iter == 2)
        Abs_Excess_ssDNA_hiCmelt = (Wa_tot_Hi - Cr_tot_Hi)*Epsil_Wa_vec_Hi;
        Abs_pairedDNA_hiCmelt = Abs_hiCmelt - Abs_Excess_ssDNA_hiCmelt;
        Abs_Excess_ssDNA_loCmelt = (Wa_tot_Lo - Cr_tot_Lo)*Epsil_Wa_vec_Lo;
        Abs_pairedDNA_loCmelt = Abs_loCmelt - Abs_Excess_ssDNA_loCmelt;

% Enthalpy was guessed at above or in previous iteration
%if (iter == 1)
        guessDS0 = (guessDH0/(guesstm+273.15)) - R*log(Wa_tot_Hi - Cr_tot_Hi/2);
        disp(['Initial guess for entropy is ',num2str(guessDS0),' cal/moleK'])
% else
%     guessDH0 = DH0;
%     guessDS0 = (DH0/(Tm_Hi+273.15)) - R*log(Wa_tot_Hi - Cr_tot_Hi/2);
%     if (debug >= 1)
%         disp(['Next guess for entropy is ',num2str(guessDS0),' cal/moleK'])
%     end
% end
        thermvals = [guessDH0 guessDS0];
% The number of evals and iterations is much larger than is needed unless
% mss and mds are being optimized as well as DH and DS
%mfe = optimget(options,'MaxFunEvals');
%thermvals = lsqcurvefit(@Apred,guess,TK,Absor,[],[],options,CT,ExW,mss,mds,Amin,Amax)
% Fit to derivative data
% Run Iter from 1:2 for just DH and DS fits, with the second one being
% adaptive limits on the all-ds and all-ss regimes. Run Iter from 1:3 to
% allow unrestrained fitting of mss and mds
%for Iter = 1:inum;
% First set reasonable conservative limits for ds, ss
% 	if Iter == 1 ;
%         dsmax = 16;
%         ssmin = L-20;
%         mds = 2*mean (Ader(1:dsmax-1));
%         mss = 2*mean (Ader(ssmin-1:L-2));
%         Amin = mean(Absor(1:5))-mds*(T(3) - T(1));
%         Amax = mean(Absor(L-5:L))+mss*(T(L)-T(L-2));
% 	end
%
%   The code here uses the estimates for DH and DS to calculate fds so that
%   we can use all appropriate data to get the best estimates for epsilon
%   dsDNA and the slope of epsilon
        fds = Fpred(thermvals,T_hiCmelt+273.15,Wa_tot_Hi,Cr_tot_Hi);
        dsmax1 = find(fds > 0.99,1,'last');
%ssmin = find(fds < 0.01,1);
%guess = thermvals;
        if (isempty(dsmax1) || dsmax1 < dslim)
            disp('ds baseline too short for Hi melt');
%            flag(1) = 2;
% Build in best guess for baseline slope here, just a short linear
% extrapolation backward (5 temperature points)
            p1 = polyfit(T_hiCmelt(1:5),Abs_pairedDNA_hiCmelt(1:5)/Cr_tot_Hi,1);
        else
% To estimate the slope of dsDNA absorbance, we need to subtract off the
% absorbance of the excess Watson
            p1 = polyfit(T_hiCmelt(1:dsmax1),Abs_pairedDNA_hiCmelt(1:dsmax1)/Cr_tot_Hi,1);
        end
        disp(['mds est from hi CT melt is ',num2str(p1(1))]);
        % Do the same thing for the lo concentration melt.
        fds = Fpred(thermvals,T_loCmelt+273.15,Wa_tot_Lo,Cr_tot_Lo);
        dsmax2 = find(fds > 0.99,1,'last');
        disp (['dsmax1 = ',num2str(dsmax1),' and dsmax2 = ',num2str(dsmax2)])
%ssmin = find(fds < 0.01,1);
%guess = thermvals;
        if (isempty(dsmax2) || dsmax2 < dslim)
            disp('ds baseline too short for Lo melt');
            disp('using mds and Abs from Hi melt');
%            flag(1) = 2;
% If the baseline is bad, just use the estimate from the higher
% concentration
            p2 = p1;
        else
            disp (['dsmax2 = ',num2str(dsmax2),' and dslim = ',num2str(dslim)])
            p2 = polyfit(T_loCmelt(1:dsmax2),Abs_pairedDNA_loCmelt(1:dsmax2)/Cr_tot_Lo,1);
        end
        disp(['mds est from lo CT melt is ',num2str(p2(1))]);
% if (debug >= 1) % && ((iter <= 5) || mod(iter,5) == 0 || one_last == 1))
%     figure(debf)
%     subplot(2,4,1)
%     hold on
%     box on
%     plot(T_hiCmelt(1:dsmax1+10),Abs_pairedDNA_hiCmelt(1:dsmax1+10)/Cr_tot_Hi,'o')
%     %plot(T_hiCmelt(1:dsmax1),Abs_pairedDNA_hiCmelt(1:dsmax1)/Cr_tot_Hi)
%     plot(T_hiCmelt(1:dsmax1+10),polyval(p1,T_hiCmelt(1:dsmax1+10)))
%     plot(T_loCmelt(1:dsmax2),Abs_pairedDNA_loCmelt(1:dsmax2)/Cr_tot_Lo)
%     plot(T_loCmelt(1:dsmax2),polyval(p2,T_loCmelt(1:dsmax2)))
%     title('DS baseline fitting')
% end
%   Average the values from the hi and low conc melts
        mds_est = 0.5*(p1(1)+p2(1));
        Epsil_WC_RefT_est = 0.5*(p1(2)+p2(2));
        disp(['mds est from average of melts is ',num2str(mds_est)]);
        disp(['Epsil dsDNA at RefT est from average of melts is ',...
            num2str(Epsil_WC_RefT_est)]);
%flag(2) = 1;        % reset the mds/mss flag because we are about to check again
        if (mds_est < 0 && mds_est*(min(T_hiCmelt)-max(T_hiCmelt)) > ...
                0.05*(max(Abs_hiCmelt) - min(Abs_hiCmelt)))
            disp(['Caution: Negative mds too large: ',num2str(mds_est)]);
%            flag(2) = 2;
        end
% To implement: sometimes fitting fails. If it sibecause of baselines, we
% could....
% We do know something about the absorbance of dsDNA. We will assume that hypochromicity
% is between 2 and 40 percent -- usually it is about 15 %. Furthermore the
% slope mds of the dsDNA extinction coefficient can be no more than a
% change of 50% over 100 °C. 
        
        
        
% if (mss < 0 && mss*(min(T)-max(T)) > 0.05*(Amax-Amin))
% 	disp([ 'Negative mss too large ',num2str(mss)]);
% 	flag(2) = 3;
% end
% This line is the guts of the whole script. It does the actual fitting
% disp(' ');
        firstfitvals = [thermvals(1), thermvals(2), mds_est, Epsil_WC_RefT_est];
        disp(['Input fitvals: ',num2str(firstfitvals)]);
%class(TK)
%Troubleshooting
% if (debug > 1)
% %    TK
% size(TK)
% size(Wa_tot_vec_all)
% size(Cr_tot_vec_all)
% size(Epsil_Wa_vec_all)
% size(Epsil_Cr_vec_all)
% size(Abs_expt_all)
% end

%
%
% if (debug > 1)
% class(fitvals)
% class(Wa_tot_vec_all)
% class(Cr_tot_vec_all)
% class(Epsil_Wa_vec_all)
% class(Epsil_Cr_vec_all)
% class(Abs_expt_all)
% class(RefT)
% end
% upper and lower bounds
%min([p1(1),p2(1)])
        options = optimset('TolFun',1e-15,'TolX',1e-15,...
            'MaxFunEvals',25000,'MaxIter',5000);
        lb = [-300000,-1000,min([p1(1),p2(1)]),min([p1(2),p2(2)])];
        ub = [0,0,max([p1(1),p2(1)]),max([p1(2),p2(2)])];
    
% else
%     lb = [-300000,-1000,min([p1(1),p2(1)]),Epsil_WC_RefT];
%     ub = [0,0,max([p1(1),p2(1)]),Epsil_WC_RefT];
% end
% if (debug >= 1)
% else
%    options = optimset('Display','off','TolFun',1e-15,'TolX',...
% 1e-15,'MaxFunEvals',15000,'MaxIter',1000);
% end
% [thermvals,resnorm, derres,exitflag,output,lambda,jacobian] = ...
% lsqnonlin(@Apred,fitvals,lb,ub,options,...
%   TK,Wa_tot_vec_all,Cr_tot_vec_all,Epsil_Wa_vec_all,Epsil_Cr_vec_all,...
%Abs_expt_all,scaling_vector,RefT);
% suppress error msg to keep identity of unused exitflag and lambda args
%
% What would it looklike if Apred took cell arrays (rather than matrices
% because data sets can be different lengths)
        [thermvals,resnorm, derres,~,output,~,jacobian] = ...
            lsqnonlin(@Apred_ca,firstfitvals,lb,ub,options,meltca,RefT);

        disp([ 'Iteration #1 Resnorm with refined baseline lims is ', num2str(resnorm),...
            ' after ',num2str(output.iterations),' lsqnonlin iterations' ]);
        disp(' ');
%
%---------------------
% nlparci provides 95% confidence limits derived from the output of
% lsqnonlin
try
        ci = nlparci(thermvals,derres,jacobian);
        Uplim4 = ci(:,1);
	    Dnlim4 = ci(:,2);
catch
    disp(['Nlparci failed for iteration ',num2str(iter),'. Setting params to best fit.'])
    Uplim4 = thermvals;
    Dnlim4 = thermvals;
end
% In reporting out the parameter matrix we want to use a defined order of strands
% because Wa and Cr are relative.
% Report out ss1 Hi, ss2 Hi, ss1 Lo, ss2 Lo
        param_matrix(1:4,(iter-1)*3+1) = thermvals;
        param_matrix(1:4,(iter-1)*3+2) = Uplim4;
        param_matrix(1:4,(iter-1)*3+3) = Dnlim4;
        if (Wa_is_ss1_Hi)
            param_matrix(5:6,(iter-1)*3+1) = [Wa_tot_Hi Cr_tot_Hi];
        else
            param_matrix(5:6,(iter-1)*3+1) = [Cr_tot_Hi Wa_tot_Hi];
        end
        if (Wa_is_ss1_Lo)
            param_matrix(7:8,(iter-1)*3+1) = [Wa_tot_Lo Cr_tot_Lo];
        else
            param_matrix(7:8,(iter-1)*3+1) = [Cr_tot_Lo Wa_tot_Lo];
        end
        param_matrix(5:8,(iter-1)*3+2) = param_matrix(5:8,(iter-1)*3+1);
        param_matrix(5:8,(iter-1)*3+3) = param_matrix(5:8,(iter-1)*3+1);
        % here upper and lower limits are just the initial fixed values
        Up_Wa_tot_Hi = Wa_tot_Hi;
        Up_Wa_tot_Lo = Wa_tot_Lo;
        Up_Cr_tot_Hi = Cr_tot_Hi;
        Up_Cr_tot_Lo = Cr_tot_Lo;
        Dn_Wa_tot_Hi = Wa_tot_Hi;
        Dn_Wa_tot_Lo = Wa_tot_Lo;
        Dn_Cr_tot_Hi = Cr_tot_Hi;
        Dn_Cr_tot_Lo = Cr_tot_Lo;
        resarray(iter) = resnorm;
        %conc_err_matrix (1:4,1) = [0 0 0 0];
        conc_err_matrix (1:4,iter) = [
            (param_matrix(5,(iter-1)*3+1) - C1_tot_Himelt_init)/C1_tot_Himelt_init ...
            (param_matrix(6,(iter-1)*3+1) - C2_tot_Himelt_init)/C2_tot_Himelt_init ...
            (param_matrix(7,(iter-1)*3+1) - C1_tot_Lomelt_init)/C1_tot_Lomelt_init ...
            (param_matrix(8,(iter-1)*3+1) - C2_tot_Lomelt_init)/C2_tot_Lomelt_init];
        verrconc = vecnorm(conc_err_matrix);
        %flag(2) = 1;
        if (thermvals(3) < 0 ) % && mds*(min(T)-max(T)) > 0.05*(Amax-Amin 
            disp([' Iteration ',num2str(iter),': Negative mds = ',num2str(thermvals(3))]);
            %flag(2) = 2;
        end
%Calcnorm = sum((Ader-Apredder).^2);
%Hinorm = sum((Ader-Hider).^2);
%Lonorm = sum((Ader-Loder).^2);
% Calculate derived parameters:
        DH0 = thermvals(1);
        DS0 = thermvals(2);
        mds = thermvals(3);
        Epsil_WC_RefT = thermvals(4);
        DG037 = thermvals(1)-310.15*thermvals(2);
        DG037_min = Uplim4(1)-310.15*Uplim4(2);
        DG037_max = Dnlim4(1)-310.15*Dnlim4(2);
        param_matrix(9,((iter-1)*3+1):((iter-1)*3+3)) = [DG037 DG037_max DG037_min];        
        Tm_Hi = -273.15+thermvals(1)/(thermvals(2)+R*log(Wa_tot_Hi - Cr_tot_Hi/2));
        Tm_Hi_max = -273.15+Uplim4(1)/(Uplim4(2)+R*log(Wa_tot_Hi - Cr_tot_Hi/2));
        Tm_Hi_min = -273.15+Dnlim4(1)/(Dnlim4(2)+R*log(Wa_tot_Hi - Cr_tot_Hi/2));
        param_matrix(10,((iter-1)*3+1):((iter-1)*3+3)) = [Tm_Hi Tm_Hi_max Tm_Hi_min];        
        Tm_Lo = -273.15+thermvals(1)/(thermvals(2)+R*log(Wa_tot_Lo - Cr_tot_Lo/2));
        Tm_Lo_max = -273.15+Uplim4(1)/(Uplim4(2)+R*log(Wa_tot_Lo - Cr_tot_Lo/2));
        Tm_Lo_min = -273.15+Dnlim4(1)/(Dnlim4(2)+R*log(Wa_tot_Lo - Cr_tot_Lo/2));
        param_matrix(11,((iter-1)*3+1):((iter-1)*3+3)) = [Tm_Lo Tm_Lo_max Tm_Lo_min];  
        % Calculate max percent error in DH0 and DS0
        dh0lim = max([Uplim4(1)/thermvals(1), thermvals(1)/Dnlim4(1)]);
        ds0lim = max([Uplim4(2)/thermvals(2), thermvals(2)/Dnlim4(2)]);

%
% This concludes the computation of parameters


%Now optimize concentrations as well. Comment out lines below if needed
%if(iter > 1)

% Do four different fits for four possible choices for Wa and Cr for Hi and
% Lo melts. (This would not scale very well for multiple melts...)
% For now at least iterate all the way through making the plots and let the
% use evaluate which is best.

%-----------------
    else  % iteration > 2
        % Launch each one independently from best estimates for iteration 2 --
        % Only one of iterations 3-6 can actually be true!
 %       size(thermvals)
 %       size(Wa_tot_Hi)
        vconc = [thermvals , Wa_tot_Hi , Cr_tot_Hi , Wa_tot_Lo , Cr_tot_Lo ]; 
        % Upper and lower bounds on epsilDS and mds are somewhat arbitrary
        % but generous I think
        lb = [-300000,-1000,0.5*min([p1(1),p2(1)]),0.5*min([p1(2),p2(2)]),...
            Wa_tot_Hi ,Cr_tot_Hi_init*0.80, Wa_tot_Lo ,Cr_tot_Lo_init*0.80];
        ub = [0,0,2*2*max([p1(1),p2(1)]),max([p1(2),p2(2)]),Wa_tot_Hi_init*1.2,...
            Cr_tot_Hi,Wa_tot_Lo_init*1.2,Cr_tot_Lo ]; 

        [morevals,resnorm, derres2,~,output,~,jacobian2] = ...
        lsqnonlin(@Apred_concvar,vconc,lb,ub,options,...
            meltca,RefT);
        disp([ 'Iteration # ',num2str(iter),' Resnorm with refined baseline lims is ',...
            num2str(resnorm),' after ',num2str(output.iterations),...
            ' lsqnonlin iterations']);

% Should collect all of these values in an array

try
        civ = nlparci(morevals,derres2,jacobian2);
        Uplim8 = civ(:,1);
        Dnlim8 = civ(:,2);
catch
    disp(['Nlparci failed for iteration ',num2str(iter),'. Setting params to best fit.'])
    Uplim8 = morevals;
    Dnlim8 = morevals;
end
        
        Wa_tot_Hi = morevals(5);
        Cr_tot_Hi = morevals(6);
        Wa_tot_Lo = morevals(7);
        Cr_tot_Lo = morevals(8);
        param_matrix(1:4,(iter-1)*3+1) = morevals(1:4);
        param_matrix(1:4,(iter-1)*3+2) = Uplim8(1:4);
        param_matrix(1:4,(iter-1)*3+3) = Dnlim8(1:4);
        Up_Wa_tot_Hi = Uplim8(5);
        Up_Cr_tot_Hi = Uplim8(6);
        Dn_Wa_tot_Hi = Dnlim8(5);
        Dn_Cr_tot_Hi = Dnlim8(6);
        Up_Wa_tot_Lo = Uplim8(7);
        Up_Cr_tot_Lo = Uplim8(8);
        Dn_Wa_tot_Lo = Uplim8(7);
        Dn_Cr_tot_Lo = Uplim8(8);
        Abs_Excess_ssDNA_hiCmelt = (Wa_tot_Hi - Cr_tot_Hi)*Epsil_Wa_vec_Hi;
        Abs_pairedDNA_hiCmelt = Abs_hiCmelt - Abs_Excess_ssDNA_hiCmelt;
        Abs_Excess_ssDNA_loCmelt = (Wa_tot_Lo - Cr_tot_Lo)*Epsil_Wa_vec_Lo;
        Abs_pairedDNA_loCmelt = Abs_loCmelt - Abs_Excess_ssDNA_loCmelt;


        if (Wa_is_ss1_Hi)
            param_matrix(5:6,(iter-1)*3+1) = [Wa_tot_Hi Cr_tot_Hi];
            param_matrix(5:6,(iter-1)*3+2) = [Uplim8(5) Uplim8(6)];
            param_matrix(5:6,(iter-1)*3+3) = [Dnlim8(5) Dnlim8(6)];
        else
            param_matrix(5:6,(iter-1)*3+1) = [Cr_tot_Hi Wa_tot_Hi];
            param_matrix(5:6,(iter-1)*3+2) = [Uplim8(6) Uplim8(5)];
            param_matrix(5:6,(iter-1)*3+3) = [Dnlim8(6) Dnlim8(5)];
        end
        if (Wa_is_ss1_Lo)
            param_matrix(7:8,(iter-1)*3+1) = [Wa_tot_Lo Cr_tot_Lo];
            param_matrix(7:8,(iter-1)*3+2) = [Uplim8(7) Uplim8(8)];
            param_matrix(7:8,(iter-1)*3+3) = [Dnlim8(7) Dnlim8(8)];
        else
            param_matrix(7:8,(iter-1)*3+1) = [Cr_tot_Lo Wa_tot_Lo];
            param_matrix(7:8,(iter-1)*3+2) = [Uplim8(8) Uplim8(7)];
            param_matrix(7:8,(iter-1)*3+3) = [Dnlim8(8) Dnlim8(7)];
        end
% Calculate derived parameters:
        DH0 = morevals(1);
        DS0 = morevals(2);
        mds = morevals(3);
        Epsil_WC_RefT = morevals(4);

        DG037 = morevals(1)-310.15*morevals(2);
        DG037_min = Uplim8(1)-310.15*Uplim8(2);
        DG037_max = Dnlim8(1)-310.15*Dnlim8(2);
        param_matrix(9,((iter-1)*3+1):((iter-1)*3+3)) = [DG037 DG037_max DG037_min];        
        Tm_Hi = -273.15+morevals(1)/(morevals(2)+R*log(Wa_tot_Hi - Cr_tot_Hi/2));
        Tm_Hi_max = -273.15+Uplim8(1)/(Uplim8(2)+R*log(Uplim8(5) - Uplim8(6)/2));
        Tm_Hi_min = -273.15+Dnlim8(1)/(Dnlim8(2)+R*log(Dnlim8(5) - Dnlim8(6)/2));
        param_matrix(10,((iter-1)*3+1):((iter-1)*3+3)) = [Tm_Hi Tm_Hi_max Tm_Hi_min];        
        Tm_Lo = -273.15+morevals(1)/(morevals(2)+R*log(Wa_tot_Lo - Cr_tot_Lo/2));
        Tm_Lo_max = -273.15+Uplim8(1)/(Uplim8(2)+R*log(Uplim8(7) - Uplim8(7)/2));
        Tm_Lo_min = -273.15+Dnlim8(1)/(Dnlim8(2)+R*log(Uplim8(7) - Uplim8(8)/2));
        param_matrix(11,((iter-1)*3+1):((iter-1)*3+3)) = [Tm_Lo Tm_Lo_max Tm_Lo_min]; 
                % Calculate max percent error in DH0 and DS0
        dh0lim = max([Uplim8(1)/morevals(1), morevals(1)/Dnlim8(1)]);
        ds0lim = max([Uplim8(2)/morevals(2), morevals(2)/Dnlim8(2)]);

%
        resarray(iter) = resnorm;
        conc_err_matrix (1:4,iter) = [
            (param_matrix(5,(iter-1)*3+1) - C1_tot_Himelt_init)/C1_tot_Himelt_init ...
            (param_matrix(6,(iter-1)*3+1) - C2_tot_Himelt_init)/C2_tot_Himelt_init ...
            (param_matrix(7,(iter-1)*3+1) - C1_tot_Lomelt_init)/C1_tot_Lomelt_init ...
            (param_matrix(8,(iter-1)*3+1) - C2_tot_Lomelt_init)/C2_tot_Lomelt_init];
        verrconc = vecnorm(conc_err_matrix);

%---------------------
% nlparci provides 95% confidence limits derived from the output of
% lsqnonlin
    end
    

    
% Call the plotting routine. At this point param_matrix is being filled
% so we can pull DH and DS etc. from it.

	plotall
% Return Fitval structure so it can be propagated up to the command line if desired
% The plot file names will have been assigned and stored during the calls to plotall
% For now we have not returned graphic handles.
Fitvals.Params = param_matrix;
Fitvals.Resnorms = resarray;
Fitvals.Conc_errors = conc_err_matrix;
    

% Save the residuals of the first iteration
%     if Iter == 1;
%         derres1 = derres;
%     end
% end
% end
% if Iter == 3 ;
% Let four parameters float, make no assumptions about
% fds. But start with previous best guess. Tried it with floating all six,
% but it fails to converge. 
%    allvals0 = [thermvals mss mds Amin Amax];
%     allvals0 = [thermvals mss mds ];
%    [allvals,resnormall, derresall,exitflag,output,lambda,jacobian3] = ...
%      lsqcurvefit(@Dpredall,allvals0,TK,Ader,[],[],options,CT,ExW); 
%     [allvals,resnorm, derresall,exitflag,output,lambda,jacobian3] = ...
%             lsqcurvefit(@Dpredall,allvals0,TK,Ader,[],[],options,CT,ExW,Amin,Amax); 
%        allvals;
%       jacobian3;
% Note that if this doesn't converge, which it often won't, then
% lsqcurvefit does not return the Jacobian? In this case nlparci will fail
% in a rather pernicious way -- it will give meaningless numbers. No, this
% doesn't seem to be the problem. The optimization for some reason just 
% doesn't give sensible conf. limits even though the residual has dropped.
%    ci = nlparci(allvals,derresall,jacobian3);
%    disp([ 'Unrestrained resnorm (', num2str(Iter), ') is ',num2str(resnorm),...
%' after ',...
%      num2str(output.iterations),' iterations']);
%    disp(' ');
%mds = thermvals(3);
%Epsil_WC_RefT = thermvals(4);
% The Tm is needed in the next iteration
%Tm_Hi = -273.15+thermvals(1)/(thermvals(2)+R*log(Wa_tot_Hi - Cr_tot_Hi/2));

%Amin = allvals(5);
%Amax = allvals(6);
% if (mss < 0 && mss*(min(T)-max(T)) > 0.05*(Amax-Amin))
%     disp([ 'Negative mss too large ',num2str(mss), 'at Iter ',num2str(Iter)]);
%     flag(2) = 3;
% end
% fdsfinHi = Fpred(thermvals,T_hiCmelt+273.15,Wa_tot_Hi,Cr_tot_Hi);
% dsnum_Hi = find(fdsfinHi > 0.98,1,'last') ;
% ssnum_Hi = find(fdsfinHi < 0.02,1,'first') ;
% fdsfinLo = Fpred(thermvals,T_loCmelt+273.15,Wa_tot_Lo,Cr_tot_Lo);
% dsnum_Lo = find(fdsfinLo > 0.98,1,'last') ;
% ssnum_Lo = find(fdsfinLo < 0.02,1,'first') ;
% if ssnum_Hi < 10
%     disp([ 'Melts too high; only  ',num2str(ssnum_Hi),' ssDNA points for Hi melt']);
%     flag(1) = 3;
% end
% if dsnum_Lo < dslim
%     disp([ 'ds baseline too short;, only ',...
%       num2str(dsnum_Lo),' dsDNA points for Lo melt']);
%     flag(1) = 2;
% end
% Put code to recalculate ssDNA concentrations in each melt here at the end
% of the for loop


%Fitvals.Iter = iter;
% The stuff below gets done only once, after we are done iterating over [Wa] and [Cr] 
%thermvals = lsqcurvefit(@Apred,guess(2),T,Absor,[],[],options,DH,CT,ExW,mss,mds,Amin,Amax)
%,CT) %,mss,mds)

%Fitvals.Name = pfile;

%
% Return DH with 95% confidence limits, DS with conf. limits, resnorm
% Also return DG37 but with calc. conf. limits from min/max DH DS
% Uncomment lines below to print. Printing from the GUI doesn't work well.
%pfile = [filestem '.eps'];
%print( '-depsc2', '-r300', '-adobecset', '-tiff', pfile)  % With tiff preview.
%
% if (~nineplots)
%     pfile = [filestem '.png'];
%     export_fig(pfile,'-m2');
%     Fitvals.Name = pfile;
% else
%     Fitvals.Name = ['none'];
% end

%Fitvals = struct('Nums',{},'Outflags',{},'Fds',fds);
% Fitvals.Nums = [thermvals(1) ci(1,:) thermvals(2) ci(2,:) DG037 DG037_max DG037_min ...
%         Tm_Hi Tm_Hi_max Tm_Hi_min Tm_Lo Tm_Lo_max Tm_Lo_min...
%Wa_tot_Hi Cr_tot_Hi Wa_tot_Lo Cr_tot_Lo ...
%         mds Epsil_WC_RefT resnorm];
% Fitvals.Outflags =  flag;
% Fitvals.Ghandle = h10;
% Fitvals.Fds = [fds_Hi ; fds_Lo];
% Fitvals.Wa_is_ss1_Hi = Wa_is_ss1_Hi;
% Fitvals.Wa_is_ss1_Lo = Wa_is_ss1_Lo;
% if (debug > 1)
%     disp('Output structure fields are ')
%     disp(fieldnames(Fitvals))
% end
end
%
% Try to guess at the truth -- this logic is not absolute
% Resnorms (in resarray) within a factor of 3 of the best are probably actually similar
% The concentration changes (in conc_err_matrix and verrconc should be minimized as
% long as this constraint is met.
% Simplest option: Choose the minimum verrconc that has resnorm within the
% factor of "tol_res_ratio." If one resnorm is by far the best it will be
% the result returned.
tol_res_ratio = 3;  % Try 3 usually
min_res = min(resarray);
ktrial = find(resarray < tol_res_ratio*min_res);
[~,best_conc] = min(verrconc(ktrial));
best_trial = ktrial(best_conc);
h12 = figure;
h12.Name = filestem;
title('Best Guess at Best Trial')
yyaxis left
semilogy(1:6,resarray,'-o','MarkerSize',12)
hold on
semilogy(best_trial,resarray(best_trial),'gs','MarkerSize',18)
xlabel('Trial number')
xlim([0.5 6.5])
xticks(1:6)
ylabel('Fit Residual Norm')
yyaxis right
plot(1:6,verrconc,'-o','MarkerSize',12)
plot(best_trial,verrconc(best_trial),'bs','MarkerSize',18)
ylabel('Concentration Error Norm')
figure(fig_handle_mat(best_trial))
subplot(3,3,7)
text(10,7,'The Chosen One','Color','red','FontSize',14)
disp('')
disp(['Best Trial Believed to Be # ',num2str(best_trial)])
if (verrconc(best_trial) > Conc_err_concern)
    disp(['However, concerned about conc err norm being ',num2str(verrconc(best_trial))])
end


% These lines just print to screen
disp('Initial concentrations')
disp(['C1_tot_Himelt_init = ',num2str(C1_tot_Himelt_init)])
disp(['C1_tot_Lomelt_init = ',num2str(C1_tot_Lomelt_init)])
disp(['C2_tot_Himelt_init = ',num2str(C2_tot_Himelt_init)])
disp(['C2_tot_Lomelt_init = ',num2str(C2_tot_Lomelt_init)])
disp('Parameter matrix: ')
disp(param_matrix)
disp('Resnorms x100: ')
disp(100*resarray)
disp('Concentration changes (%): ')
disp(100*conc_err_matrix)
disp('Vector norm of concentration changes: ')
disp(verrconc)
disp('Print files: ')
for j=1:6
    disp(Pfiles{j})
end

%
% Plotall is a nested function so it inherits access to all the variables.
% Plots one run at a time to avoid storing a whole lot of temporary stuff
%
function plotall
        %Absorfit = Apred(thermvals,T,DH,CT,mss,mds,Amin,Amax);
%
%Now regenerate predicted values using the best-fit parameters and explore
%how parameters at the 95% confidence limit fit the data
%Need to adjust this for the different numbers of parameters
Abs_hiCmelt_fit = Apred(param_matrix(1:4,(iter-1)*3+1),T_hiCmelt+273.15,...
     Wa_tot_Hi,Cr_tot_Hi,Epsil_Wa_vec_Hi,Epsil_Cr_vec_Hi,0,ones(size(T_hiCmelt)),RefT);
Abs_loCmelt_fit = Apred(param_matrix(1:4,(iter-1)*3+1),T_loCmelt+273.15,...
     Wa_tot_Lo,Cr_tot_Lo,Epsil_Wa_vec_Lo,Epsil_Cr_vec_Lo,0,ones(size(T_loCmelt)),RefT);
Uplim_Abs_Hi = Apred(param_matrix(1:4,(iter-1)*3+2),T_hiCmelt+273.15,Up_Wa_tot_Hi,...
        Up_Cr_tot_Hi,Epsil_Wa_vec_Hi,Epsil_Cr_vec_Hi,0,ones(size(T_hiCmelt)),RefT);
%isreal(Uplim_Abs_Hi)
Uplim_Abs_Lo = Apred(param_matrix(1:4,(iter-1)*3+2),T_loCmelt+273.15,Up_Wa_tot_Lo,...
        Up_Cr_tot_Lo,Epsil_Wa_vec_Lo,Epsil_Cr_vec_Lo,0,ones(size(T_loCmelt)),RefT);
%isreal(Uplim_Abs_Lo)    
Dnlim_Abs_Hi = Apred(param_matrix(1:4,(iter-1)*3+3),T_hiCmelt+273.15,Dn_Wa_tot_Hi,...
        Dn_Cr_tot_Hi,Epsil_Wa_vec_Hi,Epsil_Cr_vec_Hi,0,ones(size(T_hiCmelt)),RefT);
%isreal(Dnlim_Abs_Hi)
Dnlim_Abs_Lo = Apred(param_matrix(1:4,(iter-1)*3+3),T_loCmelt+273.15,Dn_Wa_tot_Lo,...
        Dn_Cr_tot_Lo,Epsil_Wa_vec_Lo,Epsil_Cr_vec_Lo,0,ones(size(T_loCmelt)),RefT);
% If things are well-behaved, should not trigger the lines below
if(~isreal(Abs_hiCmelt_fit) || ~isreal(Abs_loCmelt_fit))
    disp('Even best fit params give complex absorbance prediction')
    disp(['Bailing out of plotting for iteration ',num2str(iter)])
    return
end    
if(~isreal(Uplim_Abs_Hi))
    Uplim_Abs_Hi = Abs_hiCmelt_fit;
    disp('Setting Uplim_Abs_Hi to best fit due to complexity')
end
if(~isreal(Uplim_Abs_Lo))
    Uplim_Abs_Lo = Abs_loCmelt_fit;
    disp('Setting Uplim_Abs_Lo to best fit due to complexity')
end
if(~isreal(Dnlim_Abs_Hi))
    Dnlim_Abs_Hi = Abs_hiCmelt_fit;
    disp('Setting DnLim_Abs_Hi to best fit due to complexity')
end
if(~isreal(Dnlim_Abs_Lo))
    Dnlim_Abs_Lo = Abs_loCmelt_fit;
    disp('Setting DnLim_Abs_Lo to best fit due to complexity')
end

% Preallocate arrays because Matlab likes that
    L1 = length(Abs_hiCmelt_fit);
    Abs_der_pred_Hi_Opt = zeros(1,L1-2);
    Abs_der_pred_Hi_Uplim = Abs_der_pred_Hi_Opt;
    Abs_der_pred_Hi_Dnlim = Abs_der_pred_Hi_Opt;
    L2 = length(Abs_loCmelt_fit);
    Abs_der_pred_Lo_Opt = zeros(1,L2-2);
    Abs_der_pred_Lo_Uplim = Abs_der_pred_Lo_Opt;
    Abs_der_pred_Lo_Dnlim = Abs_der_pred_Lo_Opt;
%
%
for I=2:L1-1
    Abs_der_pred_Hi_Opt(I-1) = ...
        (Abs_hiCmelt_fit(I-1)-Abs_hiCmelt_fit(I+1))/(T_hiCmelt(I-1) - T_hiCmelt(I+1));
    Abs_der_pred_Hi_Uplim(I-1) = ...
        (Uplim_Abs_Hi(I-1)-Uplim_Abs_Hi(I+1))/(T_hiCmelt(I-1) - T_hiCmelt(I+1));
    Abs_der_pred_Hi_Dnlim(I-1) = ...
        (Dnlim_Abs_Hi(I-1)-Dnlim_Abs_Hi(I+1))/(T_hiCmelt(I-1) - T_hiCmelt(I+1));
end
for I=2:L2-1
    Abs_der_pred_Lo_Opt(I-1) = ...
        (Abs_loCmelt_fit(I-1)-Abs_loCmelt_fit(I+1))/(T_loCmelt(I-1) - T_loCmelt(I+1));
    Abs_der_pred_Lo_Uplim(I-1) = ...
        (Uplim_Abs_Lo(I-1)-Uplim_Abs_Lo(I+1))/(T_loCmelt(I-1) - T_loCmelt(I+1));
    Abs_der_pred_Lo_Dnlim(I-1) = ...
        (Dnlim_Abs_Lo(I-1)-Dnlim_Abs_Lo(I+1))/(T_loCmelt(I-1) - T_loCmelt(I+1));   
end
%Abs_der_pred_Lo_Uplim
%Abs_der_pred_Lo_Dnlim
%
% Plot results: 
%
%if (flag(1) ~= 1 || flag(2) ~= 1)
%     nineplots = true;
% else
%    disp('There is some issue with the baselines')
%end
%scrsz = get(0,'ScreenSize');
h10 = figure;
fig_handle_mat(iter) = h10;
% set(h9,'Position',[1+(Iter-1)*100 (scrsz(4)/3)-(Iter-1)*200 ...
%         2*scrsz(3)/3 2*scrsz(4)/3])
h10.Units = 'inches';
h10.Position = [1 1 18 13];

%if (nineplots)
    h1 = subplot(3,3,1);
%else
%    h1 = subplot(3,2,1);
%end
%
% Top left subplot
% DNA absorbance vs. T, fit vs. expt, for both concentrations and refs
%
hold on
box on
plot(T_hiCmelt,Abs_hiCmelt,'o','Color',Hi_color)
plot(T_hiCmelt,Abs_hiCmelt_fit,'k','LineWidth',1)
plot(T_loCmelt,Abs_loCmelt,'o','Color',Lo_color)
plot(T_loCmelt,Abs_loCmelt_fit,'k','LineWidth',1)
if (Wa_is_ss1_Hi)
    plot(T_hiCmelt,Abs_Cr_ref_interp_Hi,'Color','m')
    plot(T_hiCmelt,Abs_Wa_ref_interp_Hi,'Color',[0.5 0.7 0.9])
else
    plot(T_hiCmelt,Abs_Wa_ref_interp_Hi,'Color','m')
    plot(T_hiCmelt,Abs_Cr_ref_interp_Hi,'Color',[0.5 0.7 0.9])
end
Amin = min(Abs_loCmelt);
Amax = max(Abs_hiCmelt);
% Can use    'FontSize',14,'FontWeight', 'bold')
text(min(T_hiCmelt),Amax+.05*(Amax-Amin),[' DH0 = ',num2str(DH0)],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',12)
text(min(T_hiCmelt),Amax-.15*(Amax-Amin),[' DS0 = ',num2str(DS0)],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',12)
text(min(T_hiCmelt),Amax-.40*(Amax-Amin),[' DG037 = ',num2str(DG037)],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',12)
text(min(T_hiCmelt),Amax-.70*(Amax-Amin),[' TM1 = ',num2str(Tm_Hi)],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',12)
text(min(T_hiCmelt),Amax-.90*(Amax-Amin),[' TM2 = ',num2str(Tm_Lo)],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',12)
text(max(T_hiCmelt),Amin+.50*(Amax-Amin),'Strand #1 ref','Color',[0.5 0.7 0.9],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','right',...
    'FontSize',12)
text(max(T_hiCmelt),Amin+.40*(Amax-Amin),'Strand #2 ref','Color','m',...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','right',...
    'FontSize',12)
text(max(T_hiCmelt),Amin+.25*(Amax-Amin),['DH%err = ',num2str(100*(dh0lim-1))],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','right',...
    'FontSize',12)
text(max(T_hiCmelt),Amin-0.05*(Amax-Amin),['DS%err  =',num2str(100*(ds0lim-1))],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','right',...
    'FontSize',12)
% Axis statement must follow text statements or it is ignored ??!
axis(h1,[min([T_hiCmelt ; T_loCmelt])-1 ...
    max([T_hiCmelt ; T_loCmelt])+1 ...
    Amin-0.15*(Amax-Amin) Amax+0.15*(Amax-Amin)]) 
title('Experimental Absorbances vs. T','FontSize',14)
ylabel('Absorbance')
%
% Derivative plots
%
% if (nineplots)
    h2=subplot(3,3,5);
% else
%     h2=subplot(3,2,3);
% end
L = length(T_hiCmelt);
hold on
box on
plot(T_hiCmelt(2:L-1),Abs_der_pred_Hi_Opt,'g','LineWidth', 2)
plot(T_hiCmelt(2:L-1),Abs_der_Hi,'-o','Color',Hi_color)
plot(T_hiCmelt(2:L-1),Abs_der_pred_Hi_Uplim,'r--','LineWidth',1)
plot(T_hiCmelt(2:L-1),Abs_der_pred_Hi_Dnlim,'b--','LineWidth',1)
text(min(T_hiCmelt),max(Abs_der_Hi),[' Res = ',num2str(resnorm)],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',12)
text(min(T_hiCmelt),max(Abs_der_Hi)-0.70*(max(Abs_der_Hi)-min(Abs_der_Hi)),...
    [' File: ',filestem],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',12,'Interpreter','none')  % To allow underbar to be displayed
miny = min(Abs_der_Hi)-0.15*max(Abs_der_Hi);
maxy = max(Abs_der_Hi)*1.15;
plot([Tm_Hi Tm_Hi],[miny maxy],'k--')
axis(h2,[min(T_hiCmelt)-2 max(T_hiCmelt)+2 miny maxy])
title('Deriv. Absor. vs. T, High Conc','FontSize',14)
ylabel('dA260/dT')
%
%
%
% if (nineplots)
    h3=subplot(3,3,6);
% else
%     h2=subplot(3,2,3);
% end
L = length(T_loCmelt);
hold on
box on
if (debug >= 2)
    size(T_loCmelt)
    size(Abs_der_pred_Lo_Opt)
    size(Abs_der_Lo)
    size(Abs_der_pred_Lo_Uplim)
    size(Abs_der_pred_Lo_Dnlim)
end
plot(T_loCmelt(2:L-1),Abs_der_pred_Lo_Opt,'g','LineWidth', 2)
plot(T_loCmelt(2:L-1),Abs_der_Lo,'-o','Color',Lo_color)
plot(T_loCmelt(2:L-1),Abs_der_pred_Lo_Uplim,'r--','LineWidth',1)
plot(T_loCmelt(2:L-1),Abs_der_pred_Lo_Dnlim,'b--','LineWidth',1)
% text(min(T_loCmelt),max(Abs_der_Lo),[' Res = ',num2str(resnorm)],...
%     'VerticalAlignment','middle',...
%     'HorizontalAlignment','left',...
%     'FontSize',12)
% text(min(T_loCmelt),max(Abs_der_Lo)-0.25*(max(Abs_der_Lo)-min(Abs_der_Lo)),...
%     [' FILESTEM: ',filestem],...
%     'VerticalAlignment','middle',...
%     'HorizontalAlignment','left',...
%     'FontSize',12,'Interpreter','none')  % To allow underbar to be displayed
if (debug == 2)
    disp(['min(T_loCmelt) - 2 is ',num2str(min(T_loCmelt) - 2)])
    disp(['max(T_loCmelt) + 2 is ',num2str(max(T_loCmelt) + 2)])
    disp(['min([Abs_der_Lo  Abs_der_pred_Lo_Dnlim]) - 0.15*max(Abs_der_Lo) is ',...
        num2str(min([Abs_der_Lo  Abs_der_pred_Lo_Dnlim]) - 0.15*max(Abs_der_Lo))])
    disp(['max([Abs_der_Lo  Abs_der_pred_Lo_Uplim]) * 1.15 is ',...
        num2str(max([Abs_der_Lo  Abs_der_pred_Lo_Uplim]) * 1.15)])
end
%Abs_der_Lo
miny = min([Abs_der_Lo Abs_der_pred_Lo_Dnlim])-0.15*max(Abs_der_Lo);
maxy = max([Abs_der_Lo Abs_der_pred_Lo_Uplim])*1.15;
plot([Tm_Lo Tm_Lo],[miny maxy],'k--')
axis(h3,[min(T_loCmelt)-2 max(T_loCmelt)+2 miny maxy])
title('Deriv. Absor. vs. T, Low Conc','FontSize',14)
%
% Residuals
%
Der_res_Hi = Abs_der_pred_Hi_Opt - Abs_der_Hi;
Der_res_Lo = Abs_der_pred_Lo_Opt - Abs_der_Lo;
%if (nineplots)
    h4 = subplot(3,3,8);
% else
%     h7 = subplot(3,2,5);
% end
lorax = min(Der_res_Hi)-0.1*max(Der_res_Hi);
hirax = max(Der_res_Hi)*1.1;
hold on
box on
plot(T_hiCmelt(2:(length(T_hiCmelt)-1)),Der_res_Hi,'g','LineWidth', 1);
plot(T_hiCmelt(2:(length(T_hiCmelt)-1)),Abs_der_pred_Hi_Uplim - Abs_der_Hi,'r+');
plot(T_hiCmelt(2:(length(T_hiCmelt)-1)),Abs_der_pred_Hi_Dnlim - Abs_der_Hi,'b+');
xmin = min(T_hiCmelt)-2;
xmax = max(T_hiCmelt)+2;
plot([xmin xmax],[0 0],'k--')
axis (h4,[xmin xmax lorax hirax] )
title('Deriv. residuals, High CT melt','FontSize',14)
xlabel('Temperature (°C)')

% if (nineplots)
    h5 = subplot(3,3,9);
% else
%     h8 = subplot(3,2,5);
% end
lorax = min(Der_res_Lo)-0.1*max(Der_res_Lo);
hirax = max(Der_res_Lo)*1.1;
hold on
box on
plot(T_loCmelt(2:(length(T_loCmelt)-1)),Der_res_Lo,'g','LineWidth', 1);
plot(T_hiCmelt(2:(length(T_loCmelt)-1)),Abs_der_pred_Lo_Uplim - Abs_der_Lo,'r+');
plot(T_hiCmelt(2:(length(T_loCmelt)-1)),Abs_der_pred_Lo_Dnlim - Abs_der_Lo,'b+');
xmin = min(T_loCmelt)-2;
xmax = max(T_loCmelt)+2;
plot([xmin xmax],[0 0],'k--')
axis (h5,[xmin xmax lorax hirax] )

% end
% if Iter == 2;
%     plot(T(2:L-1),derres,'o-',T(2:L-1),derres1,'+','LineWidth', 1);    
%     lorax = min([min(derres1)-0.1*max(derres1) min(derres)-0.1*max(derres)]);
%     hirax = max([max(derres1)*1.1 max(derres)*1.1]);
%     axis (h7,[min(T)-2 max(T)+2 lorax hirax])
% end
% if Iter == 3;
%     plot(T(2:L-1),derresall,'-s',T(2:L-1),derres,'o',T(2:L-1),derres1,'+','LineWidth',1);    
%     lorax = min([min(derres1)-0.1*max(derres1) min(derres)-0.1*max(derres)]);
%     lorax = min([lorax min(derresall)-0.1*max(derresall)]);
%     hirax = max([max(derres1)*1.1 max(derres)*1.1]);
%     hirax = max([hirax max(derresall)*1.1]);
% end
%axis (h7,[min(T)-2 max(T)+2 lorax hirax])
title('Deriv. residuals, Low CT melt','FontSize',14)
xlabel('Temperature (°C)')
%
% Component absorbances
%
% if (nineplots)
    h6=subplot(3,3,2);
% else
%     h2=subplot(3,2,2);
% end
hold on
box on
plot(T_hiCmelt,Abs_hiCmelt,'o','Color',Hi_color)
plot(T_hiCmelt,Abs_pairedDNA_hiCmelt,'Linewidth',2)
plot(T_hiCmelt,DSBpred(T_hiCmelt,RefT,Cr_tot_Hi,mds,Epsil_WC_RefT),'b--','LineWidth',1)
plot(T_hiCmelt,Abs_hiCmelt_fit,'k','LineWidth',1)

% This line is what we would see if it were all ssDNA at the extant
% concentrations
plot(T_hiCmelt,(Epsil_Wa_vec_Hi*Wa_tot_Hi + Epsil_Cr_vec_Hi*Cr_tot_Hi),'Color',[0.5 0.5 0.8],'LineWidth',1)
%plot(T,Absor+factor*Arefi,'o',T,Absorfit,'LineWidth', 1) %+factor*Arefi
%plot(T,HiAbsor,'--',T,LoAbsor,'--','LineWidth', 0.5)
Amax = max(Abs_hiCmelt);
Amin = min(Abs_pairedDNA_hiCmelt);
text(min(T_hiCmelt),Amax,[' DHmax = ',num2str(param_matrix(1,(iter-1)*3+2))],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',12)
text(min(T_hiCmelt),0.80*Amax+0.2*Amin,...
    [' DSmax = ',num2str(param_matrix(2,(iter-1)*3+2))],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',12)
% text(min(T_hiCmelt),0.60*Amax,[' Reference spectrum method '],...
%     'VerticalAlignment','middle',...
%     'HorizontalAlignment','left',...
%     'FontSize',12)
if Wa_is_ss1_Hi
     cwa = [0.5 0.7 .9];
     ccr = 'm';
     identwa = 'ss1';
     identcr = 'ss2';
else
     cwa = 'm';
     ccr = [0.5 0.7 .9];
     identwa = 'ss2';
     identcr = 'ss1';
end
text(min(T_hiCmelt),0.60*Amax+0.4*Amin,[' [Wa(',identwa,')] = ',num2str(Wa_tot_Hi)],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',12,'Color',cwa)
text(min(T_hiCmelt),0.55*Amax+0.4*Amin,[' [Cr(',identcr,')] = ',num2str(Cr_tot_Hi)],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',12,'Color',ccr)
text(max(T_hiCmelt),0.20*(Amax-Amin)+1.05*Amin,...
    [' DHmin = ',num2str(param_matrix(1,(iter-1)*3+3))],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','right',...
    'FontSize',12)
text(max(T_hiCmelt),1.05*Amin,[' DSmin = ',num2str(param_matrix(2,(iter-1)*3+3))],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','right',...
    'FontSize',12)
title('Component Absorbances vs. T','FontSize',14)
%[min(T)-2 max(T)+2 Amin-0.2*(Amax-Amin) Amax+0.2*(Amax-Amin)]
plot([Tm_Hi Tm_Hi],[Amin*0.95 Amax+0.1*(Amax-Amin)],'k--')
axis(h6,[min(T_hiCmelt)-2 max(T_hiCmelt)+2 Amin*0.95 Amax+0.1*(Amax-Amin)]) 
%
%if (nineplots)
    h7=subplot(3,3,3);
% else
%     h3=subplot(3,2,2);
% end
hold on
box on
plot(T_loCmelt,Abs_loCmelt,'o','Color',Lo_color)
plot(T_loCmelt,Abs_pairedDNA_loCmelt,'Linewidth',2)
plot(T_loCmelt,DSBpred(T_loCmelt,RefT,Cr_tot_Lo,mds,Epsil_WC_RefT),'b--','LineWidth',1)
plot(T_loCmelt,Abs_loCmelt_fit,'k','LineWidth',1)

% plot(T_loCmelt,Abs_Cr_ref_interp_Lo)
% plot(T_loCmelt,Abs_Wa_ref_interp_Lo)
% This line is what we would see if it were all ssDNA at the extant
% concentrations
plot(T_loCmelt,(Epsil_Wa_vec_Lo*Wa_tot_Lo + Epsil_Cr_vec_Lo*Cr_tot_Lo),'Color',[0.5 0.5 0.8],'LineWidth',1)
%plot(T,Absor+factor*Arefi,'o',T,Absorfit,'LineWidth', 1) %+factor*Arefi
%plot(T,HiAbsor,'--',T,LoAbsor,'--','LineWidth', 0.5)
Amax = max(Abs_loCmelt);
Amin = min(Abs_loCmelt-Abs_Excess_ssDNA_loCmelt);
%Amin = min(Abs_Cr_ref_interp_Lo);
% text(min(T_loCmelt),Amax,[' DHmax = ',num2str(hilim(1))],...
%     'VerticalAlignment','middle',...
%     'HorizontalAlignment','left',...
%     'FontSize',12)
% text(min(T_loCmelt),0.80*Amax,[' DSmax = ',num2str(hilim(2))],...
%     'VerticalAlignment','middle',...
%     'HorizontalAlignment','left',...
%     'FontSize',12)
% text(min(T_loCmelt),0.60*Amax,[' Reference spectrum method '],...
%     'VerticalAlignment','middle',...
%     'HorizontalAlignment','left',...
%     'FontSize',12)
if Wa_is_ss1_Lo
     cwa = [0.5 0.7 0.9];
     ccr = 'm';
     identwa = 'ss1';
     identcr = 'ss2';
else
     cwa = 'm';
     ccr = [0.5 0.7 0.9];
     identwa = 'ss2';
     identcr = 'ss1';
end
text(max(T_loCmelt),0.40*Amax+0.6*Amin,[' [Wa(',identwa,')] = ',num2str(Wa_tot_Lo)],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','right',...
    'FontSize',12,'Color',cwa)
text(max(T_loCmelt),0.35*Amax+0.6*Amin,[' [Cr(',identcr,')] = ',...
    num2str(Cr_tot_Lo)],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','right',...
    'FontSize',12,'Color',ccr)
% text(max(T_loCmelt),0.25*Amax,[' DHmin = ',num2str(lolim(1))],...
%     'VerticalAlignment','middle',...
%     'HorizontalAlignment','right',...
%     'FontSize',12)
% text(max(T_loCmelt),0.05*Amax,[' DSmin = ',num2str(lolim(2))],...
%     'VerticalAlignment','middle',...
%     'HorizontalAlignment','right',...
%     'FontSize',12)
title('Component Absorbances, Low CT','FontSize',14)
%[min(T)-2 max(T)+2 Amin-0.2*(Amax-Amin) Amax+0.2*(Amax-Amin)]
plot([Tm_Lo Tm_Lo],[Amin*0.95 Amax*1.05],'k--')
axis(h7,[min(T_loCmelt)-2 max(T_loCmelt)+2 Amin*0.95 Amax*1.05]) 
%
% Compare confidence limits on derivative plot
%
% if (nineplots)
%     h5=subplot(3,3,5);
% else
%     h5=subplot(3,2,4);
% end
% plot(T(2:L-1),Ader,'o',T(2:L-1),Apredder,'LineWidth', 1)
% hold on
% plot(T(2:L-1),Hider,'--',T(2:L-1),Loder,'--','LineWidth', 0.5)
% text(min(T),max(Ader),[' Hinorm = ',num2str(Hinorm)],...
%     'VerticalAlignment','middle',...
%     'HorizontalAlignment','left',...
%     'FontSize',12)
% text(max(T),min(Ader)*0.95,[' Lonorm = ',num2str(Lonorm)],...
%     'VerticalAlignment','middle',...
%     'HorizontalAlignment','right',...
%     'FontSize',12)
% axis(h5,[min(T)-2 max(T)+2 min(Ader)-0.15*max(Ader) max(Ader)*1.15])
% title('Diff. Der. Absor. vs. T fits','FontSize',14)
%
% Calculated fraction ds and fraction ss
%
%if (nineplots)
    h8 = subplot(3,3,4);
% else
%     h6 = subplot(3,2,6);
% end
fdsexpt_Hi = Fdsexpt(Abs_hiCmelt-Abs_Excess_ssDNA_hiCmelt,T_hiCmelt,RefT,mds*Cr_tot_Hi,...
    Epsil_WC_RefT*Cr_tot_Hi,Cr_tot_Hi*Epsil_Wa_vec_Hi+Cr_tot_Hi*Epsil_Cr_vec_Hi);
fdsexpt_Lo = Fdsexpt(Abs_loCmelt-Abs_Excess_ssDNA_loCmelt,T_loCmelt,RefT,mds*Cr_tot_Lo,...
    Epsil_WC_RefT*Cr_tot_Lo,Cr_tot_Lo*Epsil_Wa_vec_Lo+Cr_tot_Lo*Epsil_Cr_vec_Lo);
% Need min and max b/c sometimes fds gets < 0 or > 1 due to noise in
% absorbance
minex = min([fdsexpt_Lo ; fdsexpt_Hi]);
maxex = max([fdsexpt_Hi ; fdsexpt_Lo]);
%Fpred does not need mds and Epsil but should ignore them
fds_Hi = Fpred(param_matrix(1:2,(iter-1)*3+1),T_hiCmelt+273.15,Wa_tot_Hi,Cr_tot_Hi);
dsmax_Hi = find(fds_Hi > 0.99,1,'last');
ssmin_Hi = find(fds_Hi < 0.01,1,'first');
fds_Lo = Fpred(param_matrix(1:2,(iter-1)*3+1),T_loCmelt+273.15,Wa_tot_Lo,Cr_tot_Lo);
dsmax_Lo = find(fds_Lo > 0.99,1,'last');
ssmin_Lo = find(fds_Lo < 0.01,1,'first');
hold on
box on
plot(T_hiCmelt,fds_Hi,'Color',Hi_color,'LineWidth',1)
plot(T_loCmelt,fds_Lo,'Color',Lo_color,'LineWidth',1)
plot (T_hiCmelt, fdsexpt_Hi,'o','Color',Hi_color)
plot (T_loCmelt, fdsexpt_Lo,'o','Color',Lo_color)

text(min(T_hiCmelt),0.85,[' All ds to ',num2str(T_loCmelt(dsmax_Lo)),' and ',...
    num2str(T_hiCmelt(dsmax_Hi))],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',12)
text(max(T_hiCmelt),0.15,[' All ss from ',num2str(T_loCmelt(ssmin_Lo)),' and ',...
    num2str(T_hiCmelt(ssmin_Hi))],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','right',...
    'FontSize',12)
% text(min(T),0.35,[' Iter ',num2str(Iter)],...
%     'VerticalAlignment','middle',...
%     'HorizontalAlignment','left',...
%     'FontSize',12)
axis(h8,[min(T_hiCmelt)-2 max(T_hiCmelt)+2 min([-0.05 minex]) max([1.05 maxex])])
%
% if Iter == 1;
%     title('Fraction ds, preset limits','FontSize',14)
% end
% if Iter == 2;
%     title('Fraction ds, adaptive limits','FontSize',14)
% end    
% if Iter == 3;
    title('Fraction double stranded','FontSize',14)
% end
xlabel('Temperature (°C)')

%fp = figure(Fitvals.Ghandle);
h10.Name = [filestem '_n' num2str(iter) '_sf' num2str(scale_factor)];
if (scale_factor == 1)
    Pfiles{iter} = h10.Name;
else
    Pfiles{iter} = [filestem '_n' num2str(iter)];
end
subplot(3,3,4)
text(Tmin+5,0.6,['Orig Tlow = ',num2str(Tlow)],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',12)
text(Tmin+5,0.45,[' TMin = ',num2str(Tmin)],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',12)
text(Tmin+5,0.3,[' TMax = ',num2str(Tmax)],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',12)
text(Tmin+5,0.15,['Orig Thigh = ',num2str(Thigh)],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',12)
subplot(3,3,7)
hold on
box on
title('Info on Fitting Process','FontSize',14)
% How could this have worked??? This function should not have access to
% optargs or numvarargs
if (scale_flag == 1)
    text(10,95,['Scale factor input as ',num2str(scale_factor)]);
else
    text(10,95,['Scale factor calculated at ',num2str(scale_factor)]);
end
% text(10,80,['Number of iterations is  ',num2str(Fitvals.Iter),...
%     ' of ',num2str(numtries),' permitted']);
text(10,86,['Nominal Strand 1 []s = ',num2str(C1_tot_Himelt_init),' and ',num2str(C1_tot_Lomelt_init),' M'])
text(10,68,['Nominal Strand 2 []s = ',num2str(C2_tot_Himelt_init),' and ',num2str(C2_tot_Lomelt_init),' M'])
if Wa_is_ss1_Hi && Wa_is_ss1_Lo
    text(10,77,['BestEst Strand 1 []s = ',num2str(Wa_tot_Hi),' and ',num2str(Wa_tot_Lo),' M'])
    text(10,59,['BestEst Strand 2 []s = ',num2str(Cr_tot_Hi),' and ',num2str(Cr_tot_Lo),' M'])
elseif Wa_is_ss1_Hi && ~Wa_is_ss1_Lo
    text(10,77,['BestEst Strand 1 []s = ',num2str(Wa_tot_Hi),' and ',num2str(Cr_tot_Lo),' M'])
    text(10,59,['BestEst Strand 2 []s = ',num2str(Cr_tot_Hi),' and ',num2str(Wa_tot_Lo),' M'])
elseif ~Wa_is_ss1_Hi && Wa_is_ss1_Lo
    text(10,77,['BestEst Strand 1 []s = ',num2str(Cr_tot_Hi),' and ',num2str(Wa_tot_Lo),' M'])
    text(10,59,['BestEst Strand 2 []s = ',num2str(Wa_tot_Hi),' and ',num2str(Cr_tot_Lo),' M'])
elseif ~Wa_is_ss1_Hi && ~Wa_is_ss1_Lo
    text(10,77,['BestEst Strand 1 []s = ',num2str(Cr_tot_Hi),' and ',num2str(Cr_tot_Lo),' M'])
    text(10,59,['BestEst Strand 2 []s = ',num2str(Wa_tot_Hi),' and ',num2str(Wa_tot_Lo),' M'])
else
    disp('Logic Fault with Was and Crs')
end
text(10,50,['Norm of conc change from nominal = ',num2str(verrconc(iter))])
text(10,41,['Epsil_WC_RefT = ',num2str(Epsil_WC_RefT),' mds = ',num2str(mds)],...
    'Interpreter','none');
text(10,32,['Hypochromicity = ',...
    num2str(100*(1-Epsil_WC_RefT/(Epsil_ss1_RefT + Epsil_ss2_RefT))),' %']);
text(10,23,['Data analysis run on ',date,'; iteration',num2str(iter)]);
text(10,14,[filestem,'_n',num2str(iter)],'Interpreter','none');
xticklabels({});
yticklabels({});
axis([0 100 0 100])
if (scale_flag == 1)
    pfile = [filestem '_n' num2str(iter) '_sf' num2str(scale_factor) '.png'];
else
    pfile = [filestem '_n' num2str(iter) '.png'];
end
export_fig(pfile,'-m2');
    end
end
%
% Section for called functions
% Some are in separate files for use in other routines
%

function dsbaseline = DSBpred(T,RefT,Cr_tot,mds,Epsil_WC_RefT)
%TminK = min(TK);
%TmaxK = max(TK);
dsbaseline = (Epsil_WC_RefT+(T - RefT)*mds).*(Cr_tot);
end

% function ssbaseline = SSBpred(x,TK,CT,ExW,mss,mds,Amin,Amax)
% %TminK = min(TK);
% TmaxK = max(TK);
% ssbaseline = 0.5*(mss.*(TK-TmaxK)+2*Amax);
% end

function fdsexp = Fdsexpt(Abs_expt,T,RefT,Abs_ds_slope,Abs_RefT,Abs_ss_ref_tot)
% TminK = min(TK);
% TmaxK = max(TK);
%Trange = TmaxK - TminK;
%Abs_ss = 0.5*(mss.*(TK-TmaxK)+2*Amax);
Abs_ds = Abs_RefT + Abs_ds_slope.*(T-RefT);
fdsexp = (Abs_expt - Abs_ss_ref_tot)./(Abs_ds - Abs_ss_ref_tot);
end

