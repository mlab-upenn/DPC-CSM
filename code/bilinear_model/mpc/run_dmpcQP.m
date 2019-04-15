% -------------------------------------------------------------------------
% Function BC_DMPC_11 is used for HVAC02 naive MPC approach 
% -------------------------------------------------------------------------
% 
% - Optimization problem formulated in standard LP format
% - simulate whole year of 2007, by simulating winter-summer-winter, 
%   i.e. Jan-April, May-Oct,Nov-Dec
% - 22 <= Troom_sum <= 26, 21 <= Troom_win <= 25
% - uses RGS computed from RG, with and without KF
% - optimization with only first two constraints
% - disturbance can be chosen to be either rawPred, KF Pred, AR2 Pred
% - fixed statistics part
% - slack variables are penalized linearly and individually
% - prediction data using KF(measurement error) und AR4 (improve predictions)
% - also, need KF to estimate prediction error due to measurement noise
% - internal gain prediction uses mean of past data
% - have instantanous correction for light if it's too dark
% - inputs u3 and u4 are removed: no heating and cooling slabs


% Author: Xiaojing Zhang
% based on works of: Frauke Oldewurtel, David Sturzenegger
% Date: 14. Aug 2012
% use discrete-time building models from BACLab software v.120

% for HVAC01 buildings:
%   States:
%     x(1)      Troom      : room temperature                                  [degC]
%     x(2:6)    Tslab(i)   : slab temperatures 1...5                           [degC]
%     x(7:9)    TiWall(i)  : inner wall temperatures 1...3                     [degC]
%     x(10:12)  ToWall(i)  : outside wall temperatures 1...3                   [degC]

%   Control inputs:
%     u(1)      bPos       : blind position [0: closed ... 1: open]            [-]
%     u(2)      eLighting  : gains electric lighting                           [W/m2]
%     -u(3)      hPowSlab   : heating power slab, positive values = heating     [W/m2]
%     -u(4)      cPowSlab   : cooling power slab, positive values = cooling     [W/m2]
%     u(5)      fcUsgFact  : free cooling usage factor [0: off ... 1: on]      [-]
%     u(6)      hPowRad    : heating power radiator, positive values = heating [W/m2]

%   Disturbance inputs:
%     v(1)      solG      : solar gains with fully closed blinds              [W/m2]
%     v(2)      dSolG     : additional solar gains with open blinds           [W/m2]
%     v(3)      illum     : daylight illuminance with fully closed blinds     [lux]
%     v(4)      dIllum    : additional daylight illuminance with open blinds  [lux]
%     v(5)      persG     : internal gains (persons)                          [W/m2]
%     v(6)      equipG    : internal gains (equipment)                        [W/m2]
%     v(7)      Tair      : outside air temperature                           [degC]
%     v(8)      TfreeCool : free cooling temperature                          [degC]

%   Outputs:
%     y(1)      Troom      : room temperature                                  [degC]
%     y(2)      Eroom      : room illuminance                                  [lux]
%     y(3)      TsrfCeil   : ceiling surface temperature                       [degC] 
% clear all
% close all
clc
%    addpath(genpath('C:/Program Files/IBM')) 

    % define list of simulations to be done, be careful with winter!
    % format: weather, setback, pa-h-wh-ih, numSample
    ctrlHzn = 6;
    list = {'MSM2007',1, 'pa', 'h', 'wh', 'ih' , ctrlHzn , 0};    

        sim.constrTight = list{8};  % constraint tightening
        sim.predhor = list{7};       % 7 days for PB = 168
        build.weather = list{1};           % 'MSM2007', 'WHW2007'
        build.KF = 'KF_AR4_';    % '', 'KF_', 'AR2_'
        build.setback = list{2};      % 1: use ymin_sb; 0: ymin_sb = ymin
        
        % setup building information
        build.sys = 'e01';                  % {e01}
        build.type = list{3};           % {Swiss Average, PAssive house}
        build.env = list{4};            % {Heavy, Light}
        build.win = list{5};            % {Window Low, Window High}
        build.int = list{6};            % {Internal gain High, Internal gain Low}
        build.dir = 'S';                    % {North, South, South East, South West}
        build.vac = 1;                      % code for vacancy
        build.occ = 5;                      % code for occupancy

        % setup simulatin environment

        total_steps = 8768-8-168;
        % Needed weather variables [Tair, TFreeCool, SolRadN, SolRadE, SolRadS, SolRadW]
        sim.wVarNeeded =  [ 1 1 0 0 1 0 ];  % must be consistent with "build.dir"
        sim.igVarNeeded = [ 1 1 ];  % Needed internal gains [persG equipG]
        cFact = 1.0000e-03;     % Cost scaling factor, from kW to W
        sim.Tolc = 1;               % number of control inputs applied in open-loop fashion
        sim.mc = 1;
        [bm] = get_building_model(build);
        xref = 21.5*ones(size(bm.A,1),1);

        % -------------------------------------------------------
        %  1. First run simulation: Jan 01 till April 30;
        % --------------------------------------------------------

        build.seas = 'win';                 % winter season
        sim.simhor = 200;      % simulate Jan 01 - April 30
        sim.xinit = [ 22 22 22 22 22 22 22 22 22 20 15 10 ]'; % start for winter
        [ucost, ymin,ymin_sb, ymax,ymax_sb, umin, umax] = getCostsAndConstraints(build); 
        build_stand = build;
        build_stand.vac = NaN;      % use simple occupancy control
        % vPred_all is cell
%         [v vPred_all] = get_disturbance(build_stand, sim, 1);  % get disturbance for Jan-April
        load v
        load vPred_all
        vPred_all = vPred_all(1:sim.simhor);
        for ii = 1:sim.simhor
           vPred_all{ii} = vPred_all{ii}(1:ctrlHzn,:); 
        end
        v = v(1:sim.simhor,:);
        [ymin_all, ymax_all] = get_constraints(sim, v, ymin, ymax, ymin_sb, ymax_sb, 1);
        ymin_all(end-sim.predhor+1:end,1) = 22*ones(sim.predhor,1);   % correct for summer
        ymax_all(end-sim.predhor+1:end,1) = 26*ones(sim.predhor,1);   % correct for summer
        ucost_pred = repmat(ucost,sim.predhor,1);
        umin_pred = repmat(umin,sim.predhor,1);
        umax_pred = repmat(umax,sim.predhor,1);
        x = sim.xinit;
        state = nan(sim.simhor+1,size(bm.A,1));
        state(1,:) = x';
        output = nan(sim.simhor,size(bm.C,1));
        input = nan(size(bm.Bu,2),sim.simhor);

        for k = 1 : sim.Tolc : sim.simhor
            disp([num2str(k),' / ',num2str(sim.simhor)])
            vPred = vPred_all{k};
            ymin_pred = ymin_all(k:k+sim.predhor-1,:);
            ymax_pred = ymax_all (k:k+sim.predhor-1,:);

%             u = get_control_inputQP(x, bm, sim, v(k,:)', vPred, ymin_pred, ymax_pred, ucost_pred, umin_pred, umax_pred, xref);
%             u = get_control_input(x, bm, sim, v(k,:)', vPred, ymin_pred, ymax_pred, ucost_pred, umin_pred, umax_pred);
            u = FakeBilinearMPC(x, bm, sim, v(k,:)', vPred, ymin_pred, ymax_pred, ucost_pred, umin_pred, umax_pred);
            [x_k, output_k] = state_update(build,bm,sim,x,v,u,k);

            % collect states, inputs, etc.
            input(:,k:k+sim.Tolc-1) = u(1:sim.Tolc,:)';
            state(k+1:k+1+sim.Tolc-1,:) = x_k';
            x = x_k(:,1);
            output(k:k+sim.Tolc-1,:) = output_k';
        end

        % store temporarily
        ymin_coll = ymin_all(1:k,:);
        ymax_coll = ymax_all(1:k,:);
        state_coll = state(1:end-1,:);
        input_coll = input;
        output_coll = output;
        numSimSteps = k;

      
%         -------------------------------------------------------
%          2. Run for summer = May 01 - Oct 31
%         --------------------------------------------------------

        build.seas = 'sum';                 % winter season
        sim.simhor = 7303-2888+1;      % simulate May 01 - Okt 31
        [ucost, ymin,ymin_sb, ymax,ymax_sb, umin, umax] = getCostsAndConstraints(build); 
        build_stand = build;
        build_stand.vac = NaN;      % use simple occupancy control
        [v vPred_all] = get_disturbance(build_stand, sim, 2888-7);  % get disturbance for May-Oct
        [ymin_all, ymax_all] = get_constraints(sim, v, ymin, ymax, ymin_sb, ymax_sb, 2888-7);
        ymin_all(end-sim.predhor+1:end,1) = 21*ones(sim.predhor,1);   % correct for winter
        ymax_all(end-sim.predhor+1:end,1) = 25*ones(sim.predhor,1);   % correct for winter
        ucost_pred = repmat(ucost,sim.predhor,1);
        umin_pred = repmat(umin,sim.predhor,1);
        umax_pred = repmat(umax,sim.predhor,1);
        state = nan(sim.simhor+1,size(bm.A,1));
        state(1,:) = x';
        output = nan(sim.simhor,size(bm.C,1));
        input = nan(size(bm.Bu,2),sim.simhor);

        for k = 1 : sim.Tolc : sim.simhor
            disp([num2str(k),' / ',num2str(sim.simhor)])
            vPred = vPred_all{k};
            ymin_pred = ymin_all(k:k+sim.predhor-1,:);
            ymax_pred = ymax_all (k:k+sim.predhor-1,:);

            u = get_control_input(x, bm, sim, v(k,:)', vPred, ymin_pred, ymax_pred, ucost_pred, umin_pred, umax_pred);

            [x_k, output_k] = state_update(build,bm,sim,x,v,u,k);
            
%           collect states, inputs, etc.
            input(:,k:k+sim.Tolc-1) = u(1:sim.Tolc,:)';
            state(k+1:k+1+sim.Tolc-1,:) = x_k';
            x = x_k(:,1);
            output(k:k+sim.Tolc-1,:) = output_k';
        end

        %store temporarily
        ymin_coll = [ymin_coll ; ymin_all(1:k,:)];
        ymax_coll = [ymax_coll ; ymax_all(1:k,:)];
        state_coll = [state_coll ; state(1:end-1,:)];
        input_coll = [input_coll input];
        output_coll = [output_coll ; output];
        numSimSteps = numSimSteps + k;

        % -------------------------------------------------------
        %  3. Run for summer = Nov 01 - Dec 31
        % --------------------------------------------------------

        build.seas = 'win';                 % winter season
        sim.simhor = 8767-7304-168+1;      % simulate May 01 - Okt 31
        [ucost, ymin,ymin_sb, ymax,ymax_sb, umin, umax] = getCostsAndConstraints(build); 
        build_stand = build;
        build_stand.vac = NaN;      % use simple occupancy control
        [v vPred_all] = get_disturbance(build_stand, sim, 7304-7);  % get disturbance for May-Oct
        [ymin_all, ymax_all] = get_constraints(sim, v, ymin, ymax, ymin_sb, ymax_sb, 7304-7);
        ucost_pred = repmat(ucost,sim.predhor,1);
        umin_pred = repmat(umin,sim.predhor,1);
        umax_pred = repmat(umax,sim.predhor,1);
        state = nan(sim.simhor+1,size(bm.A,1));
        state(1,:) = x';
        output = nan(sim.simhor,size(bm.C,1));
        input = nan(size(bm.Bu,2),sim.simhor);

        for k = 1 : sim.Tolc : sim.simhor
            disp([num2str(k),' / ',num2str(sim.simhor)])
            vPred = vPred_all{k};
            ymin_pred = ymin_all(k:k+sim.predhor-1,:);
            ymax_pred = ymax_all (k:k+sim.predhor-1,:);

            u = get_control_input(x, bm, sim, v(k,:)', vPred, ymin_pred, ymax_pred, ucost_pred, umin_pred, umax_pred);

            [x_k, output_k] = state_update(build,bm,sim,x,v,u,k);
            
            % collect states, inputs, etc.
            input(:,k:k+sim.Tolc-1) = u(1:sim.Tolc,:)';
            state(k+1:k+1+sim.Tolc-1,:) = x_k';
            x = x_k(:,1);
            output(k:k+sim.Tolc-1,:) = output_k';
        end

        % store temporarily
        ymin_coll = [ymin_coll ; ymin_all(1:k,:)];
        ymax_coll = [ymax_coll ; ymax_all(1:k,:)];
        state_coll = [state_coll ; state(1:end-1,:)];
        input_coll = [input_coll input];
        output_coll = [output_coll ; output];
        numSimSteps = numSimSteps + k;

        clc
        disp('============ bye! ============')
        close all
            
%             % create predictions / realizations
%             % occupancy
%             cd('InData/occupancy')
%             load('ig_act_real.mat')
%             load('ig_act_pred_1we.mat');
%             cd ..; cd ..;
%         
%             figure(); hold on;
%             title('occupancy information')
%             plot(repmat(ig_act_pred_1we,[52,1]),'b','LineWidth',2)
%             plot(ig_act_real,'--r','LineWidth',2)
%             legend('IG Prediction','IG Realization')
%             h_ig = gca;
%             hold off
%            
%             cd(['InData/weatherData/' build.weather]);
%             % get new data from 2007 measurements
%             load(['Vreal_' build.weather '.mat']);      % load data from weather realization in 2007
%             load(['TAreal_' build.weather '.mat']);
%             load(['TWreal_' build.weather '.mat']);
%             load(['dataTAP_' build.KF num2str(sim.predhor) '_' build.weather '.mat']);      % load data from weather realization in 2007
%             load(['dataTWP_' build.KF num2str(sim.predhor) '_' build.weather '.mat']);      % load data from weather realization in 2007
%             load(['dataRGSP_' build.KF num2str(sim.predhor) '_' build.weather '.mat']);      % load data from weather realization in 2007
%             cd ..
%             cd ..
%             cd ..
%             
%             figure(); hold on;
%             title('Solar radiation')
%             plot(V,'b')
%             plot(dataRGSP_KF_AR4(:,1),'r')
%             legend('realization','prediction 1st step')
%             h_RGS = gca;
%             hold off;
%             
%             figure(); hold on;
%             title('TA')
%             plot(TA,'b')
%             plot(dataTAP_KF_AR4(:,1),'r')
%             legend('realization','prediction 1st step')
%             h_TA = gca;
%             hold off;
% 
%             figure(); hold on;
%             title('TW')
%             plot(TW,'b')
%             plot(dataTWP_KF_AR4(:,1),'r')
%             legend('realization','prediction 1st step')
%             h_TW = gca;
%             hold off;
%             
%             linkaxes([h_troom; h_light; h_inputs; h_costs; h_ig; h_RGS; h_TA; h_TW],'x')