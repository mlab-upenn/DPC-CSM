% ----------------------------------------------------------------------------------
% Function BC_PB_09 is used for deterministic occupancy investigations 
% ----------------------------------------------------------------------------------

% - Perfekt prediction of weather data taken at Basel, 2007
% - Optimization problem formulated in standard LP format
% - two slack variables Smax, Smin on Y() to admit feasibility, penalized
%   linearly
% - data from solar radiation are obtained from BAS 2007 measurements
% - simulate whole year of 2007, by simulating winter-summer-winter, 
%   i.e. Jan-April, May-Oct,Nov-Dec
% - basically same as BC_PB_06, some parts adapted for Remote Control
% - choose weather data according to collected data
% - 22 <= Troom_sum <= 26, 21 <= Troom_sum <= 25
% - internal gain from Actelion
% - no u3 und u4 in optimization
% - slack variables individually penalized

% Author: Xiaojing Zhang
% based on works of: Frauke Oldewurtel, David Sturzenegger
% Date: 24. Mai 2012
% use discrete-time building models from BACLab software v.120

% for e01 buildings:
%   States:
%     x(1)      Troom      : room temperature                                  [degC]
%     x(2:6)    Tslab(i)   : slab temperatures 1...5                           [degC]
%     x(7:9)    TiWall(i)  : inner wall temperatures 1...3                     [degC]
%     x(10:12)  ToWall(i)  : outside wall temperatures 1...3                   [degC]

%   Control inputs:
%     u(1)      bPos       : blind position [0: closed ... 1: open]            [-]
%     u(2)      eLighting  : gains electric lighting                           [W/m2]
%     u(3)      hPowSlab   : heating power slab, positive values = heating     [W/m2]
%     u(4)      cPowSlab   : cooling power slab, positive values = cooling     [W/m2]
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


%% main function
function BC_PB_09

    serv = 0;       % if 1, set license for server, and no plots
                        % set 0 if used locally
    % define list of simulations to be done
    % format: weather, setback, pa-h-wh-ih
    list = { 
        {'WHW2007',1, 'pa', 'h', 'wh', 'ih' }, ...
        {'WHW2007',1, 'sa', 'h', 'wl', 'ih' }, ...
%         {'MSM2007' ,1, 'pa', 'h', 'wh', 'ih' }, ...
%         {'MSM2007',1, 'sa', 'h', 'wl', 'ih' }
        };

    %%
    % CPLEX license
    if ~serv
        setenv('ILOG_LICENSE_FILE','/Users/homePC/Documents/MATLAB/cplexlic/access.ilm');
    elseif serv==1
        setenv('ILOG_LICENSE_FILE','/home/xiaozhan/Desktop/MasterArbeit/MATLAB/cplexlic/access.ilm');
        addpath('/home/xiaozhan/Desktop/MasterArbeit/MATLAB/cplex12_unix/')
    elseif serv==2
        setenv('ILOG_LICENSE_FILE','C:\Users\xiaozhan\Documents\MATLAB\cplex_lic\access.ilm');
    end
    
    clear cplexint
    
    for jj = 1 : length(list)
        
        clearvars -except jj list serv;
        clc
        if ~serv
            close all;
        end
        
        t_bigbang = tic;

        
        sim.predhor = 168;       % 7 days for PB = 168
        build.weather = list{jj}{1};           % 'BAS2007', 'SMA2007'
        build.setback = list{jj}{2};

        % setup building information
        build.sys = 'e01';                 % {e01, e05}
        build.type = list{jj}{3};          % {Swiss Average, PAssive house}
        build.env = list{jj}{4};           % {Heavy, Light}
        build.win = list{jj}{5};           % {Window Low, Window High}
        build.int = list{jj}{6};           % {Internal gain High, Internal gain Low}
        build.dir = 'S';                   % {North, South, South East, South West}
        build.vac = 1;                     % code for vacancy
        build.occ = 5;                     % code for occupancy

        % setup simulatin environment
        total_steps = 8768-8-sim.predhor;
        % Needed weather variables [Tair, TFreeCool, SolRadN, SolRadE, SolRadS, SolRadW]
        sim.wVarNeeded =  [ 1 1 0 0 1 0 ];  % must be consistent with "build.dir"
        sim.igVarNeeded = [ 1 1 ];  % Needed internal gains [persG equipG]
        cFact = 1.0000e-03;     % Cost scaling factor, from kW to W
        sim.Tolc = 1;               % number of control inputs applied in open-loop fashion
        sim.mc = 1;
        [bm] = get_building_model(build);

        %% -------------------------------------------------------
        %  1. First run simulation: Jan 01 till April 30;
        % --------------------------------------------------------

        build.seas = 'win';                 % winter season
        sim.simhor = 2880;      % simulate Jan 01 - April 30
        sim.xinit = [ 22 22 22 22 22 22 22 22 22 20 15 10 ]'; % start for winter
        [ucost, ymin,ymin_sb, ymax,ymax_sb, umin, umax] = getCostsAndConstraints(build); 
        build_stand = build;
        build_stand.vac = NaN;      % use simple occupancy control
        [v] = get_disturbance(build_stand, sim, 1);  % get disturbance for Jan-April
        vPred_all = v;
        [ymin_all, ymax_all] = get_constraints(sim, vPred_all, ymin, ymax, ymin_sb, ymax_sb, 1);
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
        t_start = tic;
        t_start_thisrun = tic;
        control_cost = 0;
        control_cost_all = nan(sim.simhor,1);
        numConstrViol = 0;
        constrViol = 0;     % in [Kelvin Hours]
        numConstrViolY1 = 0;
        numConstrViolY2 = 0;
        numConstrViolY3 = 0;
        numSoftConstrViol = 0;
        numSimSteps = 0;

        for k = 1 : sim.Tolc : sim.simhor
            if (k>1)
               t_delta = toc(t_start);
               t_delta_thisrun = toc(t_start_thisrun);
               disp(['-- total elapsed time: ', num2str(t_delta), 'sec, k: ', num2str(k), ', %done of Jan-April: ', num2str(100*k/total_steps)]);
            end
            vPred = vPred_all(k:k+sim.predhor-1,:);
            ymin_pred = ymin_all(k:k+sim.predhor-1,:);
            ymax_pred = ymax_all (k:k+sim.predhor-1,:);

            [u u_all softConstr] = get_control_input(x, bm, sim, vPred, ymin_pred, ymax_pred, ucost_pred, umin_pred, umax_pred);

            if (softConstr>0)
                disp(['soft Constraints violated by: ' num2str(softConstr)]);
                numSoftConstrViol = numSoftConstrViol+1;
            end

            
            %% for debugging
%         if k >= 12
%             sim.Tolc = 168;
%             [x_k, output_k] = state_update(build,bm,sim,x,v,u,k);
% 
%             close all;
%             figure(); hold on;
%             hold on; 
%             title(['BC\_PB\_04: ' build.weather ', pred\_hor = ' num2str(sim.predhor)  ', total cost: ' num2str(round(control_cost))]);
%             plot(output_k(1,:)); plot(ymin_pred(:,1),'r'); plot(ymax_pred(:,1),'m');
%             legend('y_1 = Troom')
%             ylim([15 30])
%             hold off
% 
%             figure();
%             hold on; 
%             title(['BC\_PB\_04: ' build.weather  ', pred\_hor = ' num2str(sim.predhor)  ', total cost: ' num2str(round(control_cost))]);
%             plot(output_k(2,:)); plot(ymin_pred(:,2),'--r'); %plot(ymax_all(:,2),'--m');
%             legend('y_2 = Eroom')
%             hold off
% 
%             figure(); 
%             hold on; 
%             title(['BC\_PB\_s06: ' build.weather ', pred\_hor = ' num2str(sim.predhor)  ', total cost: ' num2str(round(control_cost))]);
%             plot(1:sim.Tolc,u(:,1), 1:sim.Tolc,u(:,2), 1:sim.Tolc,u(:,3), ...
%                  1:sim.Tolc,u(:,4) );
%             legend('bPos','eLighting','fcUsgFact','hPowRad');
%             hold off
% 
%             keyboard
%         end
            
            
            %%
            [x_k, output_k] = state_update(build,bm,sim,x,v,u,k);
            disp(['applied input: ' num2str(u(1,:))]);
            control_cost = control_cost + sum(diag(ucost_pred(1:sim.Tolc,:)*u(1:sim.Tolc,:)'));
            control_cost_all(k) = control_cost;
            disp(['accumulated cost: ' num2str(control_cost)]);

    %         if (k == 20)
    %            return 
    %         end

            % check for constraint violoations
            if (ymax_all(k:k+sim.Tolc-1,1)'  + 1e-10 - output_k(1)) < 0
                numConstrViolY1 = numConstrViolY1 + 1;
                constrViol = constrViol + norm(ymax_all(k:k+sim.Tolc-1,1)' - output_k(1));
            end
            if (ymax_all(k:k+sim.Tolc-1,2)'  + 1e-10 - output_k(2)) < 0
                numConstrViolY2 = numConstrViolY2 + 1;
            end
            if (ymax_all(k:k+sim.Tolc-1,3)'  + 1e-10 - output_k(3)) < 0
                numConstrViolY3 = numConstrViolY3 + 1;
            end
            if (output_k(1)  + 1e-10- ymin_all(k:k+sim.Tolc-1,1)') < 0
                numConstrViolY1 = numConstrViolY1 + 1;
                constrViol = constrViol + norm(ymin_all(k:k+sim.Tolc-1,1)' - output_k(1));
            end
            if (output_k(2)  + 1e-10 - ymin_all(k:k+sim.Tolc-1,2)') < 0
                numConstrViolY2 = numConstrViolY2 + 1;
            end
            if (output_k(3) + 1e-10 - ymin_all(k:k+sim.Tolc-1,3)') < 0
                numConstrViolY3 = numConstrViolY3 + 1;
            end

            if sum((ymax_all(k:k+sim.Tolc-1,:)' - output_k + 1e-10) >= zeros(3,sim.Tolc)) ~= 3
                numConstrViol = numConstrViol+1;
                disp('max constraints violated')
            elseif sum( (output_k - ymin_all(k:k+sim.Tolc-1,:)' + 1e-10 ) >= zeros(3,sim.Tolc)) ~= 3
                numConstrViol = numConstrViol+1;
                disp('min constraints violated')
            end

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
        control_cost_coll = control_cost_all;
        numSimSteps = k;

%         keyboard;

        %% -------------------------------------------------------
        %  2. Run for summer = May 01 - Oct 31
        % --------------------------------------------------------

        build.seas = 'sum';                 % winter season
        sim.simhor = 7303-2888+1;      % simulate May 01 - Okt 31
        [ucost, ymin,ymin_sb, ymax,ymax_sb, umin, umax] = getCostsAndConstraints(build); 
        build_stand = build;
        build_stand.vac = NaN;      % use simple occupancy control
        [v] = get_disturbance(build_stand, sim, 2888-7);  % get disturbance for May-Oct
        vPred_all = v;
        [ymin_all, ymax_all] = get_constraints(sim, vPred_all, ymin, ymax, ymin_sb, ymax_sb, 2888-7);
        ymin_all(end-sim.predhor+1:end,1) = 21*ones(sim.predhor,1);   % correct for winter
        ymax_all(end-sim.predhor+1:end,1) = 25*ones(sim.predhor,1);   % correct for winter
        ucost_pred = repmat(ucost,sim.predhor,1);
        umin_pred = repmat(umin,sim.predhor,1);
        umax_pred = repmat(umax,sim.predhor,1);
        state = nan(sim.simhor+1,size(bm.A,1));
        state(1,:) = x';
        output = nan(sim.simhor,size(bm.C,1));
        control_cost_all = nan(sim.simhor,1);
            input = nan(size(bm.Bu,2),sim.simhor);
        t_start = tic;
        t_start_thisrun = tic;

        for k = 1 : sim.Tolc : sim.simhor
            if (k>1)
               t_delta = toc(t_start);
               t_delta_thisrun = toc(t_start_thisrun);

               disp(['-- total elapsed time: ', num2str(t_delta), 'sec, k: ', num2str(k), ', %done of May-Oct: ', num2str(100*(k+numSimSteps)/total_steps)]);
            end
            vPred = vPred_all(k:k+sim.predhor-1,:);
            ymin_pred = ymin_all(k:k+sim.predhor-1,:);
            ymax_pred = ymax_all (k:k+sim.predhor-1,:);

            [u u_all softConstr] = get_control_input(x, bm, sim, vPred, ymin_pred, ymax_pred, ucost_pred, umin_pred, umax_pred);

%             if (k>1830)
%                 sim.Tolc = 168;
%                 [x_k, output_k] = state_update(build,bm,sim,x,v,u,k);
% 
%                 close all;
%                 figure(); hold on;
%                 hold on; 
%                 title(['BC\_PB\_04: ' build.weather ', pred\_hor = ' num2str(sim.predhor)  ', total cost: ' num2str(round(control_cost))]);
%                 plot(output_k(1,:)); plot(ymin_pred(:,1),'r'); plot(ymax_pred(:,1),'m');
%                 legend('y_1 = Troom')
%                 ylim([15 30])
%                 hold off
% 
%                 figure();
%                 hold on; 
%                 title(['BC\_PB\_04: ' build.weather  ', pred\_hor = ' num2str(sim.predhor)  ', total cost: ' num2str(round(control_cost))]);
%                 plot(output_k(2,:)); plot(ymin_pred(:,2),'--r'); %plot(ymax_all(:,2),'--m');
%                 legend('y_2 = Eroom')
%                 hold off
% 
%                 figure(); 
%                 hold on; 
%                 title(['BC\_PB\_s06: ' build.weather ', pred\_hor = ' num2str(sim.predhor)  ', total cost: ' num2str(round(control_cost))]);
%                 plot(1:sim.Tolc,u(:,1), 1:sim.Tolc,u(:,2), 1:sim.Tolc,u(:,3), ...
%                      1:sim.Tolc,u(:,4));
%                 legend('bPos','eLighting','fcUsgFact','hPowRad');
%                 hold off
% 
%                 keyboard
%             end
            
            
            if (softConstr>0)
                disp(['soft Constraints violated by: ' num2str(softConstr)]);
                numSoftConstrViol = numSoftConstrViol+1;
            end

            [x_k, output_k] = state_update(build,bm,sim,x,v,u,k);
            disp(['applied input: ' num2str(u(1,:))]);
            control_cost = control_cost + sum(diag(ucost_pred(1:sim.Tolc,:)*u(1:sim.Tolc,:)'));
            control_cost_all(k) = control_cost;
            disp(['accumulated cost: ' num2str(control_cost)]);

            % check for constraint violations
            if (ymax_all(k:k+sim.Tolc-1,1)' + 1e-10 - output_k(1)) < 0
                numConstrViolY1 = numConstrViolY1 + 1;
                constrViol = constrViol + norm(ymax_all(k:k+sim.Tolc-1,1)' - output_k(1));
            end
            if (ymax_all(k:k+sim.Tolc-1,2)' + 1e-10 - output_k(2)) < 0
                numConstrViolY2 = numConstrViolY2 + 1;
            end
            if (ymax_all(k:k+sim.Tolc-1,3)' + 1e-10 - output_k(3)) < 0
                numConstrViolY3 = numConstrViolY3 + 1;
            end
            if (output_k(1) + 1e-10 - ymin_all(k:k+sim.Tolc-1,1)') < 0
                numConstrViolY1 = numConstrViolY1 + 1;
                constrViol = constrViol + norm(ymin_all(k:k+sim.Tolc-1,1)' - output_k(1));
            end
            if (output_k(2) + 1e-10 - ymin_all(k:k+sim.Tolc-1,2)') < 0
                numConstrViolY2 = numConstrViolY2 + 1;
            end
            if (output_k(3) + 1e-10 - ymin_all(k:k+sim.Tolc-1,3)') < 0
                numConstrViolY3 = numConstrViolY3 + 1;
            end

            if sum((ymax_all(k:k+sim.Tolc-1,:)' - output_k + 1e-10) >= zeros(3,sim.Tolc)) ~= 3
                numConstrViol = numConstrViol+1;
                disp('max constraints violated')
            elseif sum( (output_k - ymin_all(k:k+sim.Tolc-1,:)' + 1e-10 ) >= zeros(3,sim.Tolc)) ~= 3
                numConstrViol = numConstrViol+1;
                disp('min constraints violated')
            end

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
        control_cost_coll = [control_cost_coll ; control_cost_all];
        numSimSteps = numSimSteps + k;

        %% -------------------------------------------------------
        %  3. Run for summer = Nov 01 - Dec 31
        % --------------------------------------------------------

        build.seas = 'win';                 % winter season
        sim.simhor = 8767-7304-168+1;      % simulate May 01 - Okt 31
        [ucost, ymin,ymin_sb, ymax,ymax_sb, umin, umax] = getCostsAndConstraints(build); 
        build_stand = build;
        build_stand.vac = NaN;      % use simple occupancy control
        [v] = get_disturbance(build_stand, sim, 7304-7);  % get disturbance for May-Oct
        vPred_all = v;
        [ymin_all, ymax_all] = get_constraints(sim, vPred_all, ymin, ymax, ymin_sb, ymax_sb, 7304-7);
        ucost_pred = repmat(ucost,sim.predhor,1);
        umin_pred = repmat(umin,sim.predhor,1);
        umax_pred = repmat(umax,sim.predhor,1);
        state = nan(sim.simhor+1,size(bm.A,1));
        state(1,:) = x';
        output = nan(sim.simhor,size(bm.C,1));
        input = nan(size(bm.Bu,2),sim.simhor);
        control_cost_all = nan(sim.simhor,1);
        t_start = tic;
        t_start_thisrun = tic;

        for k = 1 : sim.Tolc : sim.simhor
            if (k>1)
               t_delta = toc(t_start);
               t_delta_thisrun = toc(t_start_thisrun);

               disp(['-- total elapsed time: ', num2str(t_delta), 'sec, k: ', num2str(k), ', %done of Nov-December: ', num2str(100*(k+numSimSteps)/total_steps)]);
            end
            vPred = vPred_all(k:k+sim.predhor-1,:);
            ymin_pred = ymin_all(k:k+sim.predhor-1,:);
            ymax_pred = ymax_all (k:k+sim.predhor-1,:);

            [u u_all softConstr] = get_control_input(x, bm, sim, vPred, ymin_pred, ymax_pred, ucost_pred, umin_pred, umax_pred);

            if (softConstr>0)
                disp(['soft Constraints violated by: ' num2str(softConstr)]);
                numSoftConstrViol = numSoftConstrViol+1;
            end

            [x_k, output_k] = state_update(build,bm,sim,x,v,u,k);
            disp(['applied input: ' num2str(u(1,:))]);
            control_cost = control_cost + sum(diag(ucost_pred(1:sim.Tolc,:)*u(1:sim.Tolc,:)'));
            control_cost_all(k) = control_cost;
            disp(['accumulated cost: ' num2str(control_cost)]);

            % check for constraint violoations
            if (ymax_all(k:k+sim.Tolc-1,1)' + 1e-10 - output_k(1)) < 0
                numConstrViolY1 = numConstrViolY1 + 1;
                constrViol = constrViol + norm(ymax_all(k:k+sim.Tolc-1,1)' - output_k(1));
            end
            if (ymax_all(k:k+sim.Tolc-1,2)' + 1e-10 - output_k(2)) < 0
                numConstrViolY2 = numConstrViolY2 + 1;
            end
            if (ymax_all(k:k+sim.Tolc-1,3)' + 1e-10 - output_k(3)) < 0
                numConstrViolY3 = numConstrViolY3 + 1;
            end
            if (output_k(1) + 1e-10 - ymin_all(k:k+sim.Tolc-1,1)') < 0
                numConstrViolY1 = numConstrViolY1 + 1;
                constrViol = constrViol + norm(ymin_all(k:k+sim.Tolc-1,1)' - output_k(1));
            end
            if (output_k(2) + 1e-10 - ymin_all(k:k+sim.Tolc-1,2)') < 0
                numConstrViolY2 = numConstrViolY2 + 1;
            end
            if (output_k(3) + 1e-10 - ymin_all(k:k+sim.Tolc-1,3)') < 0
                numConstrViolY3 = numConstrViolY3 + 1;
            end

            if sum((ymax_all(k:k+sim.Tolc-1,:)' - output_k + 1e-10) >= zeros(3,sim.Tolc)) ~= 3
                numConstrViol = numConstrViol+1;
                disp('max constraints violated')
            elseif sum( (output_k - ymin_all(k:k+sim.Tolc-1,:)' + 1e-10 ) >= zeros(3,sim.Tolc)) ~= 3
                numConstrViol = numConstrViol+1;
                disp('min constraints violated')
            end

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
        control_cost_coll = [control_cost_coll ; control_cost_all]';

        
        time = toc(t_bigbang)/60;
        %%
        clc
        disp('============= Statistics PB_09 ===============')
%         disp(['Number of soft constraints violated: ' num2str(numSoftConstrViol)]);
%         disp(['Number of 1st constraints violated: ' num2str(numConstrViolY1)]);
        disp(['Relative Number of 1st constraints violated: ' num2str(numConstrViolY1/total_steps*100) '%']);
        disp(['Violation in Kelvin Hours (Kh): ' num2str(constrViol)]);
        disp(['Number of 2nd constraints violated (robust): ' num2str(numConstrViolY2/total_steps*100) '%']);
%         disp(['Number of 2nd constraints violated (robust): ' num2str(sum(output_coll(:,2)<500-1e-10)/total_steps*100) '%']);
%         disp(['Number of 3rd constraints violated: ' num2str(numConstrViolY3)]);
%         disp(['Number of constraints violated: ' num2str(numConstrViol)]);
        disp(['Relative Number of  constraints violated: ' num2str(numConstrViol/numSimSteps*100) '%']);
%         disp(['cost input: ' num2str(control_cost)]);
        disp(['total time elapsed: ' num2str(time) ' min'])
        
        disp('       ')
        % post processing of individual cost
        blindsCost = sum(0.01*input_coll(1,2881:7296));
        lightCost = sum(3.32*input_coll(2,:));
        fcUsgFactCost = sum(7.47*input_coll(3,:));
        hPowRadCost = sum(1.107*input_coll(4,:));
        totalCost = blindsCost + lightCost + fcUsgFactCost + hPowRadCost;
        HVACCost = fcUsgFactCost + hPowRadCost;
        disp(['blind cost: ' num2str(blindsCost)])
        disp(['light cost: ' num2str(lightCost)])
        disp(['fcUsgFact cost: ' num2str(fcUsgFactCost)])
        disp(['hPowRad cost: ' num2str(hPowRadCost)])
        disp(['HVAC cost: ' num2str(HVACCost) ])
        disp(['total cost 1: ' num2str(totalCost) ])
        disp(['total cost 2: ' num2str(control_cost) ])
        
        max_violations = find(output_coll(:,1) - 1e-10 > ymax_coll(:,1));
        min_violations = find(output_coll(:,1) < ymin_coll(:,1) - 1e-10);
        disp(['max constr violated: ' num2str(max_violations')]);
        disp(['min constr violated: ' num2str(min_violations')]);
        
        % collect cumulative individual input costs
        control_cost_ind_coll = nan(size(input_coll));
        control_cost_ind_coll(:,1) = ucost' .* input_coll(:,1);
        for ii = 2 : size(input_coll,2)
            control_cost_ind_coll(:,ii) = control_cost_ind_coll(:,ii-1) + ucost' .* input_coll(:,ii);
        end
        % correct for bPos
        for ii = 2888-7 : 7304-7
            control_cost_ind_coll(1,ii) = control_cost_ind_coll(1,ii-1) + 0.01 .* input_coll(1,ii);
        end
        control_cost_ind_coll(1,7304-6:end) = control_cost_ind_coll(1,7304-7) .* ones(1,1295);
    
    %%  
        if (build.setback)
            cd ('SimData'); 
            filename = ['BC_PB_09-predHor' num2str(sim.predhor) '-' build.type '-' build.env '-' build.win '-' build.int '-' 'withSetback' '-' build.weather '-' datestr(now, 'dd_mmmm_yyyy-HH_MM') '.mat'];
            save(filename);
            cd ..
        else
            cd ('SimData'); 
            filename = ['BC_PB_09-predHor' num2str(sim.predhor) '-' build.type '-' build.env '-' build.win '-' build.int '-' 'noSetback' '-' build.weather '-' datestr(now, 'dd_mmmm_yyyy-HH_MM') '.mat'];
            save(filename);
            cd ..
        end

        disp('============ bye! ============')

        if ~serv
           close all
            
            figure();
            hold on; 
            title(['BC\_PB\_09: ' build.weather ', pred\_hor = ' num2str(sim.predhor)  ', total cost: ' num2str(round(control_cost))]);
            plot(ymin_coll(:,1),'r','LineWidth',2); plot(ymax_coll(:,1),'r','LineWidth',2);
            plot(output_coll(:,1),'b','LineWidth',2); 
            legend('y_1 = Troom')
            xlabel('time steps [h]','FontSize', 24)
            ylabel('room temperature [degC]','FontSize',24)
            h_troom = gca;
            ylim([15 30])
            hold off
            
            %%
            
            % plot constraints
            figure();
            hold on; 
            plot(ymin_coll(1:168,1),'r','LineWidth',2);
            plot(ymin_coll(1:168,1)+tightening,'b','LineWidth',2);
            plot(ymax_coll(1:168,1)-tightening,'b','LineWidth',2);
            plot(ymax_coll(1:168,1),'r','LineWidth',2);
            xlabel('time steps [h]','FontSize', 24)
            ylabel('room temperature [degC]','FontSize',24)
            h_troom = gca;
%             ylim([15 30])
            hold off
            
            %%
            tightening = 1*ones(168,1);
            figure();
            hold on; 
            harea = area([ymin_coll(1:168,1), tightening, ...
                ymax_coll(1:168,1)-tightening-(ymin_coll(1:168,1)+tightening),...
                tightening]);
            set(harea,'BaseValue',5)
            set(harea(1),'LineWidth',2)
            set(harea(2),'LineWidth',1)
            set(harea(3),'LineWidth',1)
            set(harea(4),'LineWidth',2)
            set(harea(1),'FaceColor','w')
            set(harea(2),'FaceColor',[0.6,0.6,0.6])
            set(harea(3),'FaceColor','w')
            set(harea(4),'FaceColor',[0.6,0.6,0.6])
            xlabel('time steps [h]','FontSize', 24)
            ylabel('room temperature [degC]','FontSize',24)
            h_troom = gca;
%             ylim([15 30])
            hold off
            
            
            %%
            
            figure();
            hold on; 
            title(['BC\_PB\_09: ' build.weather ', pred\_hor = ' num2str(sim.predhor)  ', total cost: ' num2str(round(control_cost))]);
            plot(state_coll(:,6));
            legend('x6 = Tslab')
            h_x6 = gca;
            hold off

            figure();
            hold on; 
            title(['BC\_PB\_09: ' build.weather  ', pred\_hor = ' num2str(sim.predhor)  ', total cost: ' num2str(round(control_cost))]);
            plot(output_coll(:,2)); plot(ymin_coll(:,2),'--r'); %plot(ymax_all(:,2),'--m');
            legend('y_2 = Eroom')
            hold off

%             figure();
%             hold on; 
%             title(['BC\_PB\_09: ' build.weather  ', pred\_hor = ' num2str(sim.predhor)  ', total cost: ' num2str(round(control_cost))]);
%             plot(output_coll(:,3)); plot(ymin_coll(:,3),'--r'); %plot(ymax_all(:,3),'--m');
%             legend('y_3 = TsrfCeil')
%             hold off
            
            figure(); 
            hold on; 
            title(['BC\_PB\_09: ' build.weather ', pred\_hor = ' num2str(sim.predhor)  ', total cost: ' num2str(round(control_cost))]);
            plot(1:numSimSteps,input_coll(1,:), 1:numSimSteps,input_coll(2,:), 1:numSimSteps,input_coll(3,:), ...
                1:numSimSteps,input_coll(4,:) );
            legend('bPos','eLighting','fcUsgFact','hPowRad');
            h_inputs = gca;
            hold off 
            
            figure(); hold on;
            title(['BC\_PB\_09: ' build.weather ', pred\_hor = ' num2str(sim.predhor)  ', total cost: ' num2str(round(control_cost))]);
            plot(1:numSimSteps,control_cost_ind_coll(1,:), 1:numSimSteps,control_cost_ind_coll(2,:), 1:numSimSteps,control_cost_ind_coll(3,:), ...
                1:numSimSteps,control_cost_ind_coll(4,:) );
            legend('cost bPos','cost eLighting','cost fcUsgFact','cost hPowRad');
            h_costs = gca;
            hold off 
            
%             figure();
%             hold on; 
%             title(['BC\_PB\_09: sim\_hor = ' num2str(sim.predhor) ', pred\_hor = ' num2str(sim.predhor)  ', total cost: ' num2str(round(control_cost))]);
%             plot(control_cost_coll);
%             legend('control cost over hours')
%             hold off
            
            % create predictions / realizations
            % occupancy
            cd('InData/occupancy')
            load('ig_act_real.mat')
            cd ..; cd ..;
            figure(); hold on;
            title('occupancy information')
            plot(ig_act_real,'b','LineWidth',2)
            legend('IG Realization')
            h_ig = gca;
            hold off


            cd(['InData/weatherData/' build.weather]);
            % get new data from 2007 measurements
            load(['Vreal_' build.weather '.mat']);      % load data from weather realization in 2007
            load(['TAreal_' build.weather '.mat']);
            load(['TWreal_' build.weather '.mat']);
            cd ..
            cd ..
            cd ..
            
            figure(); hold on;
            title('Solar radiation')
            plot(V,'b')
            legend('realization')
            h_RGS = gca;
            hold off;
            
            figure(); hold on;
            title('TA')
            plot(TA,'b')
            legend('realization')
            h_TA = gca;
            hold off;

            figure(); hold on;
            title('TW')
            plot(TW,'b')
            legend('realization')
            h_TW = gca;
            hold off;
            
            linkaxes([h_troom; h_inputs; h_x6; h_costs; h_ig; h_RGS; h_TA; h_TW],'x')

            %%
            [v] = get_disturbance(build_stand, sim, 2888-7);  % get disturbance for May-Oct
            Bv_help = [];
            Bx_help = [];
            for j=1:size(bm.Bvu,3)  % sum over all inputs
                Bv_help = [Bv_help,bm.Bvu(:,:,j)*v(5056-2880,:)'];
                Bx_help = [Bx_help,bm.Bxu(:,:,j)*state_coll(5056,:)'];
            end                    
            x_next = bm.A*state_coll(5056,:)' + (bm.Bu + Bv_help + Bx_help)*input_coll(5056,:)' + bm.Bv*v(5056-2880,:)';
            x_next(1)
            
            
            %%
%             keyboard
        end
    end
end

%% functions called inside the main function
function [ucost, ymin,ymin_sb, ymax,ymax_sb, umin, umax] = getCostsAndConstraints(build)          
    
    if strcmp(build.sys,'e01')
        % u = [1bPos 2eLighting 3hPowSlab 4cPowSlab 5fcUsgFact 6hPowRad]
        % for NRPE
        if strcmp(build.seas, 'sum')
            ucost = [ 0.01  3.32  0.976  0.976  7.47  1.107 ];
            % y = [ 1Troom  2Eroom  3TsrfCeil ];
            ymin = [ 22  500  18 ];
            ymin_sb = [ 5  0  5 ];
%             ymin_sb = ymin;
            ymax = [ 26  1.0E+6  1000000 ];   % 26
            ymax_sb = [40  1.0E+6  1000000 ];
%             ymax_sb = ymax;
        else
            ucost = [ 0  3.32  0.976  0.976  7.47  1.107 ];
            ucost = [ 0  3.32  0.976  0.976  7.47  1.107 ];
            % y = [ 1Troom  2Eroom  3TsrfCeil ];
            ymin = [ 21  500  18 ];
            ymin_sb = [ 5  0  5 ];
%             ymin_sb = ymin;
            ymax = [ 25  1.0E+6  1000000 ];   % 25
            ymax_sb = [40  1.0E+6  1000000 ];
%             ymax_sb = ymax;
        end
        
        if (~build.setback)
            ymin_sb = ymin;
            ymax_sb = ymax;
        end
        
        if strcmp(build.type, 'sa')
            if strcmp(build.env, 'h')
                if strcmp(build.win, 'wh')
                    if strcmp(build.int, 'ih')
                        umin = [0      0           0          0          0           0       ];
                        umax = [1      7.2         0          2*40.119     1           2*37.576  ]; 
                    end
                    if strcmp(build.int, 'il')
                        umin = [0      0           0          0          0           0       ];
                        umax = [1      7.2         0          2*28.119     1           2*37.576  ]; 
                    end  
                end
                if strcmp(build.win, 'wl')
                    if strcmp(build.int, 'ih')
                        umin = [0      0           0          0          0           0       ];
                        umax = [1      7.2         0          2*40.653     1           2*22.075  ]; 
                    end
                    if strcmp(build.int, 'il')
                        umin = [0      0           0          0          0           0       ];
                        umax = [1      7.2         0          2*28.653     1           2*22.075  ]; 
                    end
                end
            end
            if strcmp(build.env, 'l')
                if strcmp(build.win, 'wh')
                    if strcmp(build.int, 'ih')
                        umin = [0      0           0          0          0           0       ];
                        umax = [1      7.2         0          2*0.126     1           2*37.375  ]; 
                    end
                    if strcmp(build.int, 'il')
                        umin = [0      0           0          0          0           0       ];
                        umax = [1      7.2         0          2*28.126     1           2*37.379  ]; 
                    end
                end
                if strcmp(build.win, 'wl')
                    if strcmp(build.int, 'ih')
                        umin = [0      0           0          0          0           0       ];
                        umax = [1      7.2         0          2*40.677     1           2*21.384  ]; 
                    end
                    if strcmp(build.int, 'il')
                        umin = [0      0           0          0          0           0       ];
                        umax = [1      7.2         0          2*28.677     1           2*21.384  ]; 
                    end
                end
            end
        end
        if strcmp(build.type, 'pa')
            if strcmp(build.env, 'h')
                if strcmp(build.win, 'wh')
                    if strcmp(build.int, 'ih')
                        umin = [0      0           0          0          0           0       ];
                        umax = [1      7.2         0          2*35.773     1           2*11.488  ]; 
                    end
                    if strcmp(build.int, 'il')
                        umin = [0      0           0          0          0           0       ];
                        umax = [1      7.2         0          2*23.773     1           2*11.488  ]; 
                    end
                end
                if strcmp(build.win, 'wl')
                    if strcmp(build.int, 'ih')
                        umin = [0      0           0          0          0           0       ];
                        umax = [1      7.2         0          2*35.922     1           2*7.163   ]; 
                    end
                    if strcmp(build.int, 'il')
                        umin = [0      0           0          0          0           0       ];
                        umax = [1      7.2         0          2*23.922     1           2*7.163   ]; 
                    end
                end
            end
            if strcmp(build.env, 'l')
                if strcmp(build.win, 'wh')
                    if strcmp(build.int, 'ih')
                        umin = [0      0           0          0          0           0       ];
                        umax = [1      7.2         0          2*35.772     1           2*11.519  ]; 
                    end
                    if strcmp(build.int, 'il')
                        umin = [0      0           0          0          0           0       ];
                        umax = [1      7.2         0          2*23.772     1           2*11.519  ]; 
                    end
                end
                if strcmp(build.win, 'wl')
                    if strcmp(build.int, 'ih')
                        umin = [0      0           0          0          0           0       ];
                        umax = [1      7.2         0          2*35.918     1           2*7.272   ]; 
                    end
                    if strcmp(build.int, 'il')
                        umin = [0      0           0          0          0           0       ];
                        umax = [1      7.2         0          2*23.918     1           2*7.272   ]; 
                    end
                end
            end
        end
    end
    % remove u4
    ucost(4) = [];
    umin(4) = [];
    umax(4) = [];
    
    % remove u3
    ucost(3) = [];
    umin(3) = [];
    umax(3) = [];
end

function [bm] = get_building_model(build)

    cd(['InData/',build.sys,'_occup'])
    multiZapprox = 0;
    filename = sprintf('%s-Or_%s_%s_%s_%s_%s_d_multiZapprox%s.txt', build.sys, build.dir, build.type, build.env, build.win, build.int, num2str(multiZapprox));
    building = filename;

    if strcmp(build.sys, 'e01')   % always have multiZapprox == 0 (implicitly)
        range = 'A11..L22';
        bm.A = dlmread(building, '\t', range);
        range = 'A25..F36';    
        bm.Bu = dlmread(building, '\t', range);
        range = 'A39..H50';
        bm.Bv = dlmread(building, '\t', range);
        range = 'A53..H64';
        bm.Bvu(:,:,1) = dlmread(building, '\t', range);
        range = 'A67..H78';
        bm.Bvu(:,:,2) = dlmread(building, '\t', range);
        range = 'A81..H92';
        bm.Bvu(:,:,3) = dlmread(building, '\t', range);
        range = 'A95..H106';
        bm.Bvu(:,:,4) = dlmread(building, '\t', range);
        range = 'A109..H120';
        bm.Bvu(:,:,5) = dlmread(building, '\t', range);
        range = 'A123..H134';
        bm.Bvu(:,:,6) = dlmread(building, '\t', range);
        range = 'A137..L148';
        bm.Bxu(:,:,1) = dlmread(building, '\t', range);
        range = 'A151..L162';
        bm.Bxu(:,:,2) = dlmread(building, '\t', range);
        range = 'A165..L176';
        bm.Bxu(:,:,3) = dlmread(building, '\t', range);
        range = 'A179..L190';
        bm.Bxu(:,:,4) = dlmread(building, '\t', range);
        range = 'A193..L204';
        bm.Bxu(:,:,5) = dlmread(building, '\t', range);
        range = 'A207..L218';
        bm.Bxu(:,:,6) = dlmread(building, '\t', range);
        range = 'A221..L223';
        bm.C = dlmread(building, '\t', range);
        range = 'A226..F228';
        bm.Du = dlmread(building, '\t', range);
        range = 'A231..H233';
        bm.Dv = dlmread(building, '\t', range);
        range = 'A236..H238';
        bm.Dvu(:,:,1) = dlmread(building, '\t', range);
        range = 'A241..H243';
        bm.Dvu(:,:,2) = dlmread(building, '\t', range);
        range = 'A246..H248';
        bm.Dvu(:,:,3) = dlmread(building, '\t', range);
        range = 'A251..H253';
        bm.Dvu(:,:,4) = dlmread(building, '\t', range);
        range = 'A256..H258';
        bm.Dvu(:,:,5) = dlmread(building, '\t', range);
        range = 'A261..H263';
        bm.Dvu(:,:,6) = dlmread(building, '\t', range);
    end

    cd ..
    cd ..
    
    % remove u4
    bm.Bu(:,4) = [];
    bm.Bxu(:,:,4) = [];
    bm.Bvu(:,:,4) = [];
    bm.Du(:,4) = [];
    bm.Dvu(:,:,4) = [];
    
    % remove u3
    bm.Bu(:,3) = [];
    bm.Bxu(:,:,3) = [];
    bm.Bvu(:,:,3) = [];
    bm.Du(:,3) = [];
    bm.Dvu(:,:,3) = [];
    
end

% ???? questions for  turbance ????
function [v] = get_disturbance(build, sim, RGS_begin)
% RGS_begin sets from where we should start reading in from RGS data

    % Read in weather variables [Tair Twetbulb SolRad]
    cd(['InData/weatherData/' build.weather]);
    % get new data from 2007 measurements
    load(['Vreal_' build.weather '.mat']);      % load data from weather realization in 2007
    load(['TAreal_' build.weather '.mat']);
    load(['TWreal_' build.weather '.mat']);
    cd ..
    cd ..
    cd ..
    
    w(:,1) = TA(RGS_begin:RGS_begin+sim.simhor+sim.predhor-1);
    w(:,2) = TW(RGS_begin:RGS_begin+sim.simhor+sim.predhor-1);
    w(:,3) = V(RGS_begin:RGS_begin+sim.simhor+sim.predhor-1)/3;

    % Read in internal gains variables [Pers  Equip]
    ig = get_internalgains(build, sim, RGS_begin);

    % concatenate disturbances to v
    v = [];
    if strcmp(build.int, 'ih')
        intGPers = 9;
        intGEqui = 15;
    end
    if strcmp(build.int, 'il')
        intGPers = 5;
        intGEqui = 7;
    end
    if strcmp(build.win, 'wh')
        AwinN = NaN;
        AwinE = NaN;
        AwinS = 14.4;
        AwinW = 14.4;
    end
    if strcmp(build.win, 'wl')
        AwinN = NaN;
        AwinE = NaN;
        AwinS = 5.4;
        AwinW = 5.4;
    end
    if strcmp(build.type, 'sa')
        solGFact0 = 0.000364508;
        solGFact1 = 0.01786;        % 56 m2 floor area
        tauVisFact0 = 0;
        tauVisFact1 = 1.025638;
    end
    if strcmp(build.type, 'pa')
        solGFact0 = 0.0002240694;
        solGFact1 = 0.0109794;
        tauVisFact0 = 0;
        tauVisFact1 = 0.9022514;
    end

    % Prepare weather variables  
    % ????? Why i=i+1 in if()?????
    % was wenn sim.wVarNeeded = [1 1 1 1 1 ] -> w overflow!
    totRad = 0;
    i = 1;
    if ( sim.wVarNeeded(1) )    % true
        Tair = w(:,i); 
        i=i+1; 
    end 
    if ( sim.wVarNeeded(2) )    % true
        TfreeCool = w(:,i); 
        i=i+1; 
    end 
    if ( sim.wVarNeeded(3) )    % false
        totRad = totRad + AwinN*w(:,i); 
        i=i+1; 
    end 
    if ( sim.wVarNeeded(4) )    % false
        totRad = totRad + AwinE*w(:,i); 
        i=i+1; 
    end 
    if ( sim.wVarNeeded(5) )    % true
        totRad = totRad + AwinS*w(:,i); 
        i=i+1; 
    end 
    if ( sim.wVarNeeded(6) )    % false
        totRad = totRad + AwinW*w(:,i); 
        i=i+1; 
    end 

    % Prepare internal gains variables
    % ???? kann er nicht vom falschen nehmen, falls sim.igVarNeeded = [0 1] ????
    i = 1;
    if ( sim.igVarNeeded(1) )
        persG  = ig(:,i); 
        i=i+1; 
    end 
    if ( sim.igVarNeeded(2) )
        equipG = ig(:,i); 
        i=i+1; 
    end 

    % Compute disturbances
    v(:,1 ) = solGFact0   * totRad;
    v(:,2 ) = solGFact1   * totRad;
    v(:,3 ) = tauVisFact0 * totRad;
    v(:,4 ) = tauVisFact1 * totRad;
    v(:,5 ) = intGPers    * persG;
    v(:,6 ) = intGEqui    * equipG;
    v(:,7 ) = Tair;
    v(:,8 ) = TfreeCool;
end  

function [ig] = get_internalgains(build, sim, RGS_begin)
    % ig_dist = 2;
    % Load saved occupancy data
    cd('InData/occupancy')
        load('ig_act_real.mat')
        cd ..; cd ..;
        
    ig_Pers = ig_act_real(RGS_begin : RGS_begin+sim.simhor+sim.predhor-1);
    ig_Equip = ig_Pers + 0.1;
    ig = [ig_Pers ig_Equip];
end

function [ymin_all, ymax_all] = get_constraints(sim, v, ymin, ymax, ymin_sb, ymax_sb, RGS_begin)
% cannot use v b/c we don't know v in advance
    cd('InData/occupancy')
    load('ig_act_constr_1we.mat');    % load empirical diurnal profile of internal gain (people)
    cd ..; cd ..;
    ig_Pers = repmat(ig_act_constr_1we,[52,1]); % all 0,1
    ig_Pers = 0.3 * ig_Pers;    % expected 0.4
    ig_Equip = ig_Pers + 0.1;   % same as in old BacLab
    ig = [ig_Pers ig_Equip];
    if (RGS_begin == 2881)  % start from Tuesday
        ig = ig(1+24:end,:);
    elseif RGS_begin == 7297
        ig = ig(1+3*24:end,:);
    end
    ig = ig(1:sim.simhor+sim.predhor,:);

    for i = 1:sim.simhor+sim.predhor
        if ig(i,1) == 0  % v(5) = internal gains (persons)
            ymin_var(i,:) = ymin_sb;
            ymax_var(i,:) = ymax_sb;
        else
            ymin_var(i,:) = ymin;   % ymin is if people are there
            ymax_var(i,:) = ymax;
        end
    end
    ymin_all = ymin_var(1:sim.simhor+sim.predhor,:);
    ymax_all = ymax_var(1:sim.simhor+sim.predhor,:);
end

% no slack variables
function [u u_opt softConstr] = get_control_input(x0, bm, sim, vPred, ymin, ymax, ucost, umin, umax)
    % Set parameters
    nx = size(bm.A,1);
    ny = size(bm.C,1);
    nu = size(bm.Bu,2);
    nv = size(bm.Bv,2);
    N_mpc  = sim.predhor;
    nslack = 2*ny*N_mpc;      % size of slackMax, slackMin

    % reformulate QP in the following, simple shape
    % X = [U slackMax slackMin]
    %       min {ucost_lin' * X} + X'*H*X
    % s.t.  [I ; -I]*U <= [umax ; umin], i.e. A_ineq1*U <= b_ineq1
    %       ymin <= y <= ymax, i.e. A_ineq2*U <= b_ineq2
    %       where, Y = tmp1*x0 + tmp2*V + A_ineq2*U
    
    ucost_lin = [reshape(ucost',N_mpc*nu,1) ; 12^4*ones(nslack,1)];   % linear part
   
    % build inequality 1 out of 5
    A_ineq1 = [eye(N_mpc*nu) ; -eye(N_mpc*nu)];
    A_ineq1 = [A_ineq1 zeros(2*nu*N_mpc,nslack)];
    b_ineq1 = [reshape(umax',N_mpc*nu,1) ; -reshape(umin',N_mpc*nu,1)];
    
    tmp1 = nan(N_mpc*ny,nx);
    for ii = 1 : N_mpc
        tmp1(1+(ii-1)*ny:ii*ny,:) = bm.C*bm.A^(ii-1);
    end

    tmp2 = zeros(N_mpc*ny,nv*N_mpc);
    for ii = 0 : N_mpc-1
        for jj = 0 : ii-1
            tmp2(1+ii*ny:(ii+1)*ny, 1+jj*nv:(jj+1)*nv) = bm.C*bm.A^(ii-1-jj)*bm.Bv;
        end
        tmp2(1+ii*ny: (ii+1)*ny, 1+ii*nv:(ii+1)*nv) = bm.Dv;
    end

    % precompute matrices   B_{u,0},...,B_{u,N-1}
    %                       D_{u,0},...,D_{u,N-1}
    
    Bvu_aux = cell(N_mpc,1);
    Bxu_aux = cell(N_mpc,1);
    Dvu_aux = cell(N_mpc,1);
    
    for jj = 1 : N_mpc  % for inputs u_0,...,u_(N_mpc-1)
        Bvu_aux{jj} = nan(nx,nu);
        Bxu_aux{jj} = nan(nx,nu);
        Dvu_aux{jj} = nan(ny,nu);
            
        for ii = 1 : nu
            Bvu_aux{jj}(:,ii) = bm.Bvu(:,:,ii) * vPred(jj,:)'; % linearized
            Bxu_aux{jj}(:,ii) = bm.Bxu(:,:,ii) * x0;   % linearized
            Dvu_aux{jj}(:,ii) = bm.Dvu(:,:,ii) * vPred(jj,:)';  % linearized
        end 
    end
    
    
    A_ineq2_tmp = zeros(N_mpc*ny,nu*N_mpc);
    for ii = 0 : N_mpc-1
        for jj = 0 : ii-1
            Bu_tmp = bm.Bu + Bvu_aux{jj+1} + Bxu_aux{jj+1};
            A_ineq2_tmp(1+ii*ny:(ii+1)*ny, 1+jj*nu:(jj+1)*nu) = bm.C*bm.A^(ii-1-jj)*Bu_tmp;
        end
        A_ineq2_tmp(1+ii*ny:(ii+1)*ny, 1+ii*nu:(ii+1)*nu) = Dvu_aux{ii+1} + bm.Du;
    end
    
    V = reshape(vPred',N_mpc*nv,1);
    
    A_ineq2 = [A_ineq2_tmp -eye(size(A_ineq2_tmp,1)) zeros(size(A_ineq2_tmp,1),nslack-size(A_ineq2_tmp,1))];
    b_ineq2 = reshape(ymax',N_mpc*ny,1) - tmp1*x0 - tmp2*V;
    
    %build inequality 3 of 4: soft constraint of lower bound
    A_ineq3 = [-A_ineq2_tmp zeros(size(A_ineq2_tmp,1),size(A_ineq2_tmp,1)) -eye(size(A_ineq2_tmp,1)) zeros(size(A_ineq2_tmp,1),nslack-2*size(A_ineq2_tmp,1))];
    b_ineq3 = -reshape(ymin',N_mpc*ny,1) + tmp1*x0 + tmp2*V;
    
    % build inequality 4 of 4: slack variable must be positive
    A_ineq4 = [zeros(nslack,N_mpc*nu) -eye(nslack)];
    b_ineq4 = zeros(nslack,1);
        
    A_ineq = [A_ineq1; A_ineq2 ; A_ineq3 ; A_ineq4];
    b_ineq = [b_ineq1; b_ineq2 ; b_ineq3 ; b_ineq4];
            
    OPTIONS.verbose = 1; % verbose output
    OPTIONS.logfile = 0; % log into file
    OPTIONS.probtype = 0; % Quadratic Problem is 1
    PARAM.int=[1063,2]; % Tipp from Martin
    [ u_opt, fmin, status, details ] = cplexint([],ucost_lin, ...
        A_ineq, b_ineq, [],[],[],[],[],PARAM,OPTIONS);
    
    u = reshape(u_opt(1:nu*N_mpc), nu, N_mpc)';
    
    softConstr = sum(u_opt(end-1:end));

end


function [x_next,output_k] = state_update(build,bm,sim,x,v,u,k)
    x_next = [];
    output_k = [];
    x_k = x;
    
    for i = k:1:k+sim.Tolc-1    % apply in open-loop fashion
        Bv_help = [];
        Bx_help = [];
        for j=1:size(bm.Bvu,3)  % sum over all inputs
            Bv_help = [Bv_help,bm.Bvu(:,:,j)*v(i,:)'];
            Bx_help = [Bx_help,bm.Bxu(:,:,j)*x_k];
        end                    

        x_next = [x_next, bm.A*x_k + (bm.Bu + Bv_help + Bx_help)*u((i-k+1),:)' + bm.Bv*v(i,:)'];
        Dv_help = [];
        for j=1:size(bm.Bvu,3)
            Dv_help = [Dv_help,bm.Dvu(:,:,j)*v(i,:)'];
        end            
        output_k = [output_k, bm.C*x_k + (bm.Du+Dv_help)*u((i-k+1),:)' + bm.Dv*v(i,:)'];    % 1 step behind
        x_k = x_next(:,end);
    end
end 