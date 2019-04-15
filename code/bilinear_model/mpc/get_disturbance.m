function [v, vPred_all] = get_disturbance(build, sim, RGS_begin)
% RGS_begin sets from where we should start reading in from RGS data

    % Read in Realizations of weather variables [Tair Twetbulb SolRad]
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
    
    % realizations of internal gains variables [Pers  Equip]
    ig = get_internalgains(build, sim, RGS_begin,'real');

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
    
    % process predictions, improved by AR4 model
    cd(['InData/weatherData/' build.weather]);
    load(['dataTAP_' build.KF num2str(sim.predhor) '_' build.weather '.mat']);      % load data from weather realization in 2007
    load(['dataTWP_' build.KF num2str(sim.predhor) '_' build.weather '.mat']);      % load data from weather realization in 2007
    load(['dataRGSP_' build.KF num2str(sim.predhor) '_' build.weather '.mat']);      % load data from weather realization in 2007
    cd ..
    cd ..
    cd ..

    dataTAP = dataTAP_KF_AR4;
    dataTWP = dataTWP_KF_AR4;
    dataRGSP = dataRGSP_KF_AR4;
    
%     predHor = 61;
%     dataV = nan(length(V)-predHor+1,predHor);
%     dataTA = nan(length(V)-predHor+1,predHor  );
%     dataTW = nan(length(V)-predHor+1,predHor);
%     for ii = 1 : length(V)-predHor+1
%         dataV(ii,:) = V(ii:ii+predHor-1);
%         dataTA(ii,:) = TA(ii:ii+predHor-1);
%         dataTW(ii,:) = TW(ii:ii+predHor-1);
%     end
    
%     dataTAP = dataTA;
%     dataTWP = dataTW;
%     dataRGSP = dataV;
    
    wP = cell(sim.simhor,1);
    for ii = 1 : sim.simhor
        wP{ii}(:,1) = dataTAP(RGS_begin+ii-1,1:sim.predhor);
        wP{ii}(:,2) = dataTWP(RGS_begin+ii-1,1:sim.predhor);
        wP{ii}(:,3) = dataRGSP(RGS_begin+ii-1,1:sim.predhor)/3;     % data from RG->RGS
    end
    
    % realizations of internal gains variables [Pers  Equip]
    igP = get_internalgains(build, sim, RGS_begin,'pred');
    
    i = 1;
    if ( sim.wVarNeeded(1) )    % true
        TairP = cell(sim.simhor,1);
        for ii = 1 : sim.simhor
            TairP{ii} = wP{ii}(:,i); 
        end
        i=i+1; 
    end 
    if ( sim.wVarNeeded(2) )    % true
        TfreeCoolP = cell(sim.simhor,1);
        for ii = 1 : sim.simhor
            TfreeCoolP{ii} = wP{ii}(:,i);
        end
        i=i+1; 
    end
    
    totRadP = cell(sim.simhor,1);
    for ii = 1 : sim.simhor
        totRadP{ii} = 0;
    end
    if ( sim.wVarNeeded(3) )    % false
        for ii = 1 : sim.simhor
            totRadP{ii} = totRad{ii} + AwinN*wP{ii}(:,i); 
        end
        i=i+1; 
    end 
    if ( sim.wVarNeeded(4) )    % false
        for ii = 1 : sim.simhor
            totRadP{ii} = totRadP{ii} + AwinE*wP{ii}(:,i); 
        end
        i=i+1; 
    end 
    if ( sim.wVarNeeded(5) )    % true
        for ii = 1 : sim.simhor
            totRadP{ii} = totRadP{ii} + AwinS*wP{ii}(:,i); 
        end
        i=i+1; 
    end 
    if ( sim.wVarNeeded(6) )    % false
        for ii = 1 : sim.simhor
            totRadP{ii} = totRadP{ii} + AwinW*wP{ii}(:,i); 
        end
        i=i+1; 
    end 
    
    % Prepare internal gains variables
    % ???? kann er nicht vom falschen nehmen, falls sim.igVarNeeded = [0 1] ????
    i = 1;
    if ( sim.igVarNeeded(1) )
        persGP  = igP(:,i); 
        i=i+1; 
    end
    if ( sim.igVarNeeded(2) )
        equipGP = igP(:,i); 
        i=i+1; 
    end 
    
    % fetch predicted internal gain
    
    vPred_all = cell(sim.simhor,1);
    for ii = 1 : sim.simhor
        vPred_all{ii}(:,1) = solGFact0 * totRadP{ii};
        vPred_all{ii}(:,2) = solGFact1 * totRadP{ii};
        vPred_all{ii}(:,3) = tauVisFact0 * totRadP{ii};
        vPred_all{ii}(:,4) = tauVisFact1 * totRadP{ii};
        vPred_all{ii}(:,5) = intGPers * persGP(ii:ii+sim.predhor-1);
        vPred_all{ii}(:,6) = intGEqui * equipGP(ii:ii+sim.predhor-1);
        vPred_all{ii}(:,7) = TairP{ii};
        vPred_all{ii}(:,8) = TfreeCoolP{ii};
    end
end  