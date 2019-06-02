function f = sdp_objfun(ymean, yvar, u, epsvars, yref, includedcontrols, Q, R, elecCost, epsc, horizon)
    f = 0;
    % Only include relevant terms in R
    R = R(includedcontrols, includedcontrols);
    for k = 1:horizon
        f = f + ...
            control_cost('blind', 1, k) + ...
            control_cost('light', 2, k) + ...
            control_cost('heating', 3, k) + ...
            control_cost('cooling', 4, k);
        
        f = f + Q*(ymean(k) - yref)^2 + control_quadcost(k) + epsc * epsvars(k); % + 1*yvar{k};
    end
    
    function c = control_cost(name, idx, k)
        c = 0;
        if includedcontrols(idx)
            % If the given control is used by the GP
            c = u.(name)(k) * elecCost(k, idx);
        end
    end
    
    function c = control_quadcost(k)
        uk = [];
        if includedcontrols(1)
            uk = [uk; u.blind(k)];
        end
        if includedcontrols(2)
            uk = [uk; u.light(k)];
        end
        if includedcontrols(3)
            uk = [uk; u.heating(k)];
        end
        if includedcontrols(4)
            uk = [uk; u.cooling(k)];
        end
        c = uk' * R * uk;
    end
end