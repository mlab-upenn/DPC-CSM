function C = sdp_userfun(u, umin, umax, epsvars, horizon)
    % Constraints only on user variables
    C = epsvars >= 0;
    add_control_bounds('blind', 1);
    add_control_bounds('light', 2);
    add_control_bounds('heating', 3);
    add_control_bounds('cooling', 4);
    
    function add_control_bounds(fn, k)
        if isfield(u, fn)
            C = [C, umin(k) <= u.(fn) <= umax(k)];
        end
    end
end