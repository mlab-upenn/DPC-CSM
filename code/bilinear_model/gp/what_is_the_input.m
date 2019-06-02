function inputs = what_is_the_input(input_model, idx, vals)
    %WHAT_IS_THE_INPUT Find the inputs corresponding to indices in idx.
    %   Used to find irrelevant inputs with ARD kernels.
    %   If vals is provided, it must be of the same length as idx and will
    %   be added to the results (e.g., to include the lengthscale of each
    %   feature).
    
    if nargin > 2 && ~isempty(vals)
        assert(length(vals) == length(idx));
    else
        vals = [];
    end
    
    % Flatten the input model
    flat_inputs = flatten(input_model, {});
    
    inputs = flat_inputs(idx, :);
    if ~isempty(vals)
        for k = 1:length(idx)
            if iscell(vals)
                inputs{k, 3} = vals{k};
            else
                inputs{k, 3} = vals(k);
            end
        end
    end
end

function outs = flatten(ins, prev)
    outs = prev;
    if iscell(ins)
        l = length(ins);
        if l == 2 && ~iscell(ins{1}) && isnumeric(ins{2}) && isvector(ins{2})
            % Atomic case: {signal, autoregressive-indices}
            for k = ins{2}(:)'
                outs(end+1,1:2) = {ins{1}, k}; %#ok<AGROW>
            end
        else
            for k = 1:l
                outs = flatten(ins{k}, outs);
            end
        end
    else
        outs(end+1,1:2) = {ins, 0};
    end
end
