function signals = gp_create_signals(X, Y, normalization)
% GP_CREATE_SIGNALS - Create SignalValues object from raw data
%   X, Y are loaded by gp_load_data.m

    if nargin < 3, normalization = true; end
    normgiven = isstruct(normalization);

    % Note that in the raw data, inputs at step k correspond to output at step k+1.

    % For Time Of Day, convert sawtooth signal to periodic signal by cosine
    X.tod = cos(X.tod/24*2*pi);
    % Same for Day of Week
    X.dow = cos(X.dow/7*2*pi);

    % Create the signal model to manage all the signals
    signals = SignalsValues();
    for f = fieldnames(X)'
        fn = f{1};
        addSignal(fn, X.(fn));
    end

    addSignal('temp', Y(2:end));

    function addSignal(fn, vals)
        do_norm = false;
        if normgiven
            do_norm = isfield(normalization, fn);
            if do_norm
                vmin = normalization.(fn).min;
                vmax = normalization.(fn).max;
            end
        else
            if normalization && ~ismember(fn, {'tod', 'dow'})   % , 'temp', 'blind', 'light', 'heating', 'cooling'
                % Do not normalize certain inputs
                vmin = min(vals);
                vmax = max(vals);
                % Normalization if the values are far enough from 0 or if (vmax-vmin) large or small enough
                do_norm = min(abs(vmin), abs(vmax)) > 5 || (vmax-vmin) > 5 || (vmax-vmin) < 0.01;
            end
        end
        if do_norm
            signals.addSignal(fn, vals, [vmin, vmax]);
        else
            signals.addSignal(fn, vals);
        end
    end
end