function m = metrics_eda(t, x, xe, theta)
% METRICS_EDA  Summary metrics of one simulated response.
%   t     : time vector of the run (h), its own grid
%   x     : states, one row per time (column 9 = eda, uM)
%   xe    : equilibrium of the scenario (before the dose)
%   theta : therapeutic threshold for M4 (uM)

    eda  = x(:,9);              % extracellular dopamine (uM)
    deda = eda - xe(9);         % increase over its own equilibrium (uM)

    % M1 - Peak: the highest increase of eda over its baseline.
    % How much the dose raises extracellular dopamine at most.
    [m.peak,index] = max(deda);

    % M2 - Time of the peak: when that maximum happens (h after the dose).
    m.t_peak = t(index);

    % M3 - AUC: area under the increase curve over the simulated window (uM*h).
    % Total extra dopamine exposure; depends on the window (here 0-24 h).
    m.auc = trapz(t, deda);

    % M4 - Time above the threshold: hours in which eda (absolute) stays
    % above theta. Reed's "therapeutic window" with the rule-based theta.
    % baseline_above = true means the baseline itself is already above theta.
    m.t_above = time_above(t, eda, theta);
    m.baseline_above = xe(9) > theta;

    % M5 - Half-peak duration: hours in which the increase stays above half of
    % its own peak. Width of the response, no external threshold needed.
    m.t_half = time_above(t, deda, 0.5*m.peak);

    % M6 - Relative rise: peak divided by the baseline (dimensionless).
    % How many times eda rises over its own baseline (comparable to rat data).
    m.rel_rise = m.peak/xe(9);
end

function T = time_above(t, y, level)
% Time (h) in which y > level, with linear interpolation inside each step.
    T = 0;
    for k = 1:numel(t)-1
        y1 = y(k);  y2 = y(k+1);  dt = t(k+1) - t(k);
        if y1 >= level && y2 >= level
            T = T + dt;                         % whole step above
        elseif y1 < level && y2 < level
            % whole step below: nothing to add
        else
            tc = (level - y1)/(y2 - y1)*dt;     % time from t(k) to the crossing
            if y1 >= level
                T = T + tc;                     % starts above, falls below
            else
                T = T + (dt - tc);              % starts below, rises above
            end
        end
    end
end