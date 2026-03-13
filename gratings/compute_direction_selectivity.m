function dsi = compute_direction_selectivity(meanResp, directions, prefIdx)
% COMPUTE_DIRECTION_SELECTIVITY Compute DSI = (Rp - Ra) / (Rp + Ra).
%
% Rp is the mean response at the preferred direction.
% Ra is the mean response at the direction 180 degrees away.
% If the exact opposite is not sampled, the nearest sampled direction is used.

    dsi = NaN;

    if isempty(prefIdx) || ~isfinite(prefIdx)
        return
    end

    Rp = meanResp(prefIdx);
    oppIdx = find_opposite_direction_idx(directions, prefIdx);
    Ra = meanResp(oppIdx);

    if ~isfinite(Rp) || ~isfinite(Ra)
        return
    end

    denom = Rp + Ra;
    if denom == 0
        return
    end

    dsi = (Rp - Ra) / denom;
end
