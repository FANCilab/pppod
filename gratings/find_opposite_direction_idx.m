function oppIdx = find_opposite_direction_idx(directions, prefIdx)
% FIND_OPPOSITE_DIRECTION_IDX Find sampled direction closest to preferred+180 deg.

    prefDir = directions(prefIdx);
    targetOpp = mod(prefDir + 180, 360);

    circularDiff = abs(mod(directions - targetOpp + 180, 360) - 180);
    [~, oppIdx] = min(circularDiff);
end
