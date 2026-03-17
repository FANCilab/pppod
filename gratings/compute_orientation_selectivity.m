function osi = compute_orientation_selectivity(meanResp, orientations)
% COMPUTE_ORIENTATION_SELECTIVITY Compute OSI as 1 - circular variance.
%
% Orientation is 180-degree periodic, so the vector sum uses exp(2i*theta).
% This implementation returns:
%   OSI = abs(sum(R(theta) * exp(2i*theta)) / sum(R(theta)))
% which equals 1 minus the circular variance.

    meanResp = meanResp(:);
    orientations = orientations(:);
    valid = isfinite(meanResp) & isfinite(orientations);

    if ~any(valid)
        osi = NaN;
        return
    end

    R = meanResp(valid);
    theta = deg2rad(orientations(valid));
    denom = sum(R);

    if denom == 0
        osi = NaN;
        return
    end

    osi = abs(sum(R .* exp(1i * 2 * theta)) / denom);
end
