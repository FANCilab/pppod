function [prefVal, prefIdx] = compute_preferred_value(meanResp, paramValues)
% COMPUTE_PREFERRED_VALUE Return parameter value with maximum mean response.

    prefVal = NaN;
    prefIdx = NaN;

    if isempty(meanResp) || all(~isfinite(meanResp))
        return
    end

    meanResp = meanResp(:);
    valid = isfinite(meanResp);
    validIdx = find(valid);

    [~, localIdx] = max(meanResp(valid));
    prefIdx = validIdx(localIdx);
    prefVal = paramValues(prefIdx);
end
