function thStack = removeBackground(stack, sig, neuSize)

if nargin < 3
    neuSize = [10, 10 , 1];
end

if nargin <2 || isempty(sig)
    sig = 1;
end

if gpuDeviceCount > 1
    useGPU = true;
else
    useGPU = false;
end

stack = mat2gray(stack);

% set nans to median
nan_list = isnan(stack(:));
stack(nan_list) = median(stack(:), 'omitnan');

[nY, nX, nImg] = size(stack);

if useGPU
    stack = gpuArray(stack);
end

fStack = imgaussfilt(stack, sig);

% if nImg > neuSize % I tried and didn't notice any advantage of 3D strel
% se = strel3D([neuSize neuSize neuSize*1.5]);
% thStack = imtophat(fStack, se);
% else
fStack = imgaussfilt(stack, sig);
se = strel('disk', neuSize);
thStack = imtophat(fStack, se);
% end

if useGPU
    thStack = gather(thStack);
end


end
