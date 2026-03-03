function movie = clipMovie(movie, dx, dy)

[h, w, ~] = size(movie);

dxMax = max(0, ceil(max(dx)));
dxMin = min(0, floor(min(dx)));
dyMax = max(0, ceil(max(dy)));
dyMin = min(0, floor(min(dy)));
validX = (1 + dxMax):(w + dxMin);
validY = (1 + dyMax):(h + dyMin);

movie = movie(validY,validX,:);

end