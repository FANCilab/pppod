function err = learnSources(scaling, data, mixGuess, waveL)

mixScale = mixGuess;
mixScale(:,2) = mixGuess(:,2).*scaling;
tot = sum(mixGuess, 2);
mixScale(:,1)  = tot - mixScale(:,2).*scaling;

[~, ~, err] = prism.unmixSources(data, mixScale, waveL);


end