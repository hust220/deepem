function r = shuffle(n, seed)
    % https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
    lcg = LCG(seed);
    m = 1:n;
    for i=n:-1:2
        j = floor(lcg.rand()*i)+1;
        t = m(i);
        m(i) = m(j);
        m(j) = t;
    end
    r = m;
end

