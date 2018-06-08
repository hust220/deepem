function r = randmat( m, n )
    lcg = LCG;
    r = zeros(m, n);
    for i=1:m
        for j=1:n
            r(i,j) = lcg.rand(((i-1)*n+j)*100);
        end
    end
end

