function r = randmat( m, n )
%   https://en.wikipedia.org/wiki/Linear_congruential_generator
r = zeros(m, n);
for i=1:m
    for j=1:n
        r(i,j) = randnum((i-1)*n+j-1);
    end
end

end

