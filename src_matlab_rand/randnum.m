function r = randnum(s)
% 此处显示有关此函数的摘要
% https://en.wikipedia.org/wiki/Linear_congruential_generator
persistent modules;
persistent a;
persistent c;
persistent seed;
if (isempty(modules)) 
    modules = 2^32-1; 
end
if (isempty(a))
    a = 1103515245;
end
if (isempty(c))
    c = 12345;
end
if (isempty(seed))
    seed = 1; 
end
seed = mod(a * s + c, modules);
r = seed/modules;
end

