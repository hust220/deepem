function result = cnntest_m(net, x)
    %  feedforward
    net = cnnff(net, x);
    %[~, h] = max(net.o);
    %[~, a] = max(y);
    %bad = find(h ~= a);
    result=net.o;
    %er = numel(bad) / size(y, 2);
end