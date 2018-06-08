function net = cnnff(net, x)
    n = numel(net.layers);
    net.layers{1}.a{1} = x;
    inputmaps = 1;

    for l = 2 : n   %  for each layer
        if strcmp(net.layers{l}.type, 'c')
            %  !!below can probably be handled by insane matrix operations
            for j = 1 : net.layers{l}.outputmaps   %  for each output map
                %  create temp output map
                %jn_infoa(net.layers{l - 1}.a{1},'net.layers{l - 1}.a{1}');
                %fprintf('net.layers{l}.kernelsize:%d\n',net.layers{l}.kernelsize);
                jn_infoa(net.layers{l - 1}.a{1},'net.layers{l - 1}.a{1}');
                z = zeros(size(net.layers{l - 1}.a{1}) - [net.layers{l}.kernelsize - 1 net.layers{l}.kernelsize - 1 0]);
                for i = 1 : inputmaps   %  for each input map
                    %  convolve with corresponding kernel and add to temp output map
                    jn_infoa(net.layers{l - 1}.a{i},'net.layers{l - 1}.a{i}');
                    jn_infoa(net.layers{l}.k{i}{j}, 'net.layers{l}.k{i}{j}');
                    z = z + convn(net.layers{l - 1}.a{i}, net.layers{l}.k{i}{j}, 'valid');
                    jn_infoa(z,'z');
                end
                %  add bias, pass through nonlinearity
                jn_infoa(z,'z');
                net.layers{l}.a{j} = sigm(z + net.layers{l}.b{j});
                jn_infoa(net.layers{l}.a{j}, 'net.layers{l}.a{j}');
            end
            %  set number of input maps to this layers number of outputmaps
            inputmaps = net.layers{l}.outputmaps;
        elseif strcmp(net.layers{l}.type, 's')
            %  downsample
            for j = 1 : inputmaps
               z = convn(net.layers{l - 1}.a{j}, ones(net.layers{l}.scale) / (net.layers{l}.scale ^ 2), 'valid');   %  !! replace with variable
               net.layers{l}.a{j} = z(1 : net.layers{l}.scale : end, 1 : net.layers{l}.scale : end, :);              
               jn_infoa(net.layers{l}.a{j}, 'net.layers{l}.a{j}');
            end
        end
    end

    %  concatenate all end layer feature maps into vector
    net.fv = [];
    for j = 1 : numel(net.layers{n}.a)
        sa = size(net.layers{n}.a{j});
        net.fv = [net.fv; reshape(net.layers{n}.a{j}, sa(1) * sa(2), sa(3))];
    end
    %  feedforward into output perceptrons
    jn_infoa(net.ffW, 'net.ffW');
    jn_infoa(net.fv, 'net.fv');

    rows = size(net.fv, 1);
    cols = size(net.fv, 2);
    fprintf('%d %d\n', rows, cols);
    for i=1:rows
        for j=1:cols
            fprintf('%f ', net.fv(i,j));
        end
        fprintf('\n');
    end

    rows = size(net.ffW, 1);
    cols = size(net.ffW, 2);
    fprintf('%d %d\n', rows, cols);
    for i=1:rows
        for j=1:cols
            fprintf('%f ', net.ffW(i,j));
        end
        fprintf('\n');
    end

    jn_infoa(net.ffb, 'net.ffb');
    temp = net.ffW * net.fv;
    jn_infoa(temp, 'multiply');
    net.o = sigm(net.ffW * net.fv + repmat(net.ffb, 1, size(net.fv, 2)));
    jn_infoa(net.o, 'net.o');
end
