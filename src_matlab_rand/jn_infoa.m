function a = jn_infoa( array, name )
%JN_INFOA 此处显示有关此函数的摘要
%   此处显示详细说明
fid = fopen('aa.txt', 'a+');
sz = size(array);
fprintf(fid, '%s %s min:%f max:%f sum:%f mean:%f\n', name, mat2str(sz), min(array(:)), max(array(:)), sum(array(:)), mean(array(:))); 
fclose(fid);

end

