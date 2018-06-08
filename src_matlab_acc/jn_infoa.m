function jn_infoa( array, name )

sz = size(array);
fprintf('%s %s min:%f max:%f sum:%f mean:%f\n', name, mat2str(sz), min(array(:)), max(array(:)), sum(array(:)), mean(array(:))); 

end

