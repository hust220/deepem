function a = jn_infoa( array, name )
%JN_INFOA �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
fid = fopen('aa.txt', 'a+');
sz = size(array);
fprintf(fid, '%s %s min:%f max:%f sum:%f mean:%f\n', name, mat2str(sz), min(array(:)), max(array(:)), sum(array(:)), mean(array(:))); 
fclose(fid);

end

