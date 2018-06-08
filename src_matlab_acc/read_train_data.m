function a=read_train_data(im,b,start_image_num,end_image_num,mrc_directory,box_directory,particle_image_size,mic_width,mic_length)
    %im indicate the prefix  of the image name
    %b indicate the bites of the number in the image name
    num=0;
    for image_num=start_image_num:end_image_num
        c=im;
        c3=num2str(image_num,['%0',num2str(b,'%01d'),'d']);
        file1=[mrc_directory,c,c3,'.mrc'];
        file2=[box_directory,c,c3,'.box'];
        boxsize=particle_image_size;

        if ~exist(file1) continue; end
        if ~exist(file2) continue; end

        mig=mrcs_read(file1,'b');
        %jn_infoa(mig, 'mig');
        mig=imrotate(flipud(mapstd(mig)),180,'nearest'); 
        %jn_infoa(mig, 'mig');

        box=load(file2);
        %jn_infoa(box, 'box');
        if(size(box,1)==0) continue; end

        for i=1:size(box,1) %In case of some boxes out of bounds
            num=num+1;
            box(i,1)=round(box(i,1));
            box(i,2)=round(box(i,2));

            if(box(i,1)<1) box(i,1)=1; end
            if(box(i,2)<1) box(i,2)=1; end

            if(box(i,1)>(mic_width-boxsize)) box(i,1)=(mic_width-boxsize); end   
            if(box(i,2)>(mic_length-boxsize)) box(i,2)=(mic_length-boxsize); end
            k = mig(box(i,1):box(i,1)+boxsize-1,box(i,2):box(i,2)+boxsize-1);
            a(1:boxsize,1:boxsize,num)=k;
        end
    end
    jn_infoa(a, 'a');
end
