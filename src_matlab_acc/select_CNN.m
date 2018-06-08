
cal=zeros(rotation_n,boxsize*boxsize);
for image_num=start_mic_num:end_mic_num
    tic;
    disp(image_num);
    c3=num2str(image_num,['%0',num2str(name_length,'%2d'),'d']);
    file=[data_path,'mic/',name_prefix,c3,'.mrc'];
    
    if ~exist(file)
        continue;
    end
    
    mig=mrcs_read(file,'b');
    mig=imrotate(flipud(mig),180,'nearest');
    jn_infoa(mig,'mig');
    rows=reshape(mig,size(mig,1)*size(mig,2),1);
    overall_var=var(rows);
    
    num=0;
    cal_r=0;
    for i=1:scan_step:dim_x-boxsize+1
        for j=1:scan_step:dim_y-boxsize+1
            
            par_img=mig(i:i+boxsize-1,j:j+boxsize-1);
            jn_infoa(par_img,'par_img');
            
            for k=1:rotation_n
                im_rot=imrotate(par_img,rotation_angel*(k-1),'nearest');
                cal(k,:)=im_rot(:);
            end
            
            par_std=std(cal,0,2);
            
            if(par_std(1)<max_std&&par_std(1)>min_std)
                cal=mapstd(cal);
                cal2= double(reshape(cal',boxsize,boxsize,4));
                jn_infoa(cal2, 'cal2');
                
                resu=cnntest_m(cnn,cal2);
                jn_infoa(resu,'resu');
                
                num=num+1;
                cal_r(num,1)=i;
                cal_r(num,2)=j;
                cal_r(num,3)=sum(resu)/4;
                cal_r(num,4)=par_std(1);
                cal_r(num,5:8)=resu(1:4);
                
            end
            
        end
    end
        
    for ll=1:9
        threhold=0.1*ll;%selection threhold
        box=0; %output the box file
        nnewbox=0;%final output the selected box file and saved in .box file
        
        j=0;
        for i=1:size(cal_r,1);
            if (cal_r(i,3)>threhold)
                j=j+1;
                box(j,1:3)=cal_r(i,1:3);
            end
        end
        
        if(size(box,1)<2)
            continue;
        end
        
        n=0;
        newbox=0;
        box_n=0;
        for i=1:size(box,1);
            m=0;
            box_m=0;
            for j=1:size(box,1);
                if ((abs(box(j,1)-box(i,1))<=range1)&&(abs(box(j,2)-box(i,2))<=range1))
                    m=m+1;
                    box_m(m,1:3)=box(j,1:3);
                end
            end
            
            n=n+1;
            box_n(n,1:3)=box_m(1,1:3);
            
            for k=2:m
                if(box_m(k,3)>box_n(n,3))
                    box_n(n,1:3)=box_m(k,1:3);
                end
                
                if(box_m(k,3)==box_n(n,3))
                    if((box_m(k,1)>box_n(n,1))||(box_m(k,2)>box_n(n,2)))
                        box_n(n,1:3)=box_m(k,1:3);
                    end
                end
            end
        end
        
        newbox=unique(box_n,'rows');
        
        n=0;
        box_n=0;
        for i=1:size(newbox,1);
            m=0;
            box_m=0;
            
            for j=1:size(newbox,1);
                if ((abs(newbox(j,1)-newbox(i,1))<=range2)&&(abs(newbox(j,2)-newbox(i,2))<=range2))
                    m=m+1;
                    box_m(m,1:3)=newbox(j,1:3);
                end
            end
            
            n=n+1;
            box_n(n,1:3)=box_m(1,1:3);
            for k=2:m
                if(box_m(k,3)>box_n(n,3))
                    box_n(n,1:3)=box_m(k,1:3);
                end
                
                if(box_m(k,3)==box_n(n,3))
                    if((box_m(k,1)>box_n(n,1))||(box_m(k,2)>box_n(n,2)))
                        box_n(n,1:3)=box_m(k,1:3);
                    end
                end
            end
        end
        
        nnewbox=unique(box_n,'rows');
        
        c4=num2str(ll,'%02d');
        
        result_path=[data_path,c4,'result'];
        
        if ~exist(result_path)
            mkdir(result_path);
        end
        
        box_name=[result_path,'/',name_prefix,c3,'.box'];
        fid=fopen(box_name,'w');
        disp(box_name);
        for i=1:size(nnewbox,1)
            fprintf('%d %d %d %d\n',nnewbox(i,1),  nnewbox(i,2),  boxsize, boxsize);
            fprintf(fid,'%d %d %d %d\n',nnewbox(i,1),  nnewbox(i,2),  boxsize, boxsize);
        end
        fclose(fid);
        
    end
    toc;
end
