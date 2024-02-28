% Divide nodes within the WM function network
tic
% load atlas of single functional network
img_nii= load_nii('label2.nii');
data=img_nii.img;
judge_label=1;
left_num=0;                        
right_num=0;
data_size=size(data);
% Statistics about the number of white matter
for i=1:data_size(1)
    for j=1:data_size(2)
        for k=1:data_size(3)
           if (data(i,j,k)==judge_label && i<=round(data_size(1)/2))
               left_num=left_num+1;
           end
           if (data(i,j,k)==judge_label && i>round(data_size(1)/2))
               right_num=right_num+1;
           end
        end
    end
end

num_wm=length(find(data==judge_label));               
position_wm_left=zeros(left_num,3);
position_wm_right=zeros(right_num,3);
position_wm=zeros(num_wm,3);                
temp_left=1;
temp_right=1;
temp_wm=1;
% The voxel positions corresponding to the left and right white matter were counted
for i=1:data_size(1)
    for j=1:data_size(2)
        for k=1:data_size(3)
           if (data(i,j,k)==judge_label && i<=round(data_size(1)/2))    % left
               position_wm_left(temp_left,:,:)=[i,j,k];
               temp_left=temp_left+1;
           end
           if (data(i,j,k)==judge_label && i>round(data_size(1)/2))     % right
               position_wm_right(temp_right,:,:)=[i,j,k];
               temp_right=temp_right+1;
           end
           if (data(i,j,k)==judge_label)                                % whole brain
               position_wm(temp_wm,:,:)=[i,j,k];
               temp_wm=temp_wm+1;
           end
        end
    end
end

disp(['***Number of voxels on the left',num2str(left_num)]);
disp(['***Number of voxels on the right',num2str(right_num)]);
grow_num=0;
%Calculate how many regions you need to grow based on the voxel proportion of the corresponding brain region
seed_left=5;   
seed_right=5;
%% 
seed_all=seed_left+seed_right;
b1=ones(seed_all,1);
minb1=0;
t=1;
grow_num_max = 0;
minb1_best = 0; 
while (grow_num<1140 || minb1<90 )  % set the minimun number to generate 
                %Generate random seed spots in the left and right brain                 
The global location index corresponding to the right
First 64-left and then 64-right random number integration


                left_rand=randperm(left_num,seed_left)';
                right_rand=randperm(right_num,seed_right)';
                right_rand=right_rand+left_num;              %The global location index corresponding to the right
                rand_all=[left_rand;right_rand];             %First 64-left and then 64-right, random number integration

                use_judge=zeros(num_wm,1);                   %Whether it was a seed point
                value_judge=zeros(num_wm,1);                 %Whether there is value
                cell_all=num2cell(rand_all,seed_all);        %Count where all seed points grow
                

                % Regional growth method
                data_grow=data;                                               
                old_num=1;
                num_grow=ones(seed_all,1);
                while(grow_num<num_wm && old_num~=grow_num)
                    old_num=grow_num;                        %Determines whether the num is the same as the previous NUM. If the NUM is the same, the growth stops
                    grow_temp=cell(128,1);                   %This time go through all the seed points grown to new seed points
                    num_temp=zeros(num_wm,1);                %Greater than 1-> There are repetitions being grown to the location
                    for i=1:seed_all
                        temp_length=length(cell_all{i});
                        for j=1:temp_length

                            if use_judge(cell_all{i}(j))==0

                                   position_seed=position_wm(cell_all{i}(j),:,:);
                                   value_judge(cell_all{i}(j))=1;            % Set the seed point to 1
                                    for l=-1:1
                                        for m=-1:1
                                            for n=-1:1
                                               seed_x=position_seed(1)+l;
                                               seed_y=position_seed(2)+m;
                                               seed_z=position_seed(3)+n;
                                               if (seed_x>0 && seed_y>0 && seed_z>0) && data(seed_x,seed_y,seed_z)==judge_label

                                                       position_row=get_position(position_wm,[seed_x,seed_y,seed_z]);
                                                       if value_judge(position_row)==0

    %                                                        data_grow(seed_x,seed_y,seed_z)=i;
                                                           if (l~=0 || m~=0 || n~=0)  && length(cell_all{i})<160
                                                                 if ~ismember(position_row,grow_temp{i})   %position_row不在grow_temp{i}中，加入grow_temp{i}
                                                                 grow_temp{i}=[grow_temp{i} position_row];
                                                                 num_temp(position_row)=num_temp(position_row)+1;
                                                                 end
    %                                                            cell_all{i}=[cell_all{i} position_row]; 
    %                                                            num_grow(i)=num_grow(i)+1;
                                                           else
                                                               data_grow(seed_x,seed_y,seed_z)=i;
                                                           end

    %                                                        value_judge(position_row)=1;
                                                       end


                                               end


                                            end
                                        end
                                    end
                                    use_judge(cell_all{i}(j))=1;

                            end


                        end

                    end


                               %Statistical processing of all seed points found
                               same_grow=find(num_temp>1);                 %Repeat the found cell
                               if ~isempty(same_grow)
                                   for iter1=1:length(same_grow)           %Repeatedly found values are processed one by one
                                       same_num=same_grow(iter1);
                                       same_list=[];
                                       seednum_list=[];

                                       for iter2=1:seed_all                 %Count the number of duplicate rows and their seeds
                                          if ismember(same_num,grow_temp{iter2})
                                              same_list=[same_list iter2];
                                              seednum_list=[seednum_list num_grow(iter2)];
                                          end 
                                       end

                                      minseed=min(seednum_list);            %mininum value     
                                      minnum=find(seednum_list==minseed);   %The row corresponding to the minimum value
                                      save_num=0;
                                      if(length(minnum)>1)                  %If the row corresponding to the minimum value is greater than 2, select one at random
                                          save_num=same_list(minnum(randperm(numel(minnum),1)));
                                      else
                                          save_num=same_list(minnum);
                                      end

                                      for iter3=1:length(same_list)         %Eliminate the same element
                                          if same_list(iter3)~=save_num
                                              temp_row=same_list(iter3);
                                              grow_temp{temp_row}=grow_temp{temp_row}(~ismember(grow_temp{temp_row},same_num));
                                          end
                                      end


                                   end
                               end
                               %Merge to cell_all
                               for iter4=1:seed_all
                                  cell_all{iter4}=[cell_all{iter4} grow_temp{iter4}];
                                  for iter5=1:length(grow_temp{iter4})
                                      value_judge(grow_temp{iter4}(iter5))=1;
                                      temp_x=position_wm(grow_temp{iter4}(iter5),1);
                                      temp_y=position_wm(grow_temp{iter4}(iter5),2);
                                      temp_z=position_wm(grow_temp{iter4}(iter5),3);
                                      data_grow(temp_x,temp_y,temp_z) =iter4;
                                  end
                               end
                    grow_num=sum(value_judge);
                end
                % move single region
                data_erase=data; 

                for i=1:length(value_judge)
                   if value_judge(i)==0
                      position_temp=position_wm(i,:,:);
                      data_erase(position_temp(1),position_temp(2),position_temp(3))=0;
                      data_grow(position_temp(1),position_temp(2),position_temp(3))=0;
                   end
                end


                nii_save=img_nii;
                nii_erase=img_nii;
                nii_save.img=data_grow;
                nii_erase.img=data_erase;

    num_test=0;
    for i=1:seed_all
        num_test=num_test+length(cell_all{i});
    end

     
    for j1=1:seed_all
      b1(j1,1)=length(cell_all{j1});
    end
    minb1=min(b1);
    
    if grow_num_max<grow_num
        grow_num_max = grow_num;
    end
    
    if minb1>minb1_best
        minb1_best = minb1;
        grow_num_minb1 = grow_num;
    end
    disp(['-------grow2---------']);
    disp(['*********',num2str(t),'********']);
    t=t+1;
    disp(['grow_num:',num2str(grow_num),'  minb1:',num2str(minb1)]);
    disp(['best————grow_num:',num2str(grow_num_minb1),'  minb1:',num2str(minb1_best),' Maximum number of calories:',num2str(grow_num_max)]);

end

% end
toc
function [position]= get_position(position_all,position_judge)
    for i=1:length(position_all)
        position_temp=position_all(i,:,:);
        if isequal(position_temp,position_judge)
            position=i;
        end
    end
end
