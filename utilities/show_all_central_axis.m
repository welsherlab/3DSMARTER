% sync of 2 structs and do some analysis
% this is a template for re-analysis of hot spots;

clear;cab

load('Z:\Manuscripts\metal particle tracking\1 usec data\200304 TR003 filopodia\res time data\res time data final sorted.mat')

data1 = sortedS;

load('Z:\Manuscripts\metal particle tracking\1 usec data\200304 TR003 filopodia\figures\12 segments params\12 seg data final sorted.mat')

data2 = sortedS;

clearvars sortedS
figure;hold on
cmap = jet(13);

for i = 1:length(data2)
   data2(i).cyl_shown = false;
    
end


for i = 1:length(data1)
   t1 = data1(i).t_begin;
   t2 = data1(i).t_end;
   
   match_found = false;
   for i2 = 1:length(data2)
       t1p = data2(i2).data_begin_time;
       t2p = data2(i2).data_end_time;
       if t1p == t1 && t2p == t2
           match_found = true;
           %data2(i2).cyl_shown = true;
           break,
       end
   end
   
   %%%% insert code here
   if ~strcmp(data2(i2).ref_fname,'default')
       load(data2(i2).ref_fname)
       
   end
   
   
   
   f_sz = 8;

main_dir = 'Z:\Manuscripts\metal particle tracking\1 usec data\200304 TR003 filopodia';
fig_dir = 'Z:\Manuscripts\metal particle tracking\1 usec data\200304 TR003 filopodia\figures';
data_dir = 'Z:\Manuscripts\metal particle tracking\1 usec data\200304 TR003 filopodia\res time data';

load(fullfile(main_dir ,data2(i2).cyl_fname))

t0 = data2(i2).data_begin_time;

clrs = Paired(12);

LW = 1;

Q1 = [mean(xx(1,:)),mean(yy(1,:)),mean(zz(1,:))];
Q2 = [mean(xx(2,:)),mean(yy(2,:)),mean(zz(2,:))];

if ~data2(i2).cyl_shown
    pcyl{i2} = plot3([Q1(1),Q2(1)],[Q1(2),Q2(2)],[Q1(3),Q2(3)]);
    pcyl{i2}.LineWidth = 2;
    pcyl{i2}.Color = cmap(i2,:);
    
    data2(i2).cyl_shown = true;
    lgd_str{i2} = [num2str(data2(i2).data_begin_time) '~' num2str(data2(i2).data_end_time) 's'];
end

   
   %%%%
    
   if~match_found
      error('no matched line found in second sturct!') 
   end
    
    
end

view(3);axis equal;
legend(lgd_str)

