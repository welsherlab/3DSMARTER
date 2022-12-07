% add path of diffusion data into 13 params data
% no need to save data in server 

clear;cab

% import data containing time stamps and other related information of
% individual manually labelled AOIs 

load(fullfile(pwd,'/data','AOI data.mat'))

%% visualization of raw data without filtering and post-processing
temp = [];
temp2 = [];

for i = 1:length(data1)
   curr_data = data1(i).res_data_comp;
   for j = 1:length(curr_data)
      temp = [temp,curr_data{j}.total_bin_val]; 
      temp2 = [temp2, curr_data{j}.res_dens];
   end
    
    
end


figure; scatter(temp./1e6,temp2,  10)
xlabel('res time (s), total')
ylabel('bin fill rate')

%% maximum residence duration of each event?
% one single residence event is defined as the particle did not leave the
% area (of 95% confidence interval) for a continuous 100-ms interval

temp = [];

for i = 1:length(data1)
   curr_data = data1(i).res_data_comp;
   curr_max_res_t = 0;
   for j = 1:length(curr_data)
       if curr_data{j}.total_bin_val>curr_max_res_t
         curr_max_res_t  = curr_data{j}.total_bin_val;
       end
      
     
   end
    
    temp = [temp, curr_max_res_t]; 
end

figure;histogram(temp./1e6)
xlabel('main peak res time, s')

% 19 spots have <1s main peak, and should not be considered

%% flag residence events that have time stamps grouped at beginning and ending of each 50~60s segment


temp = [];

for i = 1:length(data1)
   curr_data = data1(i).res_data_comp;
    sub_temp = [];
   for j = 1:length(curr_data)
sub_temp = [sub_temp,curr_data{j}.total_bin_val];
      
     
   end
    temp_idx = find(sub_temp == max(sub_temp));
    temp = temp_idx(1);
    temp_data = curr_data{temp_idx}; % get main interval
     
    curr_rg = data1(i).t_end-data1(i).t_begin;
    
    if temp_data.start_bin_idx==1||abs(temp_data.end_bin_idx-curr_rg*10) < 3
       % potential bad data 
       
       if temp_data.start_bin_idx==1
          lb = 1; % residence starts at beginning of trajectory segment
       else
           lb = 2; % residence not finished at end of segment
       end
       
       disp(['job. no. ' num2str(data1(i).job_no) ' may not be complete' ', main interval duration is ' num2str(temp_data.total_bin_val./1e6) ', label is ' num2str(lb) ])
    end
    
end

%% label duration of each entry of leaving inside each AOI

temp = [];

plt_x = {}; 
plt_y = {};

% remove idx with main pk < 1s
temp = [];

for i = 1:length(data1)
   curr_data = data1(i).res_data_comp;
   curr_max_res_t = 0;
   for j = 1:length(curr_data)
       if curr_data{j}.total_bin_val>curr_max_res_t
         curr_max_res_t  = curr_data{j}.total_bin_val;
         
       %else
           %data1(i).is_bad_data = 'main pk < 1s';
       end
      
     
   end
    
    temp = [temp, curr_max_res_t]; 
end

idx1 = find(temp<1e6);

for i = 1:length(data1)
    if ismember(i,idx1)
       data1(i).is_bad_data = 'main pk < 1s'; 
    else
        data1(i).is_bad_data = 'data ok.';
    end
    
end

% remove idx with edged main pk
idx2 = [];
for i = 1:length(data1)
   curr_data = data1(i).res_data_comp;
    sub_temp = [];
   for j = 1:length(curr_data)
sub_temp = [sub_temp,curr_data{j}.total_bin_val];
      
     
   end
    temp_idx = find(sub_temp == max(sub_temp));
    temp = temp_idx(1);
    temp_data = curr_data{temp_idx}; % get main interval
     
    curr_rg = data1(i).t_end-data1(i).t_begin;
    
    if temp_data.start_bin_idx==1||abs(temp_data.end_bin_idx-curr_rg*10) < 3
       % potential bad data 
       idx2 = [idx2, i];
       
       if temp_data.start_bin_idx==1
          lb = 1; 
          data1(i).is_bad_data = 'edged main pk at beginning';
       else
           lb = 2;
           data1(i).is_bad_data = 'edged main pk at ending';
       end
       
      
    end
    
end

remove_idx = unique([idx1,idx2]);

for i = 1:length(data1)
   if ~ismember(i,remove_idx)
       curr_data = data1(i).res_data_comp;
       
       temp = [];
       idx_needed = [];
       for j = 1:length(curr_data)
          if curr_data{j}.total_bin_val>1e6 && curr_data{j}.res_dens>.5
              temp = [temp,curr_data{j}.total_bin_val];
              idx_needed = [idx_needed, j];
          end
           
       end
       
       rounded_res_time = temp(end:-1:1);
       rounded_res_time = round(rounded_res_time./1e6./.1).*.1;
       data1(i).res_time_individual_visits = rounded_res_time;
       
       
       if length(idx_needed)==1
          time_between_intervals = nan;
          data1(i).n_visits = 1;
       else
          time_between_intervals = [];
          for i1 = 1:length(idx_needed)-1
              int1 = curr_data{idx_needed(i1)};
              int2 = curr_data{idx_needed(i1+1)};
             time_between_intervals(i1)= int2.bin_pos(1) - int1.bin_pos(end) ;
              
              
          end
          data1(i).time_between_intervals = time_between_intervals;
          data1(i).n_visits = i1+1;
       end
       
       temp= sort(temp);
       plt_x = [plt_x,1:length(temp)];
       plt_y = [plt_y,temp(end:-1:1)];
   end
    
    
end

% show filtered data
figure; hold on
clrs = jet(length(plt_x));
for i = 1:length(plt_x)
   pp = plot(plt_x{i},plt_y{i});
   pp.Color = clrs(i,:);
    
end


temp = [];
for i = 1:length(plt_y)
    temp = [temp,plt_y{i}];
end



%% how many visits per hot spot

n_visits_ct = [0,0,0,0,0,0];
for i = 1:length(data1)
   if strcmp(data1(i).is_bad_data,'data ok.')
       curr_n = data1(i).n_visits;
       n_visits_ct(curr_n) = n_visits_ct(curr_n)+1;
       
   end
    
    
end

%% add time stamps into structured data

load('Z:\Manuscripts\metal particle tracking\1 usec data\200304 TR003 filopodia\figures\12 segments params\12 seg data final sorted with diff data v2.mat')

for i = 1:length(data1)
   if strcmp(data1(i).is_bad_data,'data ok.')
      match_found = false;
      for j = 1:length(data2)
          if data1(i).t_begin == data2(j).data_begin_time && data1(i).t_end == data2(j).data_end_time
             match_found = true;
              break,
          end
          
          
      end
      if ~match_found
         error('something is wrong!') 
      end
       
   end
    
   disp([' i = ' num2str(i) ' complete.'])
    
end


