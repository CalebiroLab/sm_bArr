function [] = preprocessing_from_utrack_format(movie_list,list_groups,global_folders,list_channel,start,inter)
group_names=fieldnames(movie_list);
movie_name={};
for n_group=list_groups
    group=group_names{n_group};
    movie_name=[movie_name,strsplit(movie_list.(group),' ')];
end
movie_name=movie_name(start:inter:end);

for n_mov=1:numel(movie_name)
    for channel=list_channel
        frame_duration=30;%ms
        filename=[movie_name{n_mov},'-C',num2str(channel)];
        file_name_msk=[global_folders.rawfolder, filesep, filename, '_msk.tif'];
        file_name_detection=[global_folders.rawfolder, filesep, filename, '_detection.mat'];
        if isfile([global_folders.rawfolder, filesep, filename, '_tracking.mat'])
            disp([filename,' computing'] )
            %% Scripts to run after original utrack detection and tracking
            copyfile ([global_folders.rawfolder, filesep, filename, '_tracking.mat'],[global_folders.localfolder, filesep])
            utrack_to_OBJ_v1(filename, frame_duration, [global_folders.localfolder, filesep]);
            %% Reshape the data to remove merges and splits (saved in OBJ2)
            OBJ_reshape(filename, global_folders);
            %% Copy files to rawfolder
            copyfile ([global_folders.localfolder, filesep, filename, '_gui2.mat'], [global_folders.rawfolder, filesep]);
            copyfile ([global_folders.localfolder, filesep, filename, '_gui2_steps.mat'], [global_folders.rawfolder, filesep]);
            %% Delete temp files
            delete ([global_folders.localfolder, filesep, filename, '*.*']);
        else
            disp(['did not find tracking ',filename]);
        end
    end
end
end

