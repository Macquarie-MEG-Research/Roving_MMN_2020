for j=1:length(folders)

    SubjectID = folders(j).name;
    cd([orig,'/',SubjectID,'/ReTHM/'])
    
    coreg_output = [pwd '\\MEMES_old\\']; % where to store the output from MEMES
    
        % move the MEMES output into the coreg_output folder.
        if ~exist(coreg_output, 'dir')
            mkdir(coreg_output);
        end                
        % the files to move:
        % grad_trans, headshape, headmodel, trans_matrix, sourcemodel3d, shape
        movefile('*trans*', coreg_output);
        movefile('*shape*', coreg_output);
        movefile('*model*', coreg_output);
        movefile('*quality*', coreg_output);
        %movefile('*example*', coreg_output);
        %movefile('*scaling*', coreg_output);
        movefile('*realigned*', coreg_output);
        %movefile('*winner*', coreg_output);
        movefile('*MEMES*', coreg_output);
        movefile('*error_age*', coreg_output);
end