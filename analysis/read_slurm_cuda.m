function exetime = read_slurm_cuda(root,nline)
files = dir(fullfile(root,'slurm*.out'));
exetime = zeros(numel(files),1);
for i = 1:numel(files)
    fileID = fopen(fullfile(root,files(i).name),'r');
    for j = 1:nline
        tline = fgetl(fileID);
    end
    pos1 = strfind(tline,'time')+5;
    pos2 = strfind(tline,'s')-2;
    exetime(i) = str2double(tline(pos1:pos2(end)));
    fclose(fileID);
end

end