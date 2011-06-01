function createSCP(experimentName, filePath)

d = dir([filePath, '*.map']);
fid = fopen([experimentName '.scp'],'w');

for I = 1:size(d,1)
    fprintf(fid,[d(I).name(1:end-3) 'map\n']);
end;

fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORIGINAL STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d = dir('*.08');
% fid = fopen('largeTripletTrain3.scp','w');
% 
% for I = 1:size(d,1)
%     fprintf(fid,[d(I).name(1:end-2) 'map\n']);
% end;
% 
% fclose(fid);
