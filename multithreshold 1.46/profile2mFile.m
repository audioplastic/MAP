function profile2mFile(longTone, shortTone, gaps, BFs, TMC, offBFs, ...
    IFMCs, mFileName, location)
% profile2mFile creates an .m file suitable for use with plotProfile.
% mFileName is saved in 'location' inside the 'profiles' folder
%
% The arguments are based on matrices that can be copied and pasted 
%   from the Excel file for the participant.
%   Alternatively, they are supplied directly by multiThreshold after 
%   a 'profile' paradigm execution

restorePath=path;
addpath(['..' filesep 'utilities'])

St = ['function x = ' mFileName ];
filePath=['..'  filesep 'profiles' filesep location filesep mFileName '.m'];
fid = fopen(filePath,'w');
fprintf(fid, '%s\n', St);

St = ['%created: ' UTIL_timeStamp];
fprintf(fid, '%s\n\n', St);

St = ['x.BFs = [' num2str(BFs) '];' ];
fprintf(fid, '%s\n\n', St);

St = ['x.LongTone = [' num2str(longTone',3) '];' ];
fprintf(fid, '%s\n', St);

St = ['x.ShortTone = [' num2str(shortTone',3) '];' ];
fprintf(fid, '%s\n\n', St);

St = ['x.Gaps = [' num2str(gaps',3) '];' ];
fprintf(fid, '%s\n', St);

TMCFreq = BFs;
St = ['x.TMCFreq = [' num2str(TMCFreq) '];' ];
fprintf(fid, '%s\n', St);

fprintf(fid, '%s\n', 'x.TMC = [');
for i = 1:size(TMC,1),
    for j = 1:size(TMC,2),
        St = [num2str(TMC(i,j),3)];
        fprintf(fid, '%s\t', St);
    end
    fprintf(fid,'%s\n',' ');
end
fprintf(fid,'%s\n','];');

St = ['x.TMC = x.TMC'';' ];
fprintf(fid, '%s\n\n', St);

St = ['x.MaskerRatio = [' num2str(offBFs',2) '];' ];
fprintf(fid, '%s\n', St);

IFMCFreq = BFs;
St = ['x.IFMCFreq = [' num2str(IFMCFreq) '];' ];
fprintf(fid, '%s\n', St);

fprintf(fid, '%s\n', 'x.IFMCs = [');
for i = 1:size(IFMCs,1),
    for j = 1:size(IFMCs,2),
        St = [num2str(IFMCs(i,j),3)];
        fprintf(fid, '%s\t', St);
    end
    fprintf(fid,'%s\n',' ');
end
fprintf(fid,'%s\n','];');

St = ['x.IFMCs = x.IFMCs'';' ];
fprintf(fid, '%s\n', St);

fclose(fid);

path(restorePath)
