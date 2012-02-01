% function participant=readParticipantsList
%
% Reads the participantDetails.xls and creates a structure
%  It then finds associated data from the allParticipantsByInitials folder
%  and adds it to the structure.
%
% The individual files are then added to an allParticipants folder to form
%  the database for further analysis

restorePath=path;

disp('running...')
addpath('..\profiles\allParticipantsByInitials')

[num,txt,txt]=xlsread('participantDetails');
[r c]=size(txt);
for i=1:r
    participant(i).number=num(i,1);
    participant(i).IH_NH=txt{i,2};
    participant(i).code=txt{i,3};
    participant(i).initials=txt{i,4};
    participant(i).ear=txt{i,5};
    participant(i).gender=txt{i,6};
    participant(i).DOB=txt{i,7};
    participant(i).birthYear=num(i,8);
    participant(i).startTest=num(i,9);
    participant(i).age=num(i,10);
    participant(i).fullProfile=num(i,11);
    participant(i).longProfile=num(i,12);
    participant(i).shortProfile=num(i,13);
    participant(i).tinnitus=txt{i,14};
    participant(i).conductiveLoss=txt{i,17};

    % read profiles and add to structure
    cmd=[ 'x=profile_' participant(i).initials '_L;'];
    eval(cmd)
    code=participant(i).code;
    participant(i).leftEarData=x;

    % add file to pile
    copyFrom=['allParticipantsByInitials\profile_'...
        participant(i).initials '_L.m'];
    copyTo=['allParticipants\profile_' code '_L.m'];
    copyfile(copyFrom, copyTo)

    cmd=[ 'x=profile_' participant(i).initials '_R;'];
    eval(cmd)
    participant(i).rightEarData=x;

    % add file to pile
    copyFrom=['allParticipantsByInitials\profile_' participant(i).initials '_R.m'];
    copyTo=['allParticipants\profile_' code '_R.m'];
    copyfile(copyFrom, copyTo)

end

save participantCompendium participant

disp(['done: saved as ''participantCompendium.mat'''])
path(restorePath)
