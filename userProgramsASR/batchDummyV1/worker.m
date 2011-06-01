function worker(workFolder)

%unit testing
% close all; clear all; clc
% workFolder = 'D:\exps\_bar';


%main script
load(fullfile(workFolder,'jobObject.mat'))
x=obj;
clear obj;
x.initMAP; %Need to alert it to the path

personalWork = 0;
while(any(x.todoStatus==0))        
    x.lockJobList;
    x = x.loadSelf; %Reload incase changed
    todoNow = find(~x.todoStatus,1,'first'); %Grab 1st open job
    x.todoStatus(todoNow) = 1; %Flag it as pending
    x.storeSelf; %store pending flag as quickly as possible to minimise race condition impact
    x.unlockJobList;
        
    % ---  DO WORK  ---
    x.genFeat(x.wavList(todoNow).name);
    % --- END OF WORK ---
    
    x.lockJobList;
    x = x.loadSelf; %Reload incase changed while processing (probably has)
    x.todoStatus(todoNow) = 2; %Flag as complete
    x.storeSelf; %Update as done immediately    
    x.unlockJobList;
        
    clc
    personalWork = personalWork+1;
    disp( ['This process has completed ' num2str(personalWork) ' jobs'] )
    x.checkStatus        
end

disp('-*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*-')
disp( ' > COMPLETED CURRENT JOB' )
disp( ['  In the folder ' workFolder '  .....'] )
disp( ['  This process completed ' num2str(personalWork) ' jobs'] )
disp('-*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*-')


