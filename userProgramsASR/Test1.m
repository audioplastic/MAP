close all; clear all;  clear classes

xL = jobject('L', fullfile('testExp2','featR'));

xL.MAPopHSR = 1;
xL.MAPopMSR = 0;
xL.MAPopLSR = 0;

xL.participant = 'MPa';
xL.MAPuseEfferent = 0;
xL.numWavs = 10;
xL = xL.assignFiles;
xL.storeSelf;



% NOWdone = 0;
% tic
% while(~all(x.todoStatus==2))
%     clc              
%     x.lockJobList; %lock still needed even with a load incase it blocks a save
%     x = x.loadSelf;    
%     x.unlockJobList;
%     x.checkStatus
%     pause(10)
% end
% disp('ALL DONE!!')
% toc
