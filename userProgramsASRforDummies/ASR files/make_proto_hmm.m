% function make_proto_hmm(filename,featureType,numFeatures,numStates)
% filename is obvious
% featureType is usually USER_D_A
% numFeatures is number of features (including differences if used)
% numStates is usually 18

function make_proto_hmm(filename,featureType,numFeatures,numStates)

ofp = fopen(filename,'w');

fprintf(ofp,'<BeginHMM>\n');
fprintf(ofp,' <NumStates> %d <VecSize> %d <%s> <nullD> <diagC>\n',numStates,numFeatures,featureType');
fprintf(ofp,' <StreamInfo> 1 %d\n',numFeatures);
for state=2:numStates-1,
    fprintf(ofp,' <State> %d <NumMixes> 1\n',state);
    fprintf(ofp,'  <Stream> 1\n');
    fprintf(ofp,'  <Mixture> 1 1.0\n');
    fprintf(ofp,'    <Mean> %d\n',numFeatures);
    fprintf(ofp,'      ');
    fprintf(ofp,'%1.1f ',zeros(1,numFeatures));
    fprintf(ofp,'\n');
    fprintf(ofp,'   <Variance> %d\n',numFeatures);
    fprintf(ofp,'      ');
    fprintf(ofp,'%1.1f ',ones(1,numFeatures));
    fprintf(ofp,'\n');
end
fprintf(ofp,'<TransP> %d\n',numStates);
transp = zeros(numFeatures);
transp(1,2)=1;
for state=2:numStates-2,
    transp(state,state)=0.6;
    transp(state,state+1)=0.4;
end
transp(numStates-1,numStates-1)=0.9;
transp(numStates-1,numStates)=0.1;
for state=1:numStates,
    fprintf(ofp,'%1.3e  ',transp(state,1:numStates));
    fprintf(ofp,'\n');
end
fprintf(ofp,'<EndHMM>\n');

size(transp)

fclose(ofp);