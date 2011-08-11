function method=MAPparamsDEMO(BFlist, sampleRate)
% MAPmodelParameters establishes a complete set of MAP parameters
%  The parameters are global
%
% Parameter file names must be of the form <MAPparams> <name>
%
% input arguments
%   BFlist (optional) specifies the desired list of channel BFs
%     if none is supplied or if BFlist= -1, the values are set below.
%   sampleRate (optional), default is 20000. this value is set in the
%     global variable stimulusparameters and should be picked up by the user.

global inputStimulusParams outerMiddleEarParams DRNLParams IHC_ciliaParams
global IHC_RPParams IHCpreSynapseParams  AN_IHCsynapseParams
global MacGregorParams MacGregorMultiParams  filteredSACFParams
global experiment betweenRuns % used by multithreshold only

showParams=0;  % if set it writes all parameters to the command window

% set model sample rate using input arguments
if nargin>1
    method.dt=1/sampleRate; % i.e. as sample rate
else
    % the sample rate is set here if sampleRate is not supplied as an argument
    method.dt=1/20000;
end

if nargin<1 || (length(BFlist)==1 && BFlist<0)  % missing BFlist or BFlist==-1
    % set model BFs here if either
    %  there is no argument specifying BFlist
    %  BF is a single valued vector set to -1 (as used by multithreshold).
    BFlist=[250 500 1000 2000 4000 8000];
end
method.nonlinCF=BFlist;
method.PSTHbinWidth=1e-3;

% module #1 inputStimulus
%  no parameters required but the structure must exist
inputStimulusParams.dt=-1;
inputStimulusParams.boundaryValue=[];

% module #2 outerMiddleEar
outerMiddleEarParams.stapesScalar=	6e-008;
outerMiddleEarParams.dt=	-1;
outerMiddleEarParams.fft=	0;
% highpass filter [gain, cutoff, order]
outerMiddleEarParams.filters=	[0  900     2];
outerMiddleEarParams.boundaryValue=[];

% module #3 DRNL
DRNLParams.linOrder=	3;           % order of gammatone filters in linear path
DRNLParams.NumLinLPfilters=	0;       % not used: number of low pass filters (linear path)
DRNLParams.nonlinOrder=	3;           % order of gammatone filters in nonlinear path
DRNLParams.NumNonlinLPfilters=	0;   % not used: number of low pass filters (nonlinear path)
DRNLParams.ID=	'nonlinBW     a          b           linCF       linBW     lin GAIN    c';
DRNLParams.dt=	-1;                  % can be used to change sample rate (-1 means don't change)
DRNLParams.ATTdB=0*ones(length(BFlist),1);  % attenutation of nonlinear path to simulate efferent suppression (dB)
DRNLParams.boundaryValue=[];
%additional DRNL parameters established using interpolation are shown below

% module #4 IHC_cilia
IHC_ciliaParams.tc=	0.0003;                %  filter time constant to convert BM velocity to cilia displacement
IHC_ciliaParams.C=	0.1;                      % scalar (C_cilia )
IHC_ciliaParams.Gmax=	8e-009;      % maximum conductance (Siemens)
IHC_ciliaParams.u0=	7e-009;              % Boltzman function parameters to determine conductance re input
IHC_ciliaParams.s0=	5e-007;
IHC_ciliaParams.u1=	7e-009;
IHC_ciliaParams.s1=	5e-009;

IHC_ciliaParams.Gu0=	1.974e-009;  % resting conductance

% module #5 IHC_RP
IHC_RPParams.Cab=	16e-012;         % IHC capacitance (F)
IHC_RPParams.Et=	0.1;            % endocochlear potential (V)
IHC_RPParams.Gk=	1.8e-008;        % potassium conductance (S)
IHC_RPParams.Ek=	-0.0705;         % potassium equilibrium potential
IHC_RPParams.Rpc=	0.04;            % combined resistances
IHC_RPParams.boundaryValue=[];

% module #6 IHCpreSynapse
IHCpreSynapseParams.GmaxCa=	8e-009;   % maximum calcium conductance
IHCpreSynapseParams.ECa=	0.066;          % calcium equilibrium potential
IHCpreSynapseParams.beta=	400;            % parameters determining channel opening
IHCpreSynapseParams.gamma=	130;        % parameters determining channel opening
IHCpreSynapseParams.tauM=	0.0001;      % membrane time constant

IHCpreSynapseParams.tauCa=	0.0003;      % calcium clearance (removal)  time constant
IHCpreSynapseParams.power=	3;              %
IHCpreSynapseParams.z=	2e42;                  % scalar Ca -> vesicle release rate
IHCpreSynapseParams.boundaryValue=[];

% module #7 AN_IHCsynapse
AN_IHCsynapseParams.mode=	'spikes';
% AN_IHCsynapseParams.mode=	'probability';
AN_IHCsynapseParams.M=	12;                 % maximum number of vesicles in available store
AN_IHCsynapseParams.y=	3;                     % rate of replacement of lost vesicles (vesicles/s)
AN_IHCsynapseParams.l=	2580;               % rate of loss of vesicles from the cleft
AN_IHCsynapseParams.x=	30;                  % rate of replenishment from the re-uptake store
AN_IHCsynapseParams.r=	6580;                % rate of reuptake from the cleft into the cell
AN_IHCsynapseParams.refractory_period=	0.00075;
AN_IHCsynapseParams.numFibers=	50; % when spikes are generated this is the number of fibers at each BF
AN_IHCsynapseParams.PSTHbinWidth=	0.001;
AN_IHCsynapseParams.TWdelay=0.004; % estimate of delay between stimulus and fiber effects (not used often)
AN_IHCsynapseParams.plotSynapseContents=0; % special plot of synapse contents
AN_IHCsynapseParams.returnSynapseContents=0; % as third argument instead of PSTH. Cell array - {c, q,w, kt}
AN_IHCsynapseParams.boundaryValue=[];

MacGregorMultiType='chopper';  %'primary-like'
switch MacGregorMultiType
    case 'chopper'
        % module #8 MacGregorMulti
        MacGregorMultiParams.type = 'chopper cell';
        MacGregorMultiParams.dendriteLPfreq=	50;         % lowpass filter cutoff to simulate dendritic effects
        MacGregorMultiParams.currentPerSpike=	2e-008;  % current (A) for each spike delivered by input fibers
        MacGregorMultiParams.tauM=	0.002;                        % membrane time constant (s)
        MacGregorMultiParams.tauGk=	0.001;                        % potassium conductance time constant
        MacGregorMultiParams.tauTh=	0.02;                          % time constant of adaptation of the threshold
        MacGregorMultiParams.c=	0;                                        % threshold shift rate
        MacGregorMultiParams.b=	8000;                                 % potassium conductance shift after a spike
        MacGregorMultiParams.Th0=	0.01;                            % threshold (V)
        MacGregorMultiParams.R=	60000000;                         % cell resistance (ohms)
        MacGregorMultiParams.Ek=	-0.01;                            % potassium equilibrium potential
        MacGregorMultiParams.Eb=	0.06;                             % (cosmetic) voltage added during an actionpotential
        MacGregorMultiParams.Er=	-0.06;                             % resting potential (V)
        MacGregorMultiParams.fibersPerNeuron=	20;        % number of randomy selected input fibers
        MacGregorMultiParams.PSTHbinWidth=	0.0001;
        MacGregorMultiParams.nNeuronsPerBF=	10;           % number of macGregor neurons per BF
        MacGregorMultiParams.boundaryValue=[];
        
    case 'primary-like'
        % module #8 MacGregorMulti
        % primary-like cell
        MacGregorMultiParams.type = 'primary-like cell';
        MacGregorMultiParams.dendriteLPfreq=	200;         % lowpass filter cutoff to simulate dendritic effects
        MacGregorMultiParams.currentPerSpike=	.15e-6;% current (A) for each spike delivered by input fibers
        MacGregorMultiParams.tauM=	5e-4;                            % membrane time constant (s)
        MacGregorMultiParams.tauGk=	0.0012;                        % potassium conductance time constant
        MacGregorMultiParams.tauTh=	0.015;                          % time constant of adaptation of the threshold
        MacGregorMultiParams.c=	0.01;                                     % threshold shift rate
        MacGregorMultiParams.b=	500;                                     % potassium conductance shift after a spike
        MacGregorMultiParams.Th0=	0.01;                              % threshold (V)
        MacGregorMultiParams.R=	220000000;                         % cell resistance (ohms)
        MacGregorMultiParams.Ek=	-0.01;                             % potassium equilibrium potential
        MacGregorMultiParams.Eb=	0.06;                              % (cosmetic) voltage added during an actionpotential
        MacGregorMultiParams.Er=	-0.06;                              % resting potential (V)
        MacGregorMultiParams.fibersPerNeuron=	4;           % number of randomy selected input fibers
        MacGregorMultiParams.PSTHbinWidth=	1e-4;
        MacGregorMultiParams.nNeuronsPerBF=	10;             % number of macGregor neurons per BF
        MacGregorMultiParams.boundaryValue=[];
end

% module #9 MacGregor
MacGregorParams.comment=	'Meddis 2006 Table V new currentPerSpike';
MacGregorParams.dendriteLPfreq=	100;        % lowpass filter cutoff to simulate dendritic effects

MacGregorParams.currentPerSpike=3e-007; % current (A) for each spike delivered by input fibers
MacGregorParams.tauM=	0.0005;                  % membrane time constant (s)
MacGregorParams.tauGk=	0.001;                    % potassium conductance time constant
MacGregorParams.tauTh=	0.02;                      % time constant of adaptation of the threshold
MacGregorParams.c=	0;                                    % threshold shift rate
MacGregorParams.b=	8000;                             % potassium conductance shift after a spike
MacGregorParams.Th0=	0.01;                        % threshold (V)
MacGregorParams.R=	60000000;                     % cell resistance (ohms)
MacGregorParams.Ek=	-0.01;                            % potassium equilibrium potential
MacGregorParams.Eb=	0.06;                              % (cosmetic) voltage added during an actionpotential
MacGregorParams.Er=	-0.06;                             % resting potential (V)
MacGregorParams.PSTHbinWidth=	0.001;
MacGregorParams.boundaryValue=[];

% module #10 filteredSACF
minPitch=	300; maxPitch=	3000; numPitches=60;
pitches=100*log10(logspace(minPitch/100, maxPitch/100, numPitches));
filteredSACFParams.lags=1./pitches;         % autocorrelation lags vector
filteredSACFParams.acfTau=	.003;              % time constant of running autocorrelation
filteredSACFParams.lambda=	0.12;            % slower filter to smooth ACF
filteredSACFParams.plotFilteredSACF=1;   % 0 plots unfiltered ACFs
filteredSACFParams.plotACFs=0;                 % special plot (see code)
%  filteredSACFParams.usePressnitzer=0; % attenuates ACF at  long lags
filteredSACFParams.lagsProcedure=  'useAllLags';  %    'omitShortLags',  'useBernsteinLagWeights'
% filteredSACFParams.lagsProcedure=  'useBernsteinLagWeights';  %    'omitShortLags',  'useBernsteinLagWeights'
filteredSACFParams.criterionForOmittingLags=3;
filteredSACFParams.boundaryValue=[];


% DRNL individual channel parameters
DRNLa=10e5;                 % nonlinear path gain (below compression threshold)
DRNLa=5e5;                 % nonlinear path gain (below compression threshold)
% DRNLa=0;                 % nonlinear path gain (below compression threshold)

DRNLb=8e-5 ;           % sets compression threshold (m/s)
DRNLb=2e-5 ;           % sets compression threshold (m/s)

DRNLg=500;                  % linear path gain factor
% DRNLg=5000;                  % linear path gain factor
% DRNLg=0;                  % linear path gain factor

c=0.2;
fixedBFs=[250 500 1000 2000 4000 6000 8000 12000];
p=.000246;              %linear function (p,q) computes nonlinear gammatone bandwidths
q=1.410546;
DRNLParams.p=p;
DRNLParams.q=q;
nlBWs=fixedBFs./(p*fixedBFs+q);

%  ' nonlinBW	  a      b      linCF   linBW	linGAIN	 c  nonlinCF')
% parameters at any BF can be computed by interpolation here
interpolationTable=...
    [
    nlBWs(1)   	DRNLa	DRNLb	  250	100	    DRNLg	c	fixedBFs(1);
    nlBWs(2)    DRNLa	DRNLb	  500	400 	DRNLg	c	fixedBFs(2);
    nlBWs(3)    DRNLa	DRNLb	1000	500 	DRNLg	c	fixedBFs(3);
    nlBWs(4)    DRNLa	DRNLb	2000	700     DRNLg	c	fixedBFs(4);
    nlBWs(5)    DRNLa	DRNLb	4000	1700	DRNLg	c	fixedBFs(5);
    nlBWs(6)    DRNLa	DRNLb	6000	2200	DRNLg	c	fixedBFs(6);
    nlBWs(7)    DRNLa	DRNLb	8000	2750	DRNLg	c	fixedBFs(7);
    nlBWs(8) 	DRNLa	DRNLb	12000	3000	DRNLg	c	fixedBFs(8)
    ];

[r c] =size (interpolationTable);
if ~isequal(r, length(fixedBFs) )
    % to catch operator error if interpolationTable is changed
    error('MAPparams: the requested BF is outside of the range of allowed values')
end

% 	Create new filterbank by interpolating among manual parameters
channelParameters=zeros(length(BFlist), 8);
for i=1:length(BFlist)
    BF=BFlist(i);
    for j=1:8
        % interpolation uses the BF in column 8 as the 'guide'
        channelParameters(i,j)= interp1q(interpolationTable(:,8), interpolationTable(:,j),  BF);
        
        if isnan (channelParameters(i,j))
            error('MAPparameters***: channel parameter interpolation has failed (BF too high or too low)')
        end
    end
end

% now save them
DRNLParams.channelParameters=	channelParameters;

% and display them
if isstruct(betweenRuns) && betweenRuns.runNumber==1 && length(BFlist)<2
    % 	disp('Fixed points used for DRNL interpolation')
    % 	disp(' nonlinBW	a    b    linCF	 linBW	 linGAIN	 c  nonlinCF')
    % 	UTIL_printTabTable(interpolationTable)
    % 	disp(' ')
    disp('DRNL parameters used in this run')
    disp(' nonlinBW	a    b    linCF	 linBW	 linGAIN	 c  nonlinCF')
    UTIL_printTabTable(channelParameters);
    disp(' ')
end

% write all parameters to the command window
if showParams
    disp('OME filters')
    disp(num2str(outerMiddleEarParams.filters))
    nm=UTIL_paramsList(whos);
    for i=1:length(nm)
        eval(['UTIL_showStruct(' nm{i} ', ''' nm{i} ''')'])
    end
end




% **********************************************************************  comparison data
% store individual data here for display on the multiThreshold GUI (if used)
% the final value in each vector is an identifier (BF or duration))
if isstruct(experiment)
    switch experiment.paradigm
        case 'IFMC'
            % based on CMa
            comparisonData=[
                48	36	34	29	27	31	35	250;
                46	34	21	19	20	21	28	500;
                51	28	23	21	28	41	62	1000;
                42	33	21	17	27	35	52	2000;
                56	40	22	17	22	32	48	4000;
                61	46	27	24	31	46	33	6000;
                74	69	31	30	41	47	45	8000
                ];
            if length(BFlist)==1 && ~isempty(comparisonData)
                availableFrequencies=comparisonData(:,end)';
                findRow= find(BFlist==availableFrequencies);
                if ~isempty (findRow)
                    experiment.comparisonData=comparisonData(findRow,:);
                end
            end
            
        case 'TMC'
            % based on CMa
            comparisonData=[
                35	42	56	77	87	92	110	110	110	250;
                28	41	45	56	53	73	67	79	73	500;
                27	32	42	49	58	72	78	82	84	1000;
                20	32	44	45	53	65	83	91	88	2000;
                21	26	30	44	67	85	89	94	102	4000;
                38	41	48	57	74	84	92	91	110	6000
                ];
            if length(BFlist)==1 && ~isempty(comparisonData)
                availableFrequencies=comparisonData(:,end)';
                findRow= find(BFlist==availableFrequencies);
                if ~isempty (findRow)
                    experiment.comparisonData=comparisonData(findRow,:);
                end
            end
            
        case { 'absThreshold',  'absThreshold_8'}
            % CMa thresholds
            experiment.comparisonData=[
                22	15	15	12	12	28 0.008;
                7	-3	-7	-5	3	12 0.500
                ];
            
        otherwise
            experiment.comparisonData=[];
    end
end


