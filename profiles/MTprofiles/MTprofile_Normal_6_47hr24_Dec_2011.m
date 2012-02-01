function x = MTprofile_Normal_6_47hr24_Dec_2011
%created: 6_47hr24_Dec_2011

x.BFs = [250   500  1000  2000  4000  8000];

x.LongTone = [8.42      4.56      3.04       5.8      13.5      17.7];
x.ShortTone = [12      10.9      7.03      11.8      14.7      25.1];

x.Gaps = [0.01      0.03      0.05      0.07      0.09];
x.TMCFreq = [250   500  1000  2000  4000  8000];
x.TMC = [
36.1	51.3	20.1	44.6	48.3	41.8	 
44.8	79.5	31.7	56.9	25.2	52.4	 
90.5	84.5	60.8	67.7	36.5	60.7	 
95.4	90.8	80.7	81.3	53.3	53.5	 
101	95.2	89	68.1	66.9	78	 
];
x.TMC = x.TMC';

x.MaskerRatio = [0.5      0.7      0.9        1      1.1      1.3      1.6];
x.IFMCFreq = [250   500  1000  2000  4000  8000];
x.IFMCs = [
63.3	79	72.6	72.8	69.5	79	 
48.3	59.8	48.6	68.1	58.3	83.7	 
24	43.8	22.9	45.2	19.4	37.2	 
24.1	37.3	18.6	24.4	29.3	33.2	 
22.4	47.5	29.6	51.5	33.9	35.4	 
20.3	45.1	42.2	64.5	38.4	60.9	 
30.5	83.3	89.7	92.7	71.1	91.7	 
];
x.IFMCs = x.IFMCs';
