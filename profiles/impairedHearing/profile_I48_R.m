function x = profile_SAL_R									
x.Comments= {									
	};								
x.BFs= [   									
	250	500	1000	2000	4000	6000	8000		
	];								
x.LongTone= [ 									
	27.74	28.04	14.85	14.67	22.47	45.56	55.91		
	];								
x.ShortTone= [ 									
	35.69	34.35	18.92	20.47	26.11	48.51	56.64		
	];								
x.IFMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.MaskerRatio=[    									
	0.5	0.7	0.9	1	1.1	1.3	1.6		
	];								
x.IFMCs=[									
	NaN	NaN	NaN	75.13	73.99	76.20	85.72		
	NaN	NaN	NaN	64.79	61.79	79.06	85.05		
	NaN	NaN	NaN	40.98	37.55	53.77	55.88		
	NaN	NaN	NaN	26.27	27.52	50.73	57.85		
	NaN	NaN	NaN	46.84	44.96	52.37	70.42		
	NaN	NaN	NaN	75.87	58.46	76.36	89.85		
	NaN	NaN	NaN	NaN	62.93	NaN	NaN		
	];								
x.IFMCs= x.IFMCs';									
x.Gaps= [									
	0.01	0.02	0.03	0.04	0.05	0.06	0.07	0.08	0.09
	];								
x.TMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.TMC= [									
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	27.40	31.06	50.96	59.23		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	31.44	42.80	56.10	63.39		
	NaN	NaN	NaN	32.27	65.40	61.46	65.82		
	NaN	NaN	NaN	38.74	52.66	72.09	66.02		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	41.63	NaN	NaN	70.72		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	];								
x.TMC = x.TMC';									
