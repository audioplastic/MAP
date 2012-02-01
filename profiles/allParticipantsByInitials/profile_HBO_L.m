function x = profile_HBO_L									
x.Comments= {									
	};								
x.BFs= [   									
	250	500	1000	2000	4000	6000	8000		
	];								
x.LongTone= [ 									
	33.45	29.73	28.98	37.38	34.16	31.94	52.13		
	];								
x.ShortTone= [ 									
	37.69	39.41	35.17	40.97	38.16	39.19	53.59		
	];								
x.IFMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.MaskerRatio=[    									
	0.5	0.7	0.9	1	1.1	1.3	1.6		
	];								
x.IFMCs=[									
	NaN	NaN	NaN	82.42	85.26	82.97	94.06		
	NaN	NaN	NaN	72.05	68.13	74.24	79.96		
	NaN	NaN	NaN	59.32	60.10	62.86	79.23		
	NaN	NaN	NaN	54.32	52.28	49.67	77.92		
	NaN	NaN	NaN	74.89	60.08	68.16	NaN		
	NaN	NaN	NaN	78.64	50.65	85.54	NaN		
	NaN	NaN	NaN	76.92	66.32	NaN	NaN		
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
	NaN	NaN	NaN	NaN	59.34	63.83	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	63.76	71.84	NaN		
	NaN	NaN	NaN	NaN	68.16	80.72	NaN		
	NaN	NaN	NaN	NaN	84.75	83.23	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	75.72	84.46	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	];								
x.TMC = x.TMC';									
