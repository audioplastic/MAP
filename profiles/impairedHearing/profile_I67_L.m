function x = profile_RME_L									
x.Comments= {									
	};								
x.BFs= [   									
	250	500	1000	2000	4000	6000	8000		
	];								
x.LongTone= [ 									
	17.15	17.67	7.73	17.58	36.03	NaN	76.63		
	];								
x.ShortTone= [ 									
	25.37	23.85	16.77	20.75	39.71	NaN	80.87		
	];								
x.IFMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.MaskerRatio=[    									
	0.5	0.7	0.9	1	1.1	1.3	1.6		
	];								
x.IFMCs=[									
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	36.01	41.82	40.94	61.71	NaN	NaN		
	NaN	28.64	28.86	35.68	52.93	NaN	NaN		
	NaN	34.98	26.72	31.57	50.52	NaN	NaN		
	NaN	28.27	29.38	40.90	54.94	NaN	NaN		
	NaN	28.32	24.89	53.43	60.00	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
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
	38.72	47.20	32.45	35.55	56.06	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	61.37	60.54	51.04	45.64	71.34	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	76.20	69.58	67.71	57.55	77.49	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	];								
x.TMC = x.TMC';									
