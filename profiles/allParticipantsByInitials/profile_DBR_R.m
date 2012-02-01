function x = profile_DBR_R									
x.Comments= {									
	};								
x.BFs= [   									
	250	500	1000	2000	4000	6000	8000		
	];								
x.LongTone= [ 									
	50.09	33.20	11.74	3.40	12.85	32.75	23.68		
	];								
x.ShortTone= [ 									
	55.44	39.98	17.26	11.68	22.10	37.47	34.39		
	];								
x.IFMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.MaskerRatio=[    									
	0.5	0.7	0.9	1	1.1	1.3	1.6		
	];								
x.IFMCs=[									
	76.14	51.27	NaN	77.49	72.39	72.86	81.29		
	65.89	45.11	NaN	58.37	61.66	70.91	69.29		
	59.19	34.70	NaN	32.07	38.55	53.58	45.33		
	57.51	31.81	NaN	20.34	29.59	47.26	36.19		
	57.41	27.78	NaN	40.99	60.24	37.48	88.03		
	53.71	42.14	NaN	83.26	85.81	46.99	74.36		
	46.69	88.76	NaN	86.46	NaN	NaN	NaN		
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
	76.38	45.60	NaN	30.36	43.65	59.13	45.16		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	81.43	56.90	NaN	38.86	49.46	69.79	58.93		
	NaN	60.36	NaN	42.82	57.15	68.44	60.37		
	NaN	66.77	NaN	47.43	51.85	72.18	65.49		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	75.78	NaN	59.91	NaN	NaN	61.63		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	];								
x.TMC = x.TMC';									
