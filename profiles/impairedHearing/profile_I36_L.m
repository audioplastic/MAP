function x = profile_AFR_L									
x.Comments= {									
	};								
x.BFs= [   									
	250	500	1000	2000	4000	6000	8000		
	];								
x.LongTone= [ 									
	18.61	18.66	21.11	21.53	73.24	85.52	84.01		
	];								
x.ShortTone= [ 									
	27.32	24.56	24.37	28.16	77.62	89.33	85.11		
	];								
x.IFMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.MaskerRatio=[    									
	0.5	0.7	0.9	1	1.1	1.3	1.6		
	];								
x.IFMCs=[									
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	35.11	43.35	60.75	67.21	NaN	NaN		
	NaN	30.42	37.10	43.06	75.47	NaN	NaN		
	NaN	28.24	35.32	36.53	NaN	NaN	NaN		
	NaN	29.00	32.29	65.92	NaN	NaN	NaN		
	NaN	37.49	53.96	71.92	76.39	NaN	NaN		
	NaN	NaN	NaN	NaN	78.82	NaN	NaN		
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
	NaN	32.85	38.88	41.44	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	40.26	43.86	55.85	NaN	NaN	NaN		
	NaN	41.41	44.60	53.72	NaN	NaN	NaN		
	NaN	46.90	48.91	60.64	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	50.83	72.22	65.01	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	];								
x.TMC = x.TMC';									
