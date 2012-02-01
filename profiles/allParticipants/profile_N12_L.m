function x = profile_NCL_L									
x.Comments= {									
	};								
x.BFs= [   									
	250	500	1000	2000	4000	6000	8000		
	];								
x.LongTone= [ 									
	19.67	12.09	8.06	8.87	9.89	6.47	13.21		
	];								
x.ShortTone= [ 									
	27.18	17.81	13.24	12.48	15.40	14.29	22.07		
	];								
x.IFMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.MaskerRatio=[    									
	0.5	0.7	0.9	1	1.1	1.3	1.6		
	];								
x.IFMCs=[									
	47.34	51.94	48.26	70.38	71.78	63.75	NaN		
	39.10	39.79	41.13	67.19	53.91	47.96	NaN		
	40.35	33.61	29.57	37.80	31.32	21.96	NaN		
	36.71	33.28	28.59	25.94	26.83	21.06	NaN		
	37.16	29.36	46.27	46.10	22.22	37.29	NaN		
	39.52	31.40	59.80	50.79	18.32	31.69	NaN		
	45.49	84.68	53.03	63.79	52.48	51.58	NaN		
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
	47.13	37.78	35.51	45.41	26.48	25.99	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	61.76	41.93	54.52	63.95	35.38	38.43	NaN		
	61.81	48.74	63.21	66.17	42.39	40.72	NaN		
	64.89	58.67	77.04	87.71	55.48	67.98	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	93.33	67.45	90.67	98.00	59.55	80.13	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	];								
x.TMC = x.TMC';									
