function x = profile_BCR_L									
x.Comments= {									
	};								
x.BFs= [   									
	250	500	1000	2000	4000	6000	8000		
	];								
x.LongTone= [ 									
	15.52	13.50	7.46	33.46	58.47	67.30	78.65		
	];								
x.ShortTone= [ 									
	30.05	26.58	21.31	40.54	61.83	69.68	84.25		
	];								
x.IFMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.MaskerRatio=[    									
	0.5	0.7	0.9	1	1.1	1.3	1.6		
	];								
x.IFMCs=[									
	NaN	31.35	24.88	40.33	60.70	NaN	NaN		
	NaN	30.81	35.75	36.78	59.68	NaN	NaN		
	NaN	28.64	24.25	43.57	68.73	NaN	NaN		
	NaN	32.19	22.49	44.28	70.70	NaN	NaN		
	NaN	29.15	23.12	46.96	73.46	NaN	NaN		
	NaN	39.99	24.45	52.49	77.18	NaN	NaN		
	NaN	42.99	36.32	60.17	74.22	NaN	NaN		
	];								
x.IFMCs= x.IFMCs';									
x.Gaps= [									
	0.01	0.02	0.03	0.04	0.05	0.06	0.07	0.08	0.09
	];								
x.TMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.TMC= [									
	NaN	49.43	33.64	46.42	72.33	NaN	NaN		
	NaN	51.07	33.15	49.62	78.75	NaN	NaN		
	NaN	54.34	33.93	51.65	75.86	NaN	NaN		
	NaN	57.73	31.57	52.65	78.21	NaN	NaN		
	NaN	60.70	37.78	52.57	78.28	NaN	NaN		
	NaN	57.13	40.59	54.19	76.14	NaN	NaN		
	NaN	58.92	41.56	54.34	83.75	NaN	NaN		
	NaN	59.33	40.40	57.77	86.06	NaN	NaN		
	NaN	64.87	36.65	59.68	82.11	NaN	NaN		
	];								
x.TMC = x.TMC';									
