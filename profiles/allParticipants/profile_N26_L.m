function x = profile_DOE_L									
x.Comments= {									
	};								
x.BFs= [   									
	250	500	1000	2000	4000	6000	8000		
	];								
x.LongTone= [ 									
	10.39	3.30	-2.35	3.85	-2.69	8.29	22.33		
	];								
x.ShortTone= [ 									
	18.90	12.62	13.73	10.72	4.86	16.00	28.82		
	];								
x.IFMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.MaskerRatio=[    									
	0.5	0.7	0.9	1	1.1	1.3	1.6		
	];								
x.IFMCs=[									
	41.57	34.99	63.71	75.61	69.36	62.25	NaN		
	35.84	34.81	55.39	63.58	52.54	31.88	NaN		
	27.29	31.40	28.47	37.96	28.55	19.85	NaN		
	28.94	23.39	22.54	22.81	16.15	23.31	NaN		
	29.90	25.98	28.36	32.23	23.64	40.90	NaN		
	25.14	39.03	41.90	36.83	25.34	44.99	NaN		
	37.25	50.51	60.28	46.44	31.03	30.81	NaN		
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
	35.90	35.61	27.77	18.15	16.09	31.08	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	42.91	42.35	35.39	18.80	20.76	32.55	NaN		
	56.93	62.40	49.51	29.84	25.60	52.39	NaN		
	64.55	86.61	46.83	41.71	45.67	50.52	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	57.20	87.90	61.76	48.92	45.82	68.02	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	];								
x.TMC = x.TMC';									
