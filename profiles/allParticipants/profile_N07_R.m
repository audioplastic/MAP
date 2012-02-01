function x = profile_NBE_R									
x.Comments= {									
	};								
x.BFs= [   									
	250	500	1000	2000	4000	6000	8000		
	];								
x.LongTone= [ 									
	11.59	5.32	7.24	1.44	4.11	11.81	15.72		
	];								
x.ShortTone= [ 									
	26.86	22.68	17.68	17.37	18.05	20.28	26.62		
	];								
x.IFMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.MaskerRatio=[    									
	0.5	0.7	0.9	1	1.1	1.3	1.6		
	];								
x.IFMCs=[									
	NaN	52.95	31.57	40.85	60.61	69.63	NaN		
	NaN	42.74	28.53	53.30	49.52	35.01	NaN		
	NaN	36.30	21.84	29.50	35.91	35.95	NaN		
	NaN	33.88	25.63	16.48	23.57	21.32	NaN		
	NaN	34.81	28.86	27.59	36.34	45.03	NaN		
	NaN	41.81	35.51	44.46	41.82	51.63	NaN		
	NaN	50.71	59.55	44.62	53.90	37.23	NaN		
	];								
x.IFMCs= x.IFMCs';									
x.Gaps= [									
	0.01	0.02	0.03	0.04	0.05	0.06	0.07	0.08	0.09
	];								
x.TMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.TMC= [									
	NaN	33.87	37.00	21.42	25.77	36.75	NaN		
	NaN	39.80	40.09	33.52	29.16	37.80	NaN		
	NaN	33.73	42.39	29.14	26.00	43.66	NaN		
	NaN	38.74	47.73	29.02	32.21	43.76	NaN		
	NaN	51.04	56.18	39.75	35.04	47.06	NaN		
	NaN	64.33	58.55	31.65	40.92	55.17	NaN		
	NaN	78.52	75.17	60.42	57.40	70.78	NaN		
	NaN	87.36	78.39	76.62	75.58	85.41	NaN		
	NaN	90.20	86.42	81.20	84.93	89.51	NaN		
	];								
x.TMC = x.TMC';									
