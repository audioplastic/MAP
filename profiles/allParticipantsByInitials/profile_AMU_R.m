function x = profile_AMU_R									
x.Comments= {									
	};								
x.BFs= [   									
	250	500	1000	2000	4000	6000	8000		
	];								
x.LongTone= [ 									
	13.37	20.02	47.16	45.05	55.08	71.11	79.13		
	];								
x.ShortTone= [ 									
	35.86	38.68	59.73	59.73	71.75	82.15	85.73		
	];								
x.IFMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.MaskerRatio=[    									
	0.5	0.7	0.9	1	1.1	1.3	1.6		
	];								
x.IFMCs=[									
	45.10	54.40	83.01	86.95	63.58	NaN	NaN		
	41.85	49.87	76.39	81.00	62.65	NaN	NaN		
	34.24	47.74	70.52	75.64	69.97	NaN	NaN		
	28.66	50.37	67.84	73.07	71.67	NaN	NaN		
	26.85	52.35	71.94	73.44	75.67	NaN	NaN		
	35.63	61.21	81.89	83.90	89.31	NaN	NaN		
	57.76	71.51	92.14	84.36	NaN	NaN	NaN		
	];								
x.IFMCs= x.IFMCs';									
x.Gaps= [									
	0.01	0.02	0.03	0.04	0.05	0.06	0.07	0.08	0.09
	];								
x.TMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.TMC= [									
	NaN	NaN	72.45	69.85	NaN	90.37	NaN		
	NaN	NaN	73.73	71.58	82.73	92.18	NaN		
	NaN	NaN	79.68	78.85	87.64	NaN	NaN		
	NaN	NaN	80.52	76.99	85.11	NaN	NaN		
	NaN	NaN	80.32	81.43	91.82	NaN	NaN		
	NaN	NaN	NaN	85.31	90.05	NaN	NaN		
	NaN	NaN	NaN	85.54	89.37	NaN	NaN		
	NaN	NaN	NaN	84.63	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	88.94	NaN	NaN		
	];								
x.TMC = x.TMC';									
