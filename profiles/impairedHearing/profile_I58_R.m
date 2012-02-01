function x = profile_BPR_R									
x.Comments= {									
	};								
x.BFs= [   									
	250	500	1000	2000	4000	6000	8000		
	];								
x.LongTone= [ 									
	54.16	56.47	59.50	53.87	60.38	76.62	NaN		
	];								
x.ShortTone= [ 									
	58.27	59.38	60.69	57.19	58.44	75.59	NaN		
	];								
x.IFMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.MaskerRatio=[    									
	0.5	0.7	0.9	1	1.1	1.3	1.6		
	];								
x.IFMCs=[									
	NaN	87.05	85.40	86.58	61.32	65.41	NaN		
	NaN	82.82	81.20	72.02	63.57	72.16	NaN		
	69.81	73.85	73.70	61.78	60.03	81.26	NaN		
	70.50	71.64	74.41	63.07	64.23	77.85	NaN		
	73.07	69.91	70.82	67.07	67.65	77.98	NaN		
	71.56	73.80	72.15	66.02	70.31	94.88	NaN		
	79.71	78.78	71.99	60.81	75.40	NaN	NaN		
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
	65.97	72.94	73.00	63.85	72.86	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	67.74	73.39	72.97	66.88	74.62	NaN	NaN		
	69.35	81.34	74.01	67.35	76.42	NaN	NaN		
	68.10	76.75	76.37	66.85	74.45	81.02	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	69.77	85.06	78.17	67.39	63.46	83.97	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	];								
x.TMC = x.TMC';									
