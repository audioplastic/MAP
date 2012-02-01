function x = profile_AMO_R									
x.Comments= {									
	};								
x.BFs= [   									
	250	500	1000	2000	4000	6000	8000		
	];								
x.LongTone= [ 									
	16.42	17.71	5.67	10.12	51.92	69.85	79.29		
	];								
x.ShortTone= [ 									
	27.84	26.71	16.65	18.15	56.06	68.43	83.06		
	];								
x.IFMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.MaskerRatio=[    									
	0.5	0.7	0.9	1	1.1	1.3	1.6		
	];								
x.IFMCs=[									
	NaN	NaN	NaN	NaN	73.34	77.39	NaN		
	30.37	32.02	45.59	47.32	59.90	74.23	NaN		
	29.98	35.79	32.30	30.69	58.81	74.00	NaN		
	29.24	29.66	17.27	18.67	62.54	72.29	NaN		
	27.27	24.10	13.61	18.29	65.52	71.02	NaN		
	30.91	24.40	47.17	20.35	76.80	78.59	NaN		
	NaN	NaN	NaN	NaN	84.63	96.80	NaN		
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
	NaN	NaN	22.13	23.05	67.11	76.56	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	28.19	28.60	70.56	74.90	NaN		
	NaN	NaN	32.02	36.92	71.68	79.05	NaN		
	NaN	NaN	31.41	28.29	71.08	81.25	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	30.33	37.38	71.84	77.60	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	];								
x.TMC = x.TMC';									
