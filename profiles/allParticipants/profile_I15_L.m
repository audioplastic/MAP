function x = profile_MFR_L									
x.Comments= {									
	};								
x.BFs= [   									
	250	500	1000	2000	4000	6000	8000		
	];								
x.LongTone= [ 									
	48.35	49.70	51.05	60.25	54.25	71.33	74.15		
	];								
x.ShortTone= [ 									
	59.35	59.40	58.15	65.40	60.40	76.20	79.40		
	];								
x.IFMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.MaskerRatio=[    									
	0.5	0.7	0.9	1	1.1	1.3	1.6		
	];								
x.IFMCs=[									
	NaN	63.85	65.64	59.52	65.64	64.07	NaN		
	NaN	59.95	62.71	62.12	62.06	69.27	NaN		
	NaN	56.70	61.58	62.12	57.68	72.08	NaN		
	NaN	58.65	60.60	64.93	59.63	72.08	NaN		
	NaN	55.73	63.69	66.02	65.96	68.62	NaN		
	NaN	59.95	63.20	62.33	69.70	89.53	NaN		
	NaN	60.60	67.59	58.00	66.94	75.12	NaN		
	];								
x.IFMCs= x.IFMCs';									
x.Gaps= [									
	0.01	0.02	0.03	0.04	0.05	0.06	0.07	0.08	0.09
	];								
x.TMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.TMC= [									
	NaN	NaN	57.60	60.40	NaN	78.30	NaN		
	NaN	NaN	57.10	61.30	NaN	79.40	NaN		
	NaN	NaN	53.90	62.10	NaN	84.90	NaN		
	NaN	NaN	54.10	59.80	NaN	74.20	NaN		
	NaN	NaN	56.30	61.50	NaN	NaN	NaN		
	NaN	NaN	56.70	61.20	NaN	NaN	NaN		
	NaN	NaN	51.20	62.80	NaN	NaN	NaN		
	NaN	NaN	52.60	60.60	NaN	NaN	NaN		
	NaN	NaN	52.90	60.30	NaN	NaN	NaN		
	];								
x.TMC = x.TMC';									
