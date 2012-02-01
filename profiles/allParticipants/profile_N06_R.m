function x = profile_CMA_R									
x.Comments= {									
	};								
x.BFs= [   									
	250	500	1000	2000	4000	6000	8000		
	];								
x.LongTone= [ 									
	11.40	1.55	-13.50	-6.35	-6.40	NaN	7.45		
	];								
x.ShortTone= [ 									
	23.85	18.90	9.85	10.60	9.55	NaN	21.90		
	];								
x.IFMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.MaskerRatio=[    									
	0.5	0.7	0.9	1	1.1	1.3	1.6		
	];								
x.IFMCs=[									
	50.27	57.70	57.25	61.25	66.85	61.73	75.25		
	41.93	51.30	41.45	48.90	45.40	62.00	67.60		
	34.40	38.23	26.60	26.60	26.70	42.53	42.90		
	35.43	33.20	19.70	20.30	21.95	21.63	29.05		
	33.90	28.10	27.70	33.55	45.85	46.70	41.85		
	32.63	40.53	36.75	67.95	74.20	67.15	47.05		
	36.60	48.57	65.70	79.15	87.20	83.25	86.15		
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
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	];								
x.TMC = x.TMC';									
