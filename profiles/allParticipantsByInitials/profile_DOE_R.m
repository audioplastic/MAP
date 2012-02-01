function x = profile_DOE_R									
x.Comments= {									
	};								
x.BFs= [   									
	250	500	1000	2000	4000	6000	8000		
	];								
x.LongTone= [ 									
	10.80	7.45	-2.44	-1.16	-2.09	5.69	21.61		
	];								
x.ShortTone= [ 									
	19.92	15.71	14.53	8.77	6.96	15.51	27.06		
	];								
x.IFMCFreq= [									
	250	500	1000	2000	4000	6000	8000		
	];								
x.MaskerRatio=[    									
	0.5	0.7	0.9	1	1.1	1.3	1.6		
	];								
x.IFMCs=[									
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
	NaN	NaN	NaN	NaN	NaN	NaN	NaN		
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
