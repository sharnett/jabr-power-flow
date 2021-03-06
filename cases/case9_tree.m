function mpc = case9_tree

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 1;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	0	0	0	1	1	0	345	1	1.1	0.9;
	2	1	-163	-82.5	0	0	1	1	6.47166938	345	1	1.1	0.9;
	3	1	-85	-23.7825	0	0	1	1	8.7142251	345	1	1.1	0.9;
	4	1	0	0	0	0	1	0.956270817	-2.53581696	345	1	1.1	0.9;
	5	1	90	30	0	0	1	0.931054978	-2.74887062	345	1	1.1	0.9;
	6	1	0	0	0	0	1	0.987320689	5.82244453	345	1	1.1	0.9;
	7	1	100	35	0	0	1	0.914029602	-4.19987391	345	1	1.1	0.9;
	8	1	0	0	0	0	1	0.953893256	0.340835495	345	1	1.1	0.9;
	9	1	125	50	0	0	1	0.90987608	-5.8874998	345	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	73.4533988	77.5444488	300	-300	1	100	1	250	10	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax	Pf	Qf	Pt	Qt
mpc.branch = [
	1	4	0	0.000576	0	250	250	250	0	0	1	-360	360	73.4534	77.5444	-73.4534	-70.9731;
	8	2	0	0.000625	0	250	250	250	0	0	1	-360	360	-163.0000	-61.6404	163.0000	82.4999;
	3	6	0	0.000586	0	300	300	300	0	0	1	-360	360	85.0000	23.7825	-85.0000	-19.2172;
	4	5	0.00017	0.00092	0	250	250	250	0	0	1	-360	360	8.1642	24.7081	-8.0383	-24.0268;
	9	4	0.0001	0.00085	0	250	250	250	0	0	1	-360	360	-64.5890	-40.3133	65.2892	46.2651;
	5	6	0.00039	0.0017	0	150	150	150	0	0	1	-360	360	-81.9617	-5.9732	85.0000	19.2172;
	7	8	8.5e-05	0.00072	0	250	250	250	0	0	1	-360	360	-100.0000	-35.0000	101.1420	44.6738;
	8	9	0.00032	0.00161	0	250	250	250	0	0	1	-360	360	61.8580	16.9665	-60.4110	-9.6867;
];

%%-----  OPF Data  -----%%
%% area data
%	area	refbus
mpc.areas = [
	1	5;
];

%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	1500	0	3	0.11	5	150;
	2	2000	0	3	0.085	1.2	600;
	2	3000	0	3	0.1225	1	335;
];
