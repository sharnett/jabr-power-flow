converged
function mpc = case5_renumber_tree

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 1;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	1	1	0	0	0	1	0.853327981	-27.8151106	345	1	1.1	0.9;
	2	1	1	0	0	0	1	0.984674135	-5.82880428	345	1	1.1	0.9;
	3	1	1	0	0	0	1	0.872948484	-20.100195	345	1	1.1	0.9;
	4	1	1	0	0	0	1	0.853327981	-27.8151106	345	1	1.1	0.9;
	5	3	0	0	0	0	1	1	0	345	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	5	4.1590465	1.590465	300	-300	1	1	1	10	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax	Pf	Qf	Pt	Qt
mpc.branch = [
	3	1	0.01	0.1	0	250	250	250	0	0	1	-360	360	1.0137	0.1373	-1.0000	0.0000;
	2	5	0.01	0.1	0	250	250	250	0	0	1	-360	360	-1.0000	-0.0000	1.0103	0.1031;
	3	4	0.01	0.1	0	250	250	250	0	0	1	-360	360	1.0137	0.1373	-1.0000	0.0000;
	5	3	0.01	0.1	0	250	250	250	0	0	1	-360	360	3.1487	1.4873	-3.0275	-0.2747;
];
