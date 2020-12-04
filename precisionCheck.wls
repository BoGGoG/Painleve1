#!/usr/bin/env wolframscript
(* ::Package:: *)

(*<<"source/SeriesGenerator.wl"*)
<<FileNameJoin[{NotebookDirectory[], "source", "SeriesGenerator.wl"}]
(* Checking accuracy at x0 = 3 *)

orderN = 50;
coeffs = Painleve1PertSol[orderN];
workingPrecision = 100;
precisionGoal = 70;


(* Pade-Borel *)
BT = BorelTransform[coeffs, p];
PBT = PadeApproximant[BT, {p, 0, {orderN-1, orderN}}];
PBT = Apart[(Numerator@PBT) / FactorCompletely[Denominator@PBT, p], p];

xList = {3} (*Range[3,5,1]*);
hList = Map[InverseBorelTransform[PBT, p, (24 #)^(5/4) / 30,
	PrecisionGoal->precisionGoal, WorkingPrecision -> workingPrecision]&, xList];
dhList = Map[InverseBorelTransform[(- p) PBT, p, (24 #)^(5/4) / 30,
	PrecisionGoal->precisionGoal, WorkingPrecision -> workingPrecision]&, xList];
ddhList = Map[InverseBorelTransform[(- p)^2 PBT, p, (24 #)^(5/4) / 30,
	PrecisionGoal->precisionGoal, WorkingPrecision -> workingPrecision]&, xList];

(* formulas for y, y', y'' in terms of h and x *)
yy[x_, h_] = - Sqrt[(x/6)] ( 1 + h[t[x]] );
tt[x_] = (24 x)^(5/4)/30;
dy[x_, h_, dh_] = D[yy[x, h], x] /. {h[t_]->h, h'[t_]->dh, t -> tt};
ddy[x_, h_, dh_, ddh_] = D[yy[x, h], x,x] /. {h[t_]->h, h'[t_]->dh, h''[t_]->ddh, t -> tt};

P1DEQ[x_, h_, dh_, ddh_] = ddy[x, h, dh, ddh] - 6 (yy[x, h]/.h[t_]->h)^2 + x;

(* checking how well my solutions fulfill the DEQ *)
padeBorelErrors = P1DEQ[xList, hList, dhList, ddhList];
Print["At x = ", xList, " the error from Pade-Borel are ", padeBorelErrors // N];

(*yy[x,h]/.{h[t_]->h} /. {x->xList, h->hList}*)

(* Pade-Conformal-Borel *)
CBT = ConformalMap[BT, p, z];
TCBT = Normal@Series[CBT, {z, 0,2 orderN - 1}]; (* taylor conformal borel *)
PCBT = PadeApproximant[TCBT, {z, 0, {orderN-1, orderN}}];
(*PCBT = Apart[(Numerator@PCBT) / FactorCompletely[Denominator@PCBT, z], z];*)
CPCBT = InverseConformalMap[PCBT, z, p];

(* Laplace transform of Pade-Conformal-Borel *)
hListC = Map[InverseBorelTransform[CPCBT, p, (24 #)^(5/4) / 30,
	PrecisionGoal->precisionGoal, WorkingPrecision -> workingPrecision]&, xList];
dhListC = Map[InverseBorelTransform[(- p) CPCBT, p, (24 #)^(5/4) / 30,

		PrecisionGoal->precisionGoal, WorkingPrecision -> workingPrecision]&, xList];
ddhListC = Map[InverseBorelTransform[(- p)^2 CPCBT, p, (24 #)^(5/4) / 30,
	PrecisionGoal->precisionGoal, WorkingPrecision -> workingPrecision]&, xList];
(* yListC = - Sqrt[xList / 6] ( 1 + hListC ); (* eq. (4) *) *)

(* checking how well my solutions fulfill the DEQ *)
padeConformalBorelErrors = P1DEQ[xList, hListC, dhListC, ddhListC] ;
Print["At x = ", xList, " the error from Pade-Conformal-Borel are ", padeConformalBorelErrors // N];


(* Directly going to the origin *)

xList2 = Range[1/1000, 1/8, 1/100];
hListC2 = InverseBorelTransform[CPCBT, p, (24 xList2)^(5/4) / 30];
dhListC2 = InverseBorelTransform[(- p) CPCBT, p, (24 xList2)^(5/4) / 30];
ddhListC2 = InverseBorelTransform[(- p)^2 CPCBT, p, (24 xList2)^(5/4) / 30];


yPlot2 = ListPlot[Thread@{xList2, yy[x,h]/.{h[t_]->h} /. {x->xList2, h->hListC2}},
	PlotRange->{-0.24, 0}];
Export["yPlot2.pdf", yPlot2]

Abs@P1DEQ[xList2, hListC2, dhListC2, ddhListC2] // N


yPlot2


errorPlot2 = ListPlot[Thread@{xList2, Log10@Abs@P1DEQ[xList2, hListC2, dhListC2, ddhListC2]}]
Export["errorPlot2.pdf", errorPlot2]

(* Directly going to the origin is not really working.
	Look at the yPlot2, which is basically Figure 6 in Dunnes
	paper. It cannot work.
	I guess what they did is do a linear regression from x=0.2 or
	so down to x=0.
	*)



