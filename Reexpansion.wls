#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Reexpansion*)


<<"NumericalCalculus`"


(* ::Section:: *)
(*Explanation*)


PIDEQ = y''[x] == 6 y[x]^2 - x


D[PIDEQ, x]


yn[x_, n_/;n>1] := D[PIDEQ, {x,n-2}][[2]]


yn[x, 3]


yn[x,2]


yn[x,3]


yn[x,4]


(* ::Section:: *)
(*Development*)


Get@FileNameJoin[{NotebookDirectory[], "source", "SeriesGenerator.wl"}]


(* ::Text:: *)
(*Goal: Express all y^(n>2) in terms of y, y'*)


(*ReexpandP1Analytic[{y[x], y[x]}, {y'[x], y'[x]}, {x, 3}, 5]*)


(* ::Section:: *)
(*Implementation*)


(* ::Subsection:: *)
(*Reexpansion of Pade-Borel*)


Get@FileNameJoin[{NotebookDirectory[], "source", "SeriesGenerator.wl"}]


orderN = 50;
oderReexp = 2 orderN; (* initial expansion was in 1/t^2, new expansion in 1/x *)
coeffs = Painleve1PertSol[orderN];


BT = BorelTransform[coeffs, p];
PBT = PadeApproximant[BT, {p, 0, {orderN-1, orderN}}];
(*PBT = Apart[(Numerator@PBT) / FactorCompletely[Denominator@PBT, p], p];*)


yy[x_, h_] = - Sqrt[(x/6)] ( 1 + h[t[x]] );
tt[x_] = (24 x)^(5/4)/30;
dy[x_, h_, dh_] = D[yy[x, h], x] /. {h[t_]->h, h'[t_]->dh, t -> tt};
ddy[x_, h_, dh_, ddh_] = D[yy[x, h], x,x] /. {h[t_]->h, h'[t_]->dh, h''[t_]->ddh, t -> tt};


x0 = 3;
h0 = InverseBorelTransform[PBT, p, (24 x0)^(5/4) / 30];
dh0 = InverseBorelTransform[(-p)PBT, p, (24 x0)^(5/4) / 30];
ddh0 = InverseBorelTransform[(-p)^2PBT, p, (24 x0)^(5/4) / 30];
y0 = yy[x,h]/.{h[t_]->h} /. {x->x0, h->h0};
dy0 = dy[x0, h0, dh0];
(* ddy0 = ddy[x0, h0, dh0, ddh0]; *) (* for testing *)


h0


h0C


dh0


dh0C


yReexpanded[x_] = ReexpandP1[{y0, dy0}, {x, x0}, oderReexp];
yReexpandedPade[x_] = PadeApproximant[yReexpanded[x], {x, x0, {oderReexp/2-1, oderReexp/2}}];
yReexpandedPade[x_] = Apart[(Numerator@yReexpandedPade[x]) / FactorCompletely[Denominator@yReexpandedPade[x], x], x];


tritronquee[x_] := - 0.1875543083 - 0.3049055603 x + 0.1055298557 x^2 - 0.05229396374 x^3;
tritronqueePlot = Plot[tritronquee[x], {x,0,3.5}, PlotStyle->{Gray, Dashed}];


yReexpPlot = Plot[Re@yReexpanded[x], {x,0,7}, PlotStyle->{Dashed, Red}];
yReexpPadePlot = Plot[Re@yReexpandedPade[x], {x,0,7}, PlotStyle->Red];
yPlot = ListPlot[Thread@{xList, yList}, Joined->True, PlotStyle->{Blue, Dashed}];
Show[yReexpPlot, yReexpPadePlot, yPlot, tritronqueePlot, PlotRange->{{0,6}, {-1,0}}]


(* ::InheritFromParent:: *)
(**)


yReexpPlot = Plot[Re@yReexpanded[x], {x,0,0.2}, PlotStyle->{Dashed, Red}];
yReexpPadePlot = Plot[Re@yReexpandedPade[x], {x,0,1}, PlotStyle->Red];
yPlot = ListPlot[Thread@{xList, yList}, Joined->True, PlotStyle->{Blue, Dashed}];
Show[yReexpPadePlot, yPlot, tritronqueePlot]


(* ::Text:: *)
(*yReexpPade is pretty much spot on with Tritronquee!*)


yReexpanded[0]


yReexpandedPade[0]


(* ::Text:: *)
(*From the paper: \[Minus]0.187554308340*)


(* ::Subsection:: *)
(*Pade-Conformal-Borel*)


CBT = ConformalMap[BT, p, z];
TCBT = Normal@Chop@Series[CBT, {z, 0,2 orderN - 1}]; (* taylor conformal borel *)
PCBT = PadeApproximant[TCBT, {z, 0, {orderN-1, orderN}}];
(*PCBT = Apart[(Numerator@PCBT) / FactorCompletely[Denominator@PCBT, z], z];*)
CPCBT = InverseConformalMap[PCBT, z, p];


(*hListC = Map[InverseBorelTransform[CPCBT, p, (24 #)^(5/4) / 30]&, xList];
dhListC = Map[InverseBorelTransform[(- p) CPCBT, p, (24 #)^(5/4) / 30]&, xList];
ddhListC = Map[InverseBorelTransform[(- p)^2 CPCBT, p, (24 #)^(5/4) / 30]&, xList];*)


x0 = 3;
h0C = InverseBorelTransform[CPCBT, p, (24 x0)^(5/4) / 30];
dh0C = InverseBorelTransform[(-p)CPCBT, p, (24 x0)^(5/4) / 30];
ddh0C = InverseBorelTransform[(-p)^2 CPCBT, p, (24 x0)^(5/4) / 30];
y0C = yy[x,h]/.{h[t_]->h} /. {x->x0, h->h0C};
dy0C = dy[x0, h0C, dh0C];


yReexpandedC[x_] = ReexpandP1[{y0C, dy0C}, {x, x0}, oderReexp];
yReexpandedPadeC[x_] = PadeApproximant[yReexpandedC[x], {x, x0, {oderReexp/2-1,oderReexp/2}}];
yReexpandedPadeCApart[x_] = Chop@Apart[(Numerator@yReexpandedPadeC[x]) / FactorCompletely[Denominator@yReexpandedPadeC[x], x], x];


yReexpPlot = Plot[yReexpanded[x], {x,0,2}, PlotStyle->{Dashed, Red}];
yReexpPadePlot = Plot[Re@yReexpandedPade[x], {x,0,5}, PlotStyle->Red];
yReexpPadeCPlot = Plot[Re@yReexpandedPadeC[x], {x,0,5}, PlotStyle->{Green}];
yPlot = ListPlot[Thread@{xList, Re@yList}, Joined->True, PlotStyle->{Blue, Dashed}];
tritronquee[x_] := - 0.1875543083 - 0.3049055603 x + 0.1055298557 x^2 - 0.05229396374 x^3;
tritronqueePlot = Plot[tritronquee[x], {x,0,3.5}, PlotStyle->{Gray, Dashed, Thick}];

Show[yReexpPadeCPlot, yReexpPadePlot, tritronqueePlot, PlotRange->{{0,0.6}, {-0.35, -0.17}}];
Show[tritronqueePlot, yReexpPadeCPlot, yReexpPadePlot,yPlot, PlotRange->Automatic]


yReexpandedPadeC[0]


(* ::Section:: *)
(*Going to negative x*)


yReexpPlot = Plot[Re@yReexpanded[x], {x,-10, 10}, PlotStyle->{Dashed, Red}];
yReexpPadePlot = Plot[Re@yReexpandedPade[x], {x,-10, 10}, PlotStyle->Red];
yReexpPadeCPlot = Plot[Re@yReexpandedPadeC[x], {x,-10, 10}, PlotStyle->{Green}];
yPlot = ListPlot[Thread@{xList, yList}, Joined->True, PlotStyle->{Blue, Dashed}];
tritronquee[x_] := - 0.1875543083 - 0.3049055603 x + 0.1055298557 x^2 - 0.05229396374 x^3;
tritronqueePlot = Plot[tritronquee[x], {x,-10, 10}, PlotStyle->{Gray, Dashed}];

Show[yReexpPadePlot, yReexpPadeCPlot, yReexpPlot, tritronqueePlot, PlotRange->{{-10, 10}, {-2,5}}]


(* ::Text:: *)
(*For more terms we can get more than the first pole correct.*)


poles = x/.Solve[Denominator@yReexpandedPadeC[x] == 0, x];


(*poles = x/. Solve[1/yReexpandedPade[x]==0, x]*)


polesPlot = ListPlot[Map[{Re[#], Im[#]}&, poles], PlotStyle->Red, PlotRange->{{-12, 5}, {-8, 8}}]
