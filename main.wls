<<"source/SeriesGenerator.wl"

orderN = 10;
coeffs = Painleve1PertSol[orderN];
(* want to find prefactor in growth of coefficients.
	Have to guess the Gamma *)
Atable = Table[ coeffs[[i]] / (-1)^(i+1) / Gamma[2 i - 1/2], {i,1,Length@coeffs}] // N;

Length@coeffs

anprefactor = RichardsonExtrapolate[Atable, 5]
Print["The analytic result is ", 1/Pi Sqrt[ 6 / (5 Pi) ] // N]

BT = BorelTransform[coeffs, p];
PBT = PadeApproximant[BT, {p, 0, {orderN-1, orderN}}];
PBT = Apart[(Numerator@PBT // N) / FactorCompletely[Denominator@PBT // N, p], p];

Length@BorelTransformCoeffs[coeffs]

poles = p/.NSolve[Denominator@PBT == 0, p];
zeroes = p/.NSolve[Numerator@PBT == 0, p];

rePlot = ReImPlot[PBT /. {p-> pp}, {pp,0,2},
	AxesLabel -> {"Re[p]", "Pade-Borel"}];
imPlot = ReImPlot[PBT /. {p-> I pp}, {pp,0,2},
	AxesLabel -> {"Im[p]", "Pade-Borel"}];
plot = Show[GraphicsRow[{rePlot, imPlot}],
	PlotLabel->"Pade-Borel Transform of P1 along real and imaginary axes"];
Export["plots/Pade-Borel-ReIm.png", plot];

poles = p/.NSolve[Denominator@Together@PBT == 0, p];
zeroes = p/.NSolve[Numerator@PBT == 0, p];

polesPlot = ListPlot[Map[{Re[#], Im[#]}&, poles], PlotStyle->Red];
zeroesPlot = ListPlot[Map[{Re[#], Im[#]}&, zeroes], PlotStyle->Blue, PlotMarkers->{Automatic, 0.02}]
plot = Show[zeroesPlot, polesPlot, PlotLabel->"Fig. 1:\n Poles (red) and zeroes (blue) of the Pade-Borel transform"];
Export["plots/PadePolesAndZeroes.png", plot];

(* Laplace Transform *)
xList = Range[0,3.5,0.001];
hList = ParallelMap[InverseBorelTransform[PBT, p, (24 #)^(5/4) / 30]&, xList];
yList = - Sqrt[xList / 6] ( 1 + hList ); (* eq. (4) *)


(* plot = ListPlot[Thread@{xList, hList}];*)
yPlot = ListPlot[Thread@{xList, yList},
	AxesLabel->{"x", "y(x)"},
	PlotLabel->"Fig. 2",
	Joined->True];

tritronquee[x_] := - 0.1875543083 - 0.3049055603 x + 0.1055298557 x^2 - 0.05229396374 x^3;
tritronqueePlot = Plot[tritronquee[x], {x,0,3.5}, PlotStyle->{Gray, Dashed}];

Export["plots/yPlot.png", Show[yPlot, tritronqueePlot]];

(* Conformal Borel *)
CBT = ConformalMap[BT, p, z];
TCBT = Normal@Series[CBT, {z, 0,2 orderN - 1}]; (* taylor conformal borel *)
PCBT = PadeApproximant[TCBT, {z, 0, {orderN-1, orderN}}];
PCBT = Apart[(Numerator@PCBT // N) / FactorCompletely[Denominator@PCBT // N, z], z];
CPCBT = InverseConformalMap[PCBT, z, p];


(* stokes constants, figure 8 *)
polesC = z/.NSolve[Denominator@Together@PCBT == 0, z];
polesPlotC = ListPlot[Map[{Re[#], Im[#]}&, polesC], PlotStyle->Red,
	AxesLabel->{"Re[z]", "Im[z]"},
	PlotLabel->"Fig 8:\n Poles of conformally mapped Borel transform"];
Export["plots/conformalPoles.png", polesPlotC];

Print["We get for half the Stokes constant at p=i: ", Re[(Sqrt[p - I] CPCBT) /. p -> 0.99999 I // N]];
Print["The correct value is: ", 1 / (2 Pi) Sqrt[3/5] // N];

(* Singularities of Borel Transform Stokes Constant*)
dottedLine = Plot[1 / (2 Pi) Sqrt[3/5], {pp, 0.9, 1}, PlotStyle->{Dashed, Gray},
	PlotLabels->"analytic result"];
singPlotBT = Plot[Re[(Sqrt[p - I] BT) /. p -> I pp], {pp, 0.9, 1.0}, PlotLabels->"Borel"];
singPlotPBT = Plot[Re[(Sqrt[p - I] PBT) /. p -> I pp], {pp, 0.9, 1.0}, PlotStyle->Red,
	PlotLabels->"Pade-Borel"];
singPlotPCBT = Plot[Re[(Sqrt[p - I] CPCBT) /. p -> I pp], {pp, 0.9, 1.0}, PlotStyle->Black,
	PlotLabels->"Pade-Conformal-Borel"];
Export["plots/singPlot.png", Show[dottedLine, singPlotPBT, singPlotBT, singPlotPCBT, PlotRange->All]];

(* z plane poles (conformal) mapped back to Borel plane p *)
polesCPCBT = (2 # / (1 - #^2)&) @ polesC;
polesPlotCPCBT = ListPlot[Map[{Re[#], Im[#]}&, polesCPCBT], PlotStyle->Red,
	AxesLabel->{"Re[z]", "Im[z]"},
	PlotLabel->"Fig 8:\n Poles of CPCBT",
	PlotRange->{{-7,7}, {-4,4}}];
Export["plots/CPCBTPoles.png", polesPlotCPCBT];

(* Laplace transform of Pade-Conformal-Borel *)
hListC = ParallelMap[InverseBorelTransform[CPCBT, p, (24 #)^(5/4) / 30]&, xList];
yListC = - Sqrt[xList / 6] ( 1 + hListC ); (* eq. (4) *)


yPlotC = ListPlot[Thread@{xList, yListC},
	AxesLabel->{"x", "y(x)"},
	Joined->True,
	PlotStyle->Red,
	PlotLabels->"Pade-Comformal-Borel"];

Export["plots/yPlot2.png", Show[yPlot, tritronqueePlot, yPlotC, PlotRange->{{0,0.3}, {-0.17,-0.27}}]];

Export["plots/plot.png", Show[yPlot, yPlotC, tritronqueePlot, PlotRange->{{0,0.1},{-0.22, -0.17}}]];
