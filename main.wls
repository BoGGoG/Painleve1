<<"source/SeriesGenerator.wl"

orderN = 10;
coeffs = Painleve1PertSol[orderN];
(* want to find prefactor in growth of coefficients.
	Have to guess the Gamma *)
Atable = Table[ coeffs[[i]] / (-1)^(i+1) / Gamma[2 i - 1/2], {i,1,Length@coeffs}] // N;

anprefactor = RichardsonExtrapolate[Atable, 5]
Print["The analytic result is ", 1/Pi Sqrt[ 6 / (5 Pi) ] // N]

BorelTransformCoeffs[coeffs]

BT = BorelTransform[coeffs, p];
PBT = PadeApproximant[BT, {p, 0, {orderN-1, orderN}}];

poles = p/.NSolve[Denominator@PBT == 0, p];
zeroes = p/.NSolve[Numerator@PBT == 0, p];

rePlot = ReImPlot[PBT /. {p-> pp}, {pp,0,2},
	AxesLabel -> {"Re[p]", "Pade-Borel"}];
imPlot = ReImPlot[PBT /. {p-> I pp}, {pp,0,2},
	AxesLabel -> {"Im[p]", "Pade-Borel"}];
plot = Show[GraphicsRow[{rePlot, imPlot}],
	PlotLabel->"Pade-Borel Transform of P1 along real and imaginary axes"];
Export["plots/Pade-Borel-ReIm.png", plot];

poles = p/.NSolve[Denominator@PBT == 0, p];
zeroes = p/.NSolve[Numerator@PBT == 0, p];

polesPlot = ListPlot[Map[{Re[#], Im[#]}&, poles], PlotStyle->Red];
zeroesPlot = ListPlot[Map[{Re[#], Im[#]}&, zeroes], PlotStyle->Blue]
plot = Show[zeroesPlot, polesPlot, PlotLabel->"Fig. 1: Poles (red) and zeroes (blue) of the Pade-Borel transform"];
Export["plots/PadePolesAndZeroes.png", plot];

(* Laplace Transform *)
xList = Range[0,3.5,0.01];
hList = Map[InverseBorelTransform[PBT, p, #]&, xList];
yList = - Sqrt[xList / 6] ( 1 + hList ); (* eq. (4) *)

(* plot = ListPlot[Thread@{xList, hList}];*)
yPlot = ListPlot[Thread@{xList, yList},
	AxesLabel->{"x", "y(x)"},
	PlotLabel->"Fig. 2",
	Joined->True];

tritronquee[x_] := - 0.1875543083 - 0.3049055603 x + 0.1055298557 x^2 - 0.05229396374 x^3;
tritronqueePlot = Plot[tritronquee[x], {x,0,3.5}, PlotStyle->Gray];

Export["plots/yPlot.png", Show[yPlot, tritronqueePlot]];


(* Taylor Expansion *)
