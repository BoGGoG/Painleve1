(* ::Package:: *)

(* Package: ChebyshevSolver
Author: Marco Knipfer
Email: mknipfer+seriesgenerator @ crimson . ua . edu
Institution: University of Alabama
Date: 11/2020
Description:
Usage:
Example:
*)

BeginPackage["SeriesGenerator`"];

Unprotect["SeriesGenerator`*"];
ClearAll["SeriesGenerator`*", "SeriesGenerator`Private`*"];

(* Public Functions, Documentations *)
Painleve1PertSol::usage = "Painleve1PertSol[DEQ, t, {h, hSeries, order}] returns the coefficients gotten from plugging in the hSeries for h in DEQ and solving order by order. The series is supposed to be for large t, so a series in 1/t";
FindLowestPower::usage = "FindLowestPower[series, t] finds the lowest power of t in the series";
RichardsonExtrapolate::usage = "RichardsonExtrapolate[series, orderN]: Extrapolate to the limit n->infty of the series a_n using Richardson's extrapolation on the last orderN elements.";
BorelTransformCoeffs::usage = "BorelTransformCoeffs[seriesCoeffs]: Transforms the coefficients a_n to a_n / (2n-1)!";
BorelTransform::usage = "BorelTransform[an, p]: Transforms the series a_n x^n to the Borel transformed a_n / (2n-1)! p^(2n-1)";
InverseBorelTransform::usage = "asdf";
ConformalMap::usage = "asdf";
InverseConformalMap::usage = "asdf";
FactorCompletely::usage = "FactorCompletely[polynoimal, var] does what built in Factor[...] cannot ..., found on https://mathematica.stackexchange.com/questions/8255/factoring-polynomials-to-factors-involving-complex-coefficients]";

Begin["`Private`"];

(* MAIN FUNCTIONS *)

Painleve1PertSol[order_] := Block[
	{DEQ, hSeries, aArr, DEQPert, leadingOrder, n},

	DEQ = h''[t] + 1/t h'[t] + h[t] ( 1 + 1/2 h[t] ) - 4/(25 t^2) ( 1 + h[t] ); (* == 0 *)
	hSeries[ord_] := Sum[a[n] / (#^(2 n) ), {n,1,ord}]&;
	aArr = {};

	leadingOrder = -2;
	Do[
		DEQPert = DEQ /. {h->hSeries[n+2]} /. aArr;
		LeadingOrderCoeff = Coefficient[DEQPert, t, leadingOrder-2n+2];
		aSol = Solve[LeadingOrderCoeff == 0, a[n]][[1,1]];
		AppendTo[aArr, aSol];
		,{n, 1, order}];

	Values@aArr

];

RichardsonExtrapolate[series_, orderN_] := Block[{n},
	If[orderN>Length@series, Print@"orderN is larger than the length of the series"];
	n = Length@series;
	Sum[ series[[n-orderN+k]] (n - orderN + k)^orderN (-1)^(k+orderN) / (k! (orderN - k)! ), {k, 0, orderN}]
];

LeadingOrder[series_, t_] := Block[{},
	Exponent[series, t]
];

BorelTransformCoeffs[seriesCoeffs_] := Block[{},
	Table[seriesCoeffs[[n]] / ((2 n - 1)!), {n,1,Length@seriesCoeffs}]
];

BorelTransform[seriesCoeffs_, p_] := Block[{borelCoeffs, n},
	borelCoeffs = BorelTransformCoeffs[seriesCoeffs];
	Sum[borelCoeffs[[n]] p^(2 n - 1), {n,1,Length@seriesCoeffs}]
];

InverseBorelTransform[borelTransform_, var_, newVar_] := Block[{},
	NIntegrate[Exp[-var newVar] borelTransform, {var, 0, Infinity}]
];

ConformalMap[borelTransform_, var_, newVar_] := Block[{},
	borelTransform /.var -> 2 newVar / ( 1 - newVar^2 )
]

InverseConformalMap[padeConformalBorelTranform_, var_, newVar_] := Block[{},
	padeConformalBorelTranform /. var -> newVar / (1 + Sqrt[1 + newVar^2])
];

FactorCompletely[poly_, x_] := Module[
  {solns, lcoeff},
  solns = Solve[poly == 0, x, Cubics -> False, Quartics -> False];
  lcoeff = Coefficient[poly, x^Exponent[poly, x]];
  lcoeff*(Times @@ (x - (x /. solns)))
  ]

(* END OF FUNCTIONS *)

Scan[SetAttributes[#, {Protected, ReadProtected}]&,
     Select[Symbol /@ Names["SeriesGenerator`*"], Head[#] === Symbol &]];
End[];
EndPackage[];
