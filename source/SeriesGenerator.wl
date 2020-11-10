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

Begin["`Private`"];

(* MAIN FUNCTIONS *)

Painleve1PertSol[order_] := Block[
	{DEQ, hSeries, aArr, DEQPert, leadingOrder, n},

	DEQ = h''[t] + 1/t h'[t] + h[t] ( 1 + 1/2 h[t] ) - 4/(25 t^2) ( 1 + h[t] ); (* == 0 *)
	hSeries[ord_] := Sum[a[n] / (#^(2 n) ), {n,1,ord}]&;
	aArr = {};

	leadingOrder = -2;
	Do[
		DEQPert = DEQ /. {h->hSeries[n+1]} /. aArr;
		LeadingOrderCoeff = Coefficient[DEQPert, t, leadingOrder-2n+2];
		aSol = Solve[LeadingOrderCoeff == 0, a[n]][[1,1]];
		AppendTo[aArr, aSol];
		,{n, 1, order}];

	Values@aArr

];

LeadingOrder[series_, t_] := Block[{},
	Exponent[series, t]
];

(* END OF FUNCTIONS *)

Scan[SetAttributes[#, {Protected, ReadProtected}]&,
     Select[Symbol /@ Names["SeriesGenerator`*"], Head[#] === Symbol &]];
End[];
EndPackage[];
