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
DEQPerturbativeSolution::usage = "DEQPerturbativeSolution[DEQ, t, {h, hSeries, order}] returns the coefficients gotten from plugging in the hSeries for h in DEQ and solving order by order. The series is supposed to be for large t, so a series in 1/t";
FindLowestPower::usage = "FindLowestPower[series, t] finds the lowest power of t in the series";

Begin["`Private`"];

(* MAIN FUNCTIONS *)

DEQPerturbativeSolution[DEQ_, t_, {h_, hSeries_, a_, order_}] := Block[
	{DEQPert, leadingOrder, aArr, n},
	aArr = {};

	(* leadingOrder = LeadingOrder[DEQPert, t]; *)
	leadingOrder = -2;
	n = 1;

	DEQPert = DEQ /. {h->hSeries[n+1]};
	LeadingOrderCoeff = Coefficient[DEQPert, t, leadingOrder];
	aSol = Solve[LeadingOrderCoeff == 0, a[n]];
	AppendTo[aArr, aSol];

	n = 2;
	DEQPert = DEQ /. {h->hSeries[n+1]} /. aSol[[1]];
	LeadingOrderCoeff = Coefficient[DEQPert, t, leadingOrder-2n+2];
	aSol = a[n] /. Solve[LeadingOrderCoeff == 0, a[n]]

];

LeadingOrder[series_, t_] := Block[{},
	Exponent[series, t]
];

(* END OF FUNCTIONS *)

Scan[SetAttributes[#, {Protected, ReadProtected}]&,
     Select[Symbol /@ Names["SeriesGenerator`*"], Head[#] === Symbol &]];
End[];
EndPackage[];
