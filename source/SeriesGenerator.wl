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

DEQPerturbativeSolution[DEQ_, t_, {h_, hSeries_, a_, order_}] := Block[{DEQPert, leadingOrder},
	DEQPert = DEQ /. {h->hSeries[order]};
	leadingOrder = LeadingOrder[DEQPert, t];
	LeadingOrderCoeff = Coefficient[DEQPert, t, leadingOrder];
	aSol = a[1] /. Solve[LeadingOrderCoeff == 0, a[1]]

];

LeadingOrder[series_, t_] := Block[{},
	Exponent[series, t]
];

(* END OF FUNCTIONS *)

Scan[SetAttributes[#, {Protected, ReadProtected}]&,
     Select[Symbol /@ Names["SeriesGenerator`*"], Head[#] === Symbol &]];
End[];
EndPackage[];
