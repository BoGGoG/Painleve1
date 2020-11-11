<<"source/SeriesGenerator.wl"

coeffs = Painleve1PertSol[10];
(* want to find prefactor in growth of coefficients.
	Have to guess the Gamma *)
Atable = Table[ coeffs[[i]] / (-1)^(i+1) / Gamma[2 i - 1/2], {i,1,Length@coeffs}] // N;

Atable2 = Map[ RichardsonExtrapolate[Atable, #]&, Range[1,9]] // N


Print["The analytic result is ", 1/Pi Sqrt[ 6 / (5 Pi) ] // N]
