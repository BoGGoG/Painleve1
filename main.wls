<<"source/SeriesGenerator.wl"


DEQ = h''[t] + 1/t h'[t] + h[t] ( 1 + 1/2 h[t] ) - 4/(25 t^2) ( 1 + h[t] ); (* == 0 *)

hSeries[order_] := Sum[a[n] / (#^(2 n) ), {n,1,order}]&;

hSeries[5][t]

Collect[DEQ /. h->hSeries[5], t]

DEQPerturbativeSolution[DEQ, t, {h, hSeries, a, 3}]

Exponent[tmp, t]

Exponent[tmp, 1/t]

FindLowestPower[tmp, t]
