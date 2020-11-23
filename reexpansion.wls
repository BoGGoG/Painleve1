<<"source/SeriesGenerator.wl"
(* Checking accuracy at x0 = 3 *)

orderN = 10;
coeffs = Painleve1PertSol[orderN];


(* Pade-Borel *)
BT = BorelTransform[coeffs, p];
PBT = PadeApproximant[BT, {p, 0, {orderN-1, orderN}}];
PBT = Apart[(Numerator@PBT // N) / FactorCompletely[Denominator@PBT // N, p], p];

xList = Range[3.,5.,1.];
hList = Map[InverseBorelTransform[PBT, p, (24 #)^(5/4) / 30]&, xList];
dhList = Map[InverseBorelTransform[(- p) PBT, p, (24 #)^(5/4) / 30]&, xList];
ddhList = Map[InverseBorelTransform[(- p)^2 PBT, p, (24 #)^(5/4) / 30]&, xList];

(* formulas for y, y', y'' in terms of h and x *)
yy[x_, h_] = - Sqrt[(x/6)] ( 1 + h[t[x]] );
tt[x_] = (24 x)^(5/4)/30;
dy[x_, h_, dh_] = D[yy[x, h], x] /. {h[t_]->h, h'[t_]->dh, t -> tt};
ddy[x_, h_, dh_, ddh_] = D[yy[x, h], x,x] /. {h[t_]->h, h'[t_]->dh, h''[t_]->ddh, t -> tt};

P1DEQ[x_, h_, dh_, ddh_] = ddy[x, h, dh, ddh] - 6 (yy[x, h]/.h[t_]->h)^2 + x;

P1DEQ[xList, hList, dhList, ddhList] (* checking how well my solutions fulfill the DEQ *)



(* Pade-Conformal-Borel *)
CBT = ConformalMap[BT, p, z];
TCBT = Normal@Series[CBT, {z, 0,2 orderN - 1}]; (* taylor conformal borel *)
PCBT = PadeApproximant[TCBT, {z, 0, {orderN-1, orderN}}];
PCBT = Apart[(Numerator@PCBT // N) / FactorCompletely[Denominator@PCBT // N, z], z];
CPCBT = InverseConformalMap[PCBT, z, p];

(* Laplace transform of Pade-Conformal-Borel *)
hListC = Map[InverseBorelTransform[CPCBT, p, (24 #)^(5/4) / 30]&, xList];
dhListC = Map[InverseBorelTransform[(- p) CPCBT, p, (24 #)^(5/4) / 30]&, xList];
ddhListC = Map[InverseBorelTransform[(- p)^2 CPCBT, p, (24 #)^(5/4) / 30]&, xList];
(* yListC = - Sqrt[xList / 6] ( 1 + hListC ); (* eq. (4) *) *)

P1DEQ[xList, hListC, dhListC, ddhListC] (* checking how well my solutions fulfill the DEQ *)
