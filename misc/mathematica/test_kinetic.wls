#!/ccc/products1/mathematica-10.0.0/RedHat-7-x86_64/bin/WolframScript -script

(*==============================================================================*)
(* HFB3                                                                         *)
(* Copyright CEA, DAM F-91297 Arpajon, France                                   *)
(*                                                                              *)
(* This program is free software: you can redistribute it and/or modify         *)
(* it under the terms of the GNU General Public License as published by         *)
(* the Free Software Foundation, either version 3 of the License, or            *)
(* (at your option) any later version.                                          *)
(*                                                                              *)
(* This program is distributed in the hope that it will be useful,              *)
(* but WITHOUT ANY WARRANTY; without even the implied warranty of               *)
(* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *)
(* GNU General Public License for more details.                                 *)
(*                                                                              *)
(* You should have received a copy of the GNU General Public License            *)
(* along with this program.  If not, see <http://www.gnu.org/licenses/>.        *)
(*==============================================================================*)

ClearAll["Global`*"];
zpart[z_, nz_, bz_] := 1/Sqrt[bz*Sqrt[Pi]]*1/Sqrt[2^nz*nz!]*Exp[-z^2/2/bz^2]* HermiteH[nz, z/bz];
rpart[r_, m_, n_, br_] := 1/(br*Sqrt[Pi])*Sqrt[n!/(n + m)!]*Exp[-r^2/2/br^2]*(r/br)^m* LaguerreL[n, m, r^2/br^2];
mpart[\[Theta]_, m_] := Exp[I*m*\[Theta]];

totfunc[r_, \[Theta]_, z_, m_, n_, nz_, d_, br_, bz_, d0_] := rpart[r, m, n, br]*zpart[z + (d - 1/2)*d0, nz, bz]* mpart[\[Theta], m];

br = 2.3878704919130063;
bz = 3.1087425890539850;
d0 = 5.0;

hbarc = 197.32698045930246;
nmass = hbarc * hbarc / 41.47;
(* factor = hbarc * hbarc / nmass; *)

massP = 938.27;
massN = 939.57;

factor = hbarc * hbarc / massP;

iz = 20.0;
in = 22.0;

funcBra[r_, \[Theta]_, z_] := Conjugate[totfunc[r, \[Theta], z, 2, 0, 5, 1, br, bz, d0]];
funcKet[r_, \[Theta]_, z_] :=           totfunc[r, \[Theta], z, 1, 2, 0, 0, br, bz, d0];
Print[-factor * (0.5 - 0.5 / (in + iz)) * Integrate[ funcBra[r, \[Theta], z]* Laplacian[funcKet[r, \[Theta], z], {r, \[Theta], z}, "Cylindrical"]*r, {\[Theta], 0, 2*Pi}, {r, 0, Infinity}, {z, -Infinity, Infinity}]]

funcBra[r_, \[Theta]_, z_] := Conjugate[totfunc[r, \[Theta], z, 2, 0, 2, 0, br, bz, d0]];
funcKet[r_, \[Theta]_, z_] :=           totfunc[r, \[Theta], z, 2, 0, 4, 1, br, bz, d0];
Print[-factor * (0.5 - 0.5 / (in + iz)) * Integrate[ funcBra[r, \[Theta], z]* Laplacian[funcKet[r, \[Theta], z], {r, \[Theta], z}, "Cylindrical"]*r, {\[Theta], 0, 2*Pi}, {r, 0, Infinity}, {z, -Infinity, Infinity}]]
