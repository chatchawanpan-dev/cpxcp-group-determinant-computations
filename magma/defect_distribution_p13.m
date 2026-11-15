p := 13; B := StringToInteger(B);
assert B ge 2;
SetSeed(1);
SetColumns(1024);
print "Magma version:", GetVersion();
print "input p, B, seed:", p, B, 1;
K<z> := CyclotomicField(p);
O := RingOfIntegers(K);
assert ClassNumber(O : Proof := "Full", Bound := MinkowskiBound(O)) eq 1;
print "certified class number:", 1;
if assigned records then
  recordfile := Open(records, "w");
  fprintf recordfile, "ell;norm;power_basis_coefficients;higher_odd_coordinates;admitted\n";
end if;
F := GF(p);
PS<T> := PowerSeriesRing(F, p-1);
E := &+[ (F!Factorial(n))^-1 * T^n : n in [0..p-2] ];
Epow := [PS!1];
for t in [1..p-2] do Append(~Epow, Epow[t]*E); end for;
function residue(y)
  c := Eltseq(K!y);
  return &+[ (F!c[t+1]) * Epow[t+1] : t in [0..p-2] ];
end function;
function logs(y)
  r := residue(y);
  c0 := Coefficient(r, 0);
  X := r/c0 - 1;
  L := PS!0; P := PS!1;
  for k in [1..p-2] do
    P := P*X;
    L := L + (F!((-1)^(k+1)) / F!k) * P;
  end for;
  return [ Integers()!Coefficient(L, j) : j in [1..p-2] ];
end function;
H := AssociativeArray();
n := 0; adm := [];
for ell in PrimesUpTo(B) do
  if ell mod p ne 1 then continue; end if;
  I := Decomposition(O, ell)[1][1];
  ok, y := IsPrincipal(I);
  assert ok and Abs(Norm(y)) eq ell;
  L := logs(y);
  n +:= 1;
  key := L[3] + 13*L[5] + 169*L[7];
  if IsDefined(H, key) then H[key] +:= 1; else H[key] := 1; end if;
  accepted := L[3] eq 0 and L[5] eq 0 and L[7] eq 0 and (L[9] ne 0 or L[11] ne 0);
  if assigned records then
    coefficients := [Integers()!c : c in Eltseq(K!y)];
    fprintf recordfile, "%o;%o;%o;%o;%o\n", ell, Abs(Norm(y)), coefficients,
      [L[j] : j in [3..p-2 by 2]], accepted select 1 else 0;
  end if;
  if accepted then Append(~adm, ell); end if;
end for;
vals := [ IsDefined(H,k) select H[k] else 0 : k in [0..2196] ];
print "split", n, "cells", #Keys(H), "count(0,0,0)=", vals[1], "mean", RealField(6)!(n/2197);
print "count distribution (count:frequency):";
M := Max(vals); print [ <c, #[v : v in vals | v eq c]> : c in [0..M] ];
print "counts of (x,0,0), x=0..12:", [vals[1+x] : x in [0..12]];
print "counts of (0,y,0):", [vals[1+13*y] : y in [0..12]];
print "counts of (0,0,w):", [vals[1+169*w] : w in [0..12]];
print "admitted:", #adm, adm[1..Min(#adm,10)];
print "all admitted:", adm;
print "full tuple histogram by key x+13*y+169*z:", vals;
for checkpoint in [10^7, B] do
  if checkpoint le B then
    print "membership checkpoint:", checkpoint, #[ell : ell in adm | ell le checkpoint];
  end if;
end for;
if assigned records then delete recordfile; end if;
quit;
