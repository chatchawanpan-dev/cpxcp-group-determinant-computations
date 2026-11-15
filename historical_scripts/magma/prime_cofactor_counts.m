// Admitted prime cofactors (valuation p+4) with counts at checkpoints.
p := StringToInteger(p); B := StringToInteger(B);
K<z> := CyclotomicField(p);
O := RingOfIntegers(K);
F := GF(p);
PS<T> := PowerSeriesRing(F, p-1);
Ex := &+[ (F!Factorial(n))^-1 * T^n : n in [0..p-2] ];   // exp(T)
Epow := [PS!1];
for t in [1..p-2] do Append(~Epow, Epow[t]*Ex); end for;
function residue(y)
  c := Eltseq(K!y);
  return &+[ (F!c[t+1]) * Epow[t+1] : t in [0..p-2] ];
end function;
function oddlogs(y)
  r := residue(y); c0 := Coefficient(r, 0); X := r/c0 - 1;
  L := PS!0; P := PS!1;
  for k in [1..p-2] do P := P*X; L := L + (F!((-1)^(k+1)) / F!k) * P; end for;
  return [ Coefficient(L, j) : j in [3..p-2 by 2] ];
end function;
h := 1;  // class number of Q(zeta_p) for p = 11, 13 (all ideals principal; checked below)
dens := (p+1)/(h*p^((p-3) div 2));
checkpoints := [10^k : k in [5..10] | 10^k le B] cat [B];
adm := []; nprimes := 0; nsplit := 0; ci := 1;
bound := p^(p-4)/2^(p-1);
for ell in PrimesUpTo(B) do
  nprimes +:= 1;
  while ci le #checkpoints and ell gt checkpoints[ci] do
    printf "x=%o pi(x)=%o admitted=%o predicted=%o\n", checkpoints[ci], nprimes-1, #adm, RealField(6)!(dens*(nprimes-1));
    ci +:= 1;
  end while;
  if ell mod p ne 1 then continue; end if;
  nsplit +:= 1;
  I := Decomposition(O, ell)[1][1];
  ok, y := IsPrincipal(I);
  assert ok and Abs(Norm(y)) eq ell;
  d := oddlogs(y); m := #d;
  if &and[ d[k] eq 0 : k in [1..m-2] ] and (d[m-1] ne 0 or d[m] ne 0) then
    Append(~adm, ell);
  end if;
end for;
printf "x=%o pi(x)=%o admitted=%o predicted=%o\n", B, nprimes, #adm, RealField(6)!(dens*nprimes);
printf "size bound p^(p-4)/2^(p-1) = %o ; least admitted = %o\n", RealField(8)!bound, #adm gt 0 select adm[1] else 0;
printf "first admitted: %o\n", adm[1..Min(#adm,12)];
quit;
