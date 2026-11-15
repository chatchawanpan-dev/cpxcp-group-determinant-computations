function OddCoordinates(p,K,O,ell)
 fac:=Factorization(ideal<O|ell>); ok,y:=IsPrincipal(fac[1][1]);
 if not ok then return false,[]; end if;
 F:=GF(p); S<T>:=PowerSeriesRing(F,p-1); E:=&+[T^k/(F!Factorial(k)):k in [0..p-2]];
 cs:=Eltseq(K!y); Y:=&+[(S!F!cs[i])*E^(i-1):i in [1..#cs]];
 V:=Y/Coefficient(Y,0)-1; L:=&+[(-1)^(k+1)*V^k/(F!k):k in [1..p-2]];
 return true,[Coefficient(L,j):j in [3..p-2 by 2]];
end function;
K<z>:=CyclotomicField(11); O:=MaximalOrder(K); hits:=[]; candidates:=0;
ell:=2;
while ell le 100000 do
 if ell mod 11 eq 1 then
  candidates+:=1;ok,c:=OddCoordinates(11,K,O,ell);
  if ok and c[1] eq 0 and c[2] eq 0 and c[3..4] ne [0,0] then Append(~hits,<ell,c>);end if;
 end if;
 ell:=NextPrime(ell);
end while;
print "p11_candidates",candidates,"hits",hits;
K<z>:=CyclotomicField(13);O:=MaximalOrder(K);ok,c:=OddCoordinates(13,K,O,6285241);print "p13_first_asserted",ok,c;
quit;
