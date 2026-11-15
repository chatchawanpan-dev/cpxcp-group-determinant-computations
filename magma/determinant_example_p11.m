// End-to-end check of the prime-cofactor theorem at p=11: an explicit 121x121 integer
// group matrix of C_11 x C_11 with determinant 11^15 * ell for the least admitted prime ell.
p := 11; ell := StringToInteger(ell);
assert IsPrime(ell) and ell mod p eq 1;
SetSeed(1);
SetColumns(1024);
print "Magma version:", GetVersion();
print "input p, ell, seed:", p, ell, 1;
K<z> := CyclotomicField(p);
O := RingOfIntegers(K);
F := GF(p);
PS<T> := PowerSeriesRing(F, p-1);
E := &+[ (F!Factorial(n))^-1 * T^n : n in [0..p-2] ];
Epow := [PS!1];
for t in [1..p-2] do Append(~Epow, Epow[t]*E); end for;
function residue(y)
  c := Eltseq(K!y);
  return &+[ (F!c[t+1]) * Epow[t+1] : t in [0..p-2] ];
end function;
function LOG(r) // normalized log of a unit series r
  c0 := Coefficient(r, 0); X := r/c0 - 1; L := PS!0; P := PS!1;
  for k in [1..p-2] do P := P*X; L := L + (F!((-1)^(k+1))/F!k)*P; end for;
  return L;
end function;
m := (p-1) div 2;
U := [ O!(&+[z^k : k in [0..a-1]]) : a in [2..m] ];    // u_2..u_m
UL := [ LOG(residue(u)) : u in U ];
g := PrimitiveRoot(p);
function unit_with_residue(target) // target: series with nonzero constant, zero odd coords >=3
  L := LOG(target);
  for j in [3..p-2 by 2] do assert Coefficient(L, j) eq 0; end for;
  rows := [ [ Coefficient(UL[a], 2*i) : a in [1..m-1] ] : i in [1..m-1] ];
  A := Matrix(F, rows); rhs := Vector(F, [ Coefficient(L, 2*i) : i in [1..m-1] ]);
  sol := Solution(Transpose(A), rhs);
  n := [ Integers()!sol[a] : a in [1..m-1] ];
  u := &*[ U[a]^n[a] : a in [1..m-1] ];
  lin := Coefficient(LOG(residue(u)), 1);
  k := Integers()!(Coefficient(L,1) - lin);
  u := u * (O!z)^k;
  c := Coefficient(residue(u), 0); want := Coefficient(target, 0);
  j := Log(F!g, want/c);
  u := u * (U[g-1])^(p*j);
  assert residue(u) eq target;   // full residue modulo pR
  return u;
end function;
I := Decomposition(O, ell)[1][1];
ok, y := IsPrincipal(I); assert ok and Abs(Norm(y)) eq ell;
print "prime generator coefficients:", [Integers()!c : c in Eltseq(K!y)];
print "prime generator norm:", Norm(y);
d := [ Coefficient(LOG(residue(y)), j) : j in [3..p-2 by 2] ];   // d3,d5,d7,d9
print "defects:", d;
assert d[1] eq 0 and d[2] eq 0;
beta := d[3]; assert beta ne 0;      // use the exceptional packet (d_{p-4} != 0)
a := -beta^-1;
s := p;
Eu := &+[ (F!Factorial(n+1))^-1 * T^n : n in [0..p-2] ];  // (exp T - 1)/T
Einv := 1/Eu;
rays := AssociativeArray();
function ray(f, s)
  c := Eltseq(K!f); c := [Integers()!x : x in c];
  k := (s - &+c) div p; assert (s - &+c) mod p eq 0;
  return [ c[t] + k : t in [1..p-1] ] cat [k];
end function;
for b in [0..p-1] do
  target := -(1 + a*(F!b)^3*T^2) * Einv;          // z_b residue, constant -1
  ub := unit_with_residue(target);
  fb := (1 - O!z) * ub;
  assert Norm(fb) eq p;
  rays[b] := ray(fb, s);
end for;
// exceptional: y*u_inf == -a(1+beta T^7)/E^3 mod T^(p-3); extend by multiplying exact series
tgt := -a*(1 + beta*T^(p-4)) * Einv^3 / residue(y);
tgt := ChangePrecision(tgt, p-1);
// zero the odd coordinates beyond the visible precision p-1-3 = 7 (only degree 9 here)
L := LOG(tgt); 
Lfix := L - Coefficient(L, p-2)*T^(p-2);
tgt2 := Coefficient(tgt,0) * &+[ Lfix^k / Factorial(k) : k in [0..p-2] ];
uinf := unit_with_residue(tgt2);
finf := (1 - O!z)^3 * y * uinf;
print "norm f_inf / (p^3 ell) =", Norm(finf) / (p^3*ell);
rays[-1] := ray(finf, s);
v := [ [ 0 : j in [1..p] ] : i in [1..p] ];
for i in [0..p-1] do for j in [0..p-1] do
  num := rays[-1][j+1] + &+[ rays[b][((i + b*j) mod p)+1] : b in [0..p-1] ] - s;
  assert num mod p eq 0;
  v[i+1][j+1] := num div p;
end for; end for;
print "augmentation:", &+[ &+r : r in v ];
// Index (i,j) by p*i+j, with i major and j minor.
// Row w=(i,j), column x=(k,l): entry v[x-w].
Mx := ZeroMatrix(Integers(), p^2, p^2);
for i in [0..p-1] do
  for j in [0..p-1] do
    for k in [0..p-1] do
      for l in [0..p-1] do
        Mx[p*i+j+1,p*k+l+1] := v[((k-i) mod p)+1][((l-j) mod p)+1];
      end for;
    end for;
  end for;
end for;
D := Determinant(Mx);
print "generic character norms:", [ Norm(&+[v[i+1][j+1]*z^((i+b*j) mod p) : i,j in [0..p-1]]) : b in [0..p-1] ];
print "exceptional character norm:", Norm(&+[v[i+1][j+1]*z^j : i,j in [0..p-1]]);
print "determinant =", D;
print "equals 11^15 * ell:", D eq p^15*ell, " max |coefficient| =", Max([Abs(x) : x in &cat v]);
PrintFile("determinant_example_p11_vector.txt", v : Overwrite := true);
matrixfile := Open("determinant_example_p11_matrix_entries.csv", "w");
fprintf matrixfile, "row_index;column_index;value\n";
for i in [1..p^2] do
  for j in [1..p^2] do
    fprintf matrixfile, "%o;%o;%o\n", i-1, j-1, Mx[i,j];
  end for;
end for;
delete matrixfile;
quit;
