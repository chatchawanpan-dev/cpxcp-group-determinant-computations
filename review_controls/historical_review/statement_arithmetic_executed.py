import sympy as s
x=s.symbols('x');Y=-2-2*x-3*x*x-2*x**3-x**4-2*x**5-x**6-3*x**7
print('GEN1 resultant',s.resultant(s.cyclotomic_poly(11,x),Y,x),'prime',s.isprime(77419))
for p,x in [(11,100000),(11,1000000),(11,10000000),(11,30000000),(13,10000000),(13,40000000)]:
 n=int(s.primepi(x));print('density comparison',p,x,n,round(n*(p+1)/p**((p-3)//2),1))
for p in [11,13]:print('bounds',p,p**(p-4)/2**(p-1),p**(p-2)/2**(p-1))
