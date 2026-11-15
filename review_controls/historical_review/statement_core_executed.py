import json,math,random,collections,sys,platform
from pathlib import Path
from flint import nmod_poly,fmpz_mat
print('environment',sys.version.split()[0],platform.platform(),'python-flint0.9.0')
PACK=Path('<review-packets>')
def mul(a,b,p,w):
 c=[0]*w
 for i,x in enumerate(a):
  for j,y in enumerate(b):
   if i+j<w:c[i+j]=(c[i+j]+x*y)%p
 return c
def logseries(a,p,w):
 v=[x*pow(a[0],-1,p)%p for x in a]+[0]*w;v=v[:w];v[0]=0
 ans=[0]*w;power=[1]+[0]*(w-1)
 for k in range(1,w):
  power=mul(power,v,p,w);scale=pow(k,-1,p)*(1 if k%2 else -1)
  ans=[(x+scale*y)%p for x,y in zip(ans,power)]
 return ans
def epoly(coeff,p,w):
 invfac=[pow(math.factorial(n),-1,p) for n in range(w)]
 return [sum(c*pow(t,n,p) for t,c in enumerate(coeff))*invfac[n]%p for n in range(w)]
def reduce_phi(f,p):
 f=list(f)
 while len(f)>=p:
  c=f.pop();start=len(f)-(p-1)
  for j in range(start,len(f)):f[j]-=c
 return f+[0]*(p-1-len(f))
def norm(f,p):
 cols=[reduce_phi([0]*k+list(f),p) for k in range(p-1)]
 return int(fmpz_mat([[cols[j][i] for j in range(p-1)] for i in range(p-1)]).det())
def interp(v,p):
 return [v[0]]+[-sum(v[b]*pow(b,p-1-k,p) for b in range(p))%p for k in range(1,p)]
rng=random.Random(20261108);fail=[];f3n=f4n=f2n=0
for p in [5,7,11]:
 for e in range(1,p):
  for trial in range(4):
   z=[0]*(p-1) if trial==0 else [rng.randrange(-3,4) for _ in range(p-1)]
   raw=[0]*(p-1+e)
   for i in range(e+1):
    for j,c in enumerate(z):raw[i+j]+=(-1)**i*math.comb(e,i)*c
   rho=reduce_phi(raw,p)+[0];k=(p-sum(rho))//p;rho=[c+k for c in rho]
   H=epoly(raw,p,p);lhs=epoly(rho,p,p);wanted=H[:];wanted[p-1]=(wanted[p-1]+1)%p
   gamma=interp([c%p for c in rho],p)
   ok=lhs==wanted and gamma[p-1]==0 and gamma[0]==(1+H[p-1])%p
   ok=ok and all(gamma[p-1-n]==-math.factorial(n)*H[n]%p for n in range(1,p-1))
   if not ok:fail.append(('F3',p,e,trial))
   f3n+=1
 for trial in range(16):
  v=[[rng.randrange(-3,4) for j in range(p)] for i in range(p)]
  v[0][0]-=sum(map(sum,v))%p;s=sum(map(sum,v))
  pi=[[sum(v[i][j] for i in range(p) for j in range(p) if (i+b*j)%p==t) for t in range(p)] for b in range(p)]
  inf=[sum(v[i][j] for i in range(p)) for j in range(p)]
  recovered=[[(inf[j]+sum(pi[b][(i+b*j)%p] for b in range(p))-s)//p for j in range(p)] for i in range(p)]
  if recovered!=v:fail.append(('F2',p,trial))
  f2n+=1
  for compatible in [True,False]:
   pp=[r[:] for r in pi];ii=inf[:]
   if not compatible:pp[0][0]+=1;pp[0][1]-=1
   tau=(s//p)%p;H=[epoly(r,p,p) for r in pp];Hi=epoly(ii,p,p)
   for hh in H+[Hi]:hh[p-1]=(hh[p-1]-tau)%p
   hp=[interp([H[b][n] for b in range(p)],p) for n in range(p)]
   moments=all(not any(hp[n][n+1:]) and Hi[n]==hp[n][n] for n in range(1,p-1))
   moments=moments and (sum(H[b][p-1] for b in range(p))+Hi[p-1]+tau)%p==0
   cells=all((ii[j]+sum(pp[b][(i+b*j)%p] for b in range(p))-s)%p==0 for i in range(p) for j in range(p))
   if moments!=cells or cells!=compatible:fail.append(('F4',p,trial,compatible))
   f4n+=1
print('foundation_controls',{'F3':f3n,'F2':f2n,'F4':f4n,'failures':fail})
p=11;Y=[-2,-2,-3,-2,-1,-2,-1,-3]
logs=logseries(epoly(Y,p,p-1),p,p-1)
print('GEN1',{'norm':norm(Y,p),'Y1':sum(Y),'odd_coords':[logs[j] for j in [3,5,7,9]],'conjugacy_t3':[6*pow(3,7,11)%11,5*pow(3,9,11)%11]})
src=(PACK/'numerical_statement_addendum.txt').read_text();v=json.loads(src[src.index('\n[\n'):])
A=fmpz_mat([[v[(i-k)%p][(j-l)%p] for k in range(p) for l in range(p)]for i in range(p) for j in range(p)])
print('DET1',{'augmentation':sum(map(sum,v)),'max_abs':max(map(abs,sum(v,[]))),'determinant':int(A.det()),'expected':11**15*77419})
for p,r,I in [(37,1,[16]),(59,1,[22]),(67,1,[29]),(103,1,[12]),(103,2,[12]),(103,3,[12]),(103,4,[12]),(131,1,[11]),(157,1,[31,55])]:
 V=[i//r for i in I if i%r==0];w=max(V)+1;P=[1,1]+[0]*(w-2)
 for nu in V:P[nu]=(-logseries(P,p,w)[nu])%p
 H=[[0]*p for b in range(p)];Hi=[0]*p;Hi[2*r+1]=1
 for nu,c in enumerate(P):
  n=1+2*r*nu
  if c and n<p:
   for b in range(p):H[b][n]=c*pow(b,(2*r+1)*nu,p)%p
 failed=[]
 for n in range(1,p-1):
  pol=interp([H[b][n] for b in range(p)],p)
  degree=max([j for j,c in enumerate(pol) if c]+[-1])
  if degree>n or Hi[n]!=pol[n]:failed.append((n,degree))
 if p in [37,59,67]:print('IR4',p,[(n,c) for n,c in enumerate(P) if c],failed)
 else:
  def projection(h):
   return [(h[p-1]-sum(math.factorial(n)*h[n]*pow(t,p-1-n,p) for n in range(1,p-1)))%p for t in range(p)]
  pi=[projection(h) for h in H];ii=projection(Hi)
  wrong=[]
  for i in range(p):
   for j in range(p):
    if (ii[j]+sum(pi[b][(i+b*j)%p] for b in range(p)))%p:wrong.append((i,j))
  print('OBS1',p,r,'failed_moments',failed,'failed_cells',len(wrong),'exact_xy_nonzero',set(wrong)=={(i,j) for i in range(1,p) for j in range(1,p)})
sieve=bytearray([1])*20000;sieve[0]=sieve[1]=0
for i in range(2,math.isqrt(19999)+1):
 if sieve[i]:sieve[i*i:20000:i]=bytearray(len(range(i*i,20000,i)))
count=collections.Counter();irregular=0;small=[];bounds=[]
for p in range(7,20000):
 if not sieve[p]:continue
 invfacts=[1];fac=1
 for j in range(1,p):fac=fac*j%p;invfacts.append(pow(fac,-1,p))
 ber=nmod_poly(invfacts[1:],p).inverse_series_trunc(p-2)
 I=[i for i in range(1,(p-3)//2+1) if int(ber[2*i])==0]
 r=next(r for r in range(1,p) if all(i%r or (2*r+1)*i>r*(p-1) for i in I))
 count[r]+=1;irregular+=bool(I)
 if I and p<400:small.append((p,[2*i for i in I],r,3+2*r))
 if r>(p-5)//2 or p+3+2*r>7*(p+1)//4:bounds.append((p,r))
print('CEN1',{'number_primes':sum(count.values()),'irregular':irregular,'frequencies':dict(sorted(count.items())),'r_gt_1':sum(c for r,c in count.items() if r>1),'bound_failures':bounds})
print('CEN1_below400',small)
