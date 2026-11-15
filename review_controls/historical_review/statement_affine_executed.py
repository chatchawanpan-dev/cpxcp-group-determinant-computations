import random, math, json
def mul(a,b,p,w):
 c=[0]*w
 for i,x in enumerate(a):
  for j,y in enumerate(b):
   if i+j<w:c[i+j]=(c[i+j]+x*y)%p
 return c
def logseries(a,p,w):
 assert a[0]%p
 v=[x*pow(a[0],-1,p)%p for x in a]+[0]*w;v=v[:w];v[0]=0
 ans=[0]*w;power=[1]+[0]*(w-1)
 for k in range(1,w):
  power=mul(power,v,p,w)
  scale=pow(k,-1,p)*(1 if k%2 else -1)
  ans=[(x+scale*y)%p for x,y in zip(ans,power)]
 return ans
def op(q,d,x,p,w):
 l=[0]*w;l[1]=q
 for j in range(3,w,2):l[j]=d.get(j,0)
 psi=[1]+[0]*(w-1)
 for n in range(1,w):psi[n]=sum(j*l[j]*psi[n-j] for j in range(1,n+1))*pow(n,-1,p)%p
 eps=[0]*w;eps[0]=1
 for i in range(2,w,2):eps[i]=(x.get(i,0)-sum(eps[j]*psi[i-j] for j in range(0,i,2)))%p
 return mul(psi,eps,p,w)
def interp(v,p):
 return [v[0]]+[-sum(v[b]*pow(b,p-1-k,p) for b in range(p))%p for k in range(1,p)]
def evalpoly(c,b,p):return sum(x*pow(b,k,p) for k,x in enumerate(c))%p
def residuals(p,e,D,dinf,tau,a,coeffs):
 M=(p-3)//2
 ds=interp([z.get(3,0) for z in D],p)
 eta=((-ds[5]*pow(a,-2,p)) if e==2 else ds[4]*pow(a,-1,p))%p
 h=[[0]*p for _ in range(p)]
 for b in range(p):
  x={2*i:evalpoly(coeffs.get(2*i,[]),b,p) for i in range(1,M+1)}
  A=op(a*b*b%p if e==2 else 0,D[b],x,p,p-1)
  for n in range(1,p):h[n][b]=A[n-1]
 hp=[interp(v,p) for v in h];t=[hp[n][n] for n in range(p)]
 out=[hp[n][k] for n in range(4,p-2,2) for k in range(n+1,p)]
 if e==2:
  out += [(t[2]-a)%p,(t[3]-a*eta)%p]
  Q=op(eta,dinf,{2*i:t[2*i+2]*pow(a,-1,p)%p for i in range(1,M)},p,p-3)
  out += [(t[j+2]-a*Q[j])%p for j in range(3,p-3,2)]
 else:
  out += [(t[3]-a)%p,(t[4]-a*eta)%p]
  Q=op(eta,dinf,{2*i:t[2*i+3]*pow(a,-1,p)%p for i in range(1,M)},p,p-3)
  out += [(t[j+3]-a*Q[j])%p for j in range(3,p-5,2)]
  out += [(tau+sum(h[p-1])+a*Q[p-4])%p]
 return out
def consistent(p,e,D,dinf,tau):
 M=(p-3)//2;ds=interp([z.get(3,0) for z in D],p)
 degree=max([k for k,x in enumerate(ds) if x]+[-1])
 if e==2:
  if degree!=6:return False
  aa=[a for a in range(1,p) if (a**3-3*ds[6])%p==0]
 else:
  if degree>4:return False
  aa=range(1,p)
 basis=[(2*i,k) for i in range(1,M+1) for k in range(2*i+2)]
 for a in aa:
  c=residuals(p,e,D,dinf,tau,a,{})
  cols=[]
  for i,k in basis:
   z=[0]*(k+1);z[k]=1
   col=residuals(p,e,D,dinf,tau,a,{i:z})
   cols.append([(x-y)%p for x,y in zip(col,c)])
  A=[[cols[j][i] for j in range(len(cols))]+[(-c[i])%p] for i in range(len(c))]
  row=0
  for col in range(len(cols)):
   pivot=next((i for i in range(row,len(A)) if A[i][col]),None)
   if pivot is None:continue
   A[row],A[pivot]=A[pivot],A[row];scale=pow(A[row][col],-1,p)
   A[row]=[x*scale%p for x in A[row]]
   for i in range(len(A)):
    if i!=row and A[i][col]:
     scale=A[i][col];A[i]=[(x-scale*y)%p for x,y in zip(A[i],A[row])]
   row+=1
  if all(any(r[:-1]) or not r[-1] for r in A):return True
 return False
rng=random.Random(20261108);tests=[];bad=[]
for p in [7,11,13]:
 for e in [2,3]:
  for trial in range(8):
   a=rng.randrange(1,p);polys={1:[1],2:[0,0,a] if e==2 else [0]}
   for n in range(3,p):polys[n]=[rng.randrange(p) for _ in range(n+1)]
   if e==3:polys[3][3]=a
   hv={n:[evalpoly(polys[n],b,p) for b in range(p)] for n in range(1,p)}
   D=[{j:logseries([hv[n][b] for n in range(1,p)],p,p-1)[j] for j in range(3,p-1,2)} for b in range(p)]
   hinf=[0]*p
   for n in range(1,p-1):hinf[n]=polys[n][n] if len(polys[n])>n else 0
   tau=rng.randrange(p);hinf[p-1]=(-tau-sum(hv[p-1]))%p
   norm=[hinf[e+k]*pow(a,-1,p)%p for k in range(p-e)]
   ld=logseries(norm,p,p-e)
   dinf={j:ld[j] for j in range(3,p-e,2)}
   coeffs={2*i:polys[2*i+1] for i in range(1,(p-3)//2+1)}
   res=residuals(p,e,D,dinf,tau,a,coeffs)
   if any(res):bad.append(('constructed',p,e,trial,res))
   tests.append((p,e))
 print('A2_counts',p,(p*p-9)//4, (p-3)**2//4+1)
print('affine_constructed_moment_packets',len(tests),'failures',bad)
results=[]
for p in [11,13]:
 js=list(range(3,p-1,2));coords=[tuple([0]*len(js))]
 for k in range(len(js)):
  for x in [1,p-1]:coords.append(tuple(x if i==k else 0 for i in range(len(js))))
 for t in range(8):coords.append(tuple([0]*(len(js)-2)+[rng.randrange(p),rng.randrange(p)]))
 for t in range(5):coords.append(tuple(rng.randrange(p) for _ in js))
 for c in coords:
  expected=not any(c[:-2]) and any(c[-2:]);admitted=False
  # DF3 permits all conjugates. Fixed-ideal leading-one tests alone are insufficient.
  for t in range(1,p):
   dc={j:c[k]*pow(t,j,p)%p for k,j in enumerate(js)}
   D=[dc]+[{} for _ in range(p-1)]
   generic=consistent(p,3,D,{},1)
   exceptional=consistent(p,3,[{} for _ in range(p)],dc,1)
   two=consistent(p,2,D,{},1)
   if two:bad.append(('lowest_prime',p,c,t))
   if generic and (dc[js[-2]] or dc[js[-1]]!=p-1):bad.append(('generic_normalization',p,c,t))
   admitted=admitted or generic or exceptional
   if admitted:break
  if admitted!=expected:bad.append(('prime_orbit',p,c,admitted,expected))
  results.append((p,c,admitted))
 print('prime_orbit_controls',p,len(coords))
print('affine_and_prime_control_failures',bad)
print('zero_unit_tau_controls',[(p,consistent(p,3,[{} for _ in range(p)],{},1),consistent(p,3,[{} for _ in range(p)],{},0)) for p in [7,11,13]])
