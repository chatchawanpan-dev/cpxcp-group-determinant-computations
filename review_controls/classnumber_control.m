SetSeed(1);
SetColumns(1024);
print "Magma version:",GetVersion();
for p in [11,13] do
 K:=CyclotomicField(p);
 print p,ClassNumber(K);
end for;
quit;
