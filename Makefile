
INCS	=	

DEFS    =	-Wno-deprecated -DNORANDOM -DPATCH -DSPHERICAL #-Wno-deprecated -DSTRIPED -DSEPARATOR -DGENERATOR -DMEMBRANE -DSQUARE -DSYMMETRY -DVOLUME 

LIBS	=	-lm

LIBPATH =

COMP	=	g++ #-O4 #g++ 

all: wetS # wetP

String.o: String.cc
	$(COMP) -c String.cc $(DEFS)

RT_Timer.o: RT_Timer.cc
	$(COMP) -c RT_Timer.cc $(DEFS)

Timer.o: Timer.cc
	$(COMP) -c Timer.cc $(DEFS)

wetP: wet.cc String.o RT_Timer.o Timer.o
	$(COMP) wet.cc -o wetP String.o RT_Timer.o Timer.o $(INCS) $(DEFS) -DPARALLEL $(LIBPATH) $(LIBS)

wetS: wet.cc String.o RT_Timer.o Timer.o
	$(COMP) wet.cc -o wetS String.o RT_Timer.o Timer.o $(INCS) $(DEFS) $(LIBPATH) $(LIBS)

clean:
	rm -f *.o wetP wetS wetS.exe *~ *.bak

depend:
	makedepend $(INCS) $(DEFS) wet.cc

# DO NOT DELETE

wet.o: /usr/include/math.h /usr/include/features.h /usr/include/sys/cdefs.h
wet.o: /usr/include/gnu/stubs.h /usr/include/bits/huge_val.h
wet.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
wet.o: /usr/include/stdlib.h
wet.o: /usr/lib/gcc-lib/i386-redhat-linux/2.96/include/stddef.h
wet.o: /usr/include/sys/time.h /usr/include/bits/types.h /usr/include/time.h
wet.o: /usr/include/bits/time.h /usr/include/sys/select.h
wet.o: /usr/include/bits/select.h /usr/include/bits/sigset.h RT_Timer.hh
wet.o: Timer.hh Bool.hh String.hh wet.hh
