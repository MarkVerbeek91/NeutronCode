IDIR =../include
CC=gcc
CFLAGS=-I$(IDIR)

ODIR=obj
LDIR =../lib

LIBS=-lm

_DEPS = includes.h Auxiliary.h CathodeCurrents.h CrossSections.h CrossSectionsCont.h EnergySpectrumIons.h IonParticleFlux.h MathFunctions.h NeutralParticleFlux.h NeutronProductionRate.h PotentialFunctions.h SourceFunction.h SurvivalFunctions.h constants.h

DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o Auxiliary.o CathodeCurrents.o CrossSections.o EnergySpectrumIons.o IonParticleFlux.o MathFunctions.o NeutralParticleFlux.o NeutronProductionRate.o PotentialFunctions.o SourceFunction.o SurvivalFunctions.o constants.o

OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.o $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

NeutronCode: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

main.o:

Auxiliary.o:

CathodeCurrents.o:

CrossSections.o:

EnergySpectrumIons.o:

IonParticleFlux.o:

MathFunctions.o:

NeutralParticleFlux.o:

NeutronProductionRate.o:

PotentialFunctions.o:

SourceFunction.o:

SurvivalFunctions.o:

constants.o:
	
.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 