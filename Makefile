vpath %.f src
vpath %.F src
vpath %.f90 src
vpath %.o OBJ
FC=ifort 
FC1=gfortran -ffixed-line-length-none

OBJS = MainProgram.o sacio.o green.o bessel.o comvar.o dwim.o fft.o Gbasic.o Ginput.o Grtcoefs.o G-s-r.o\
	 Gsubs.o integrang_calc.o MTX_SUB.o iptam.o SOURCE.o UKO.o Ydumtx.o

OBJ_PATH = OBJ
SRC_PATH = src
BIN_PATH = bin
#gfortran=gfortran
phony_list := all grtm trav
.PHONY: ${phony_list}

all: grtm trav
grtm:$(BIN_PATH)/grtm
trav:$(BIN_PATH)/trav

$(BIN_PATH)/grtm : $(OBJS)
	${FC} $^ -o $@

$(BIN_PATH)/trav : $(OBJ_PATH)/tau_p.o $(OBJ_PATH)/trav.o
	${FC1} $^ -o $@

$(OBJ_PATH)/MainProgram.o : src/MainProgram.f
	${FC} -c $< -o $@

$(OBJ_PATH)/sacio.o : src/sacio.c
	cc -c $< -o $@

$(OBJ_PATH)/green.o : src/green.f90
	${FC} -c $< -o $@

$(OBJ_PATH)/bessel.o : src/bessel.f
	${FC} -c $< -o $@

$(OBJ_PATH)/comvar.o : src/comvar.f90
	${FC} -c $< -o $@

$(OBJ_PATH)/dwim.o : src/dwim.f90
	${FC} -c $< -o $@

$(OBJ_PATH)/fft.o : src/fft.f
	${FC} -c $< -o $@

$(OBJ_PATH)/Gbasic.o : src/Gbasic.F
	${FC} -c $< -o $@

$(OBJ_PATH)/Ginput.o : src/Ginput.f90
	${FC} -c $< -o $@

$(OBJ_PATH)/Grtcoefs.o : src/Grtcoefs.F
	${FC} -c $< -o $@

$(OBJ_PATH)/G-s-r.o : src/G-s-r.f
	${FC} -c $< -o $@

$(OBJ_PATH)/Gsubs.o : src/Gsubs.f
	${FC} -c $< -o $@

$(OBJ_PATH)/integrang_calc.o : src/integrang_calc.f90
	${FC} -c $< -o $@

$(OBJ_PATH)/MTX_SUB.o : src/MTX_SUB.F
	${FC} -c $< -o $@

$(OBJ_PATH)/iptam.o : src/iptam.f90
	${FC} -c $< -o $@

$(OBJ_PATH)/SOURCE.o : src/SOURCE.F
	${FC} -c $< -o $@

$(OBJ_PATH)/UKO.o : src/UKO.F
	${FC} -c $< -o $@

$(OBJ_PATH)/Ydumtx.o  : src/Ydumtx.f
	${FC} -c $< -o $@

$(OBJ_PATH)/tau_p.o :src/tau_p.f
	${FC1} -c $< -o $@

$(OBJ_PATH)/trav.o : src/trav.f
	${FC1} -c $< -o $@

clean:
	rm OBJ/*.o bin/grtm bin/trav
