FLAGS = -O2 -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic
output_dir=.

ABCDIR = $(HOME)/work/AbcSmc
GSL_PATH = $(HOME)/work/AbcSmc/gsl_local
SQLDIR = $(ABCDIR)/sqdb
ABC_LIB = -L$(ABCDIR) -labc -ljsoncpp -lsqdb $(ABCDIR)/sqlite3.o
GSL_LIB = -lm -L$(GSL_PATH)/lib/ -lgsl -lgslcblas -lpthread -ldl

INCLUDE = -I$(ABCDIR) -I$(GSL_PATH)/include/
default: cholera

cholera: main.cpp DiffEq_Sim.h Cholera_Sim.h
	g++ $(FLAGS) $(INCLUDE) main.cpp -o cholera $(GSL_LIB)

current_version: main_EventDrivenIBM.cpp EventDriven_Sim_Teunis_waning.hpp EventDriven_parameters.hpp | $(output_dir)/polio_data
	g++ $(FLAGS) main_EventDrivenIBM.cpp -o polio

debug: main_Gillespie.cpp
	g++ $(FLAGS) main_Gillespie.cpp -o debug

polio: main.cpp
	g++ $(FLAGS) main.cpp -o polio

edma: EventDriven_MassAction_Sim.h mass_action_SIR.cpp
	g++ $(FLAGS) mass_action_SIR.cpp -o edma

pde: PDE_Simulator.cpp
	g++ $(FLAGS) PDE_Simulator.cpp -o pde

pde_abc: libabc PDE_Simulator_abc.cpp
	g++ $(FLAGS) $(INCLUDE) -I$(SQLDIR) PDE_Simulator_abc.cpp -o pde_abc $(ABC_LIB) $(GSL_LIB)

libabc:
	$(MAKE) -C $(ABCDIR) -f Makefile

polio_boxcar: Polio_boxcar_simulator.cpp Polio_boxcar_model_extended_waning.h DiffEq_Sim.h Polio_boxcar_parameters.h
	g++ $(FLAGS) $(INCLUDE) Polio_boxcar_simulator.cpp -o polio_boxcar $(GSL_LIB)

polio_boxcar_cel: Polio_boxcar_simulator.cpp Polio_boxcar_model_extended_waning.h DiffEq_Sim.h Polio_boxcar_parameters.h
	g++ $(FLAGS) Polio_boxcar_simulator.cpp -o polio_boxcar_cel -lgsl

$(output_dir)/polio_data:
	mkdir $@
