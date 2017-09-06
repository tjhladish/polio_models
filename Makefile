FLAGS=-O2 --std=c++11
output_dir=.

default: current_version edma polio

current_version: EventDrivenIBM.cpp EventDriven_Sim_nonexp_recovery_waningfn_ENV.hpp EventDriven_parameters.hpp | $(output_dir)/polio_data
	g++ $(FLAGS) EventDrivenIBM.cpp -o polio

debug: main_Gillespie.cpp
	g++ $(FLAGS) main_Gillespie.cpp -o debug

polio: main.cpp
	g++ $(FLAGS) main.cpp -o polio

edma: EventDriven_MassAction_Sim.h mass_action_SIR.cpp
	g++ $(FLAGS) mass_action_SIR.cpp -o edma

$(output_dir)/polio_data:
	mkdir $@
