FLAGS=-O2 --std=c++11
output_dir=.

default: current_version 

current_version: main_EventDrivenIBM.cpp EventDriven_Sim_Teunis_waning.hpp EventDriven_parameters.hpp | $(output_dir)/polio_data
	g++ $(FLAGS) main_EventDrivenIBM.cpp -o polio

debug: main_Gillespie.cpp
	g++ $(FLAGS) main_Gillespie.cpp -o debug

polio: main.cpp
	g++ $(FLAGS) main.cpp -o polio

edma: EventDriven_MassAction_Sim.h mass_action_SIR.cpp
	g++ $(FLAGS) mass_action_SIR.cpp -o edma

$(output_dir)/polio_data:
	mkdir $@
