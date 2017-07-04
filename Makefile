FLAGS=-O2 --std=c++11

default: edma

debug: main_Gillespie.cpp
	g++ $(FLAGS) main_Gillespie.cpp -o debug

polio: main.cpp
	g++ $(FLAGS) main.cpp -o polio

edma: EventDriven_MassAction_Sim.h mass_action_SIR.cpp
	g++ $(FLAGS) mass_action_SIR.cpp -o edma
