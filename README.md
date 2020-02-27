# LEF-Simulation

# SIMULATION CODE FOR THE SLIDING MODEL OF COMPOSITE LEF SIMULATION

# Main.m: Main program of loop extrusion simulation. The 1st part is parameter assginment. The 2nd part is parameter initiation. The 3rd part is event execution. The 4th part is contact map generation. The LEF model has the following features: 1 one-sided extrusion. 2 LEF loading at gene start sites. 3. LEF cannot go across each other.

# do_event.m: Event execution file. Four events are programmed: 1 LEF translocation. 2 LEF initiation at gene start sites. 3 LEF unbinding. 4 LEF loading to gene start sites.

# regenerate_event.m: Event regeneration. Every time a new event is considered, this program calculates the event time using exponential distribution and puts the event into the heap structure.

# regenerate_neighbours.m: Every time an event is executed, it checks whether there are neighbouring LEFs whose events need updating.

# LEFSystem.m: Class object that saves all the variables throughout the simulation. All the variables are saved as properties in the object. In addition, it sets up the permeability of all the positions on the chromatin and initial positions of the LEFs. It also provides basic functions for LEF actions.

# Event_heap.m: Heap structure. Events are saved as entries in this object. Each event has two properties: event id and time.

# Event_t.m: Event property assignment.

# heatmap_3d.m: Contact map generation. The inputs include the time series of loop configurations and the number of maps averaged. The output is the averaged contact map.

# N_eff_TwoPoints.m: Effective genomic distance calculation. Given two points and the loop configuration, it returns the effective genomic distance by using partial Gaussian integral with the inversed covariance matrix.

# FindPairs.m: Given two points between which effective genomic distance is calculated, this program finds all the LEF pairs that need to be invloved in the calculation.

# Cal_InvCovM.m: Given all the nodes in the loop configuration, this program returns the inverse of the covariance matrix.

# PolymerSim.m: The evolution of loop configurations is taken as the input. Using Rouse model chain, the polymers with and without loops are simulated. MSDs of different beads are then calculated.
