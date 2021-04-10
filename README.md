# NeuroSim-TM
Simulator of spiking dynamics in networks of Leaky Integrate-and-fire (LIF) neurons with Tsodyks-Uziel-Markram short-term synaptic plasticity

## A Brief User Guide for NeuroSim-TM-2.1 Simulator

### 1. Introductory Overview

The custom-made software `NeuroSim-TM-2.1` written in the C programming language is used for simulation of spiking dynamics of generative networks composed of Leaky Integrate-and-Fire (LIF) point-like neurons with synaptic plasticity. The software has been specifically designed for numerical studies of spatiotemporal properties of population spikes (or population bursts) [1].

Synaptic plasticity may be only the short-term one, according to the Tsodyks-Uziel-Markram (TUM) model [2], or it can be complemented by ‘long-term’ Spike-Timing-Dependent Plasticity (STDP). Ordinary differential equations for the dynamics of neuron potentials and synaptic variables are solved using the standard Euler method. The default initial conditions are that all neurons have the same ’resting’ value of the membrane potential, and all synapses have the same initial values of synaptic variables.

Both metric (i.e. spatially dependent) and metric-free network connectomes can be generated. For the metric networks, the software allows various spatial arrangements of neurons. In this case, the probability of forming a unilateral connection between every two neurons decreases exponentially with increasing their relative distance. The maximum size of the generated network is limited by the amount of operative memory available to a single-thread program, as multi-thread parallelization is not yet implemented.

After compilation with a standard C compiler (e.g., GCC), the executable file for `NeuroSim-TM-2.1` can be used for all available variants of simulation without recompiling. All parameters required for choosing and performing simulations are set through the single input file input.txt, which should be placed in the same folder as the executable file.

The software also supports several protocols for specific numerical experiments on influence of various forms of neuron stimulation on network activity, as well as a simple model for global homeostasis. These and other capabilities unused in the reported study will be omitted here.

Based on the data in input.txt (see `Table 1` below), the program simulates neuronal network dynamics during a specified time and writes the resulting data in a set of output files. All simulation parameters of individual neurons and synapses are also exported in the output files allowing an accurate reproducibility of the simulation, except for the optional mode of stochastic spontaneous spike generation. In total, there are 18 output files (including 4 files in subfolder `\syn_resources_dynamics`). There are three sorts of these files: primary output files of network dynamics, files of parameters of neuronal network (connectome, etc.) that are required for reproducing the simulation, and auxiliary files made for fast analysis of the results and other purposes.



The primary output files are as follows:
-	`activity.txt` contains normalized spiking activity of the network neurons as a function of chunked time,
-	`raster.txt` contains all spike generation times of individual neurons,
-	`x.txt`, `y.txt`, `z.txt`, `u.tx`t (in subfolder `\syn_resources_dynamics`) contain the dynamics of network-averaged fractions of synaptic resources in the TUM model [2].

The files of neuronal network parameters are as follows:

-	`coordinates.txt` contains spatial Cartesian coordinates of neurons on two-dimensional plane,
-	`connections.txt` contains the network connectome (a directed and unweighted graph), i.e. complete information on all connections between neurons,
-	`synaptic_parameters_distribution.txt` contains the parameters of every synapse of the network,
-	`exc_I_distribution.txt` and `inh_I_distribution.txt` contain the values of constant ‘background’ current (which determines neuronal excitability and the fraction of steady pacemakers among the neurons) for every excitatory and inhibitory neuron, respectively,
-	or, alternatively to the last, `p_sp_distribution.txt` contains the values of the probability of spontaneous spike generation for every neuron.

Finally, the auxiliary files are as follows: `active_connections.txt`, `M.txt`, `w_initital.txt`,
-	`connections_per_neuron.txt` contains the number of outgoing synapses for a neuron with the given ID (the neuron IDs range from zero to N-1, where N is the total number of neurons), i.e. there are two integer numbers per row.
-	`burst_times.txt` contains the population spike (or population burst, PB) times, if the PB amplitude exceeds the threshold value burst_threshold stated in input.txt,
-	`IBI.txt` contains successive time intervals between the detected PBs,
-	`info.txt` contains the overview of the mode and key parameters for the particular simulation run.

To interactively control the network activity during the simulation, Gnuplot scripts plot_activity.plt and, for small networks (see Figure 1 in the main article), plot_raster.plt have been used. These scripts visualize the current content of the output files activity.txt and raster.txt, respectively, which are updated regularly during the simulation.

[1]: A.V. Paraskevov, D.K. Zendrikov, A spatially resolved network spike in model neuronal cultures reveals nucleation centers, circular traveling waves and drifting spiral waves, Phys. Biol. 14 (2017) 026003.
[2]: M. Tsodyks, A. Uziel, H. Markram, Synchrony generation in recurrent networks with frequency-dependent synapses, J. Neurosci. 20 (2000) RC50.


### 2. Compiling the simulator

The software can be standardly compiled with the GCC compiler in both Windows and Linux operating systems using console command:
```
gcc NeuroSim-TM-2.1.c -o NeuroSim-TM-2.1
```
### 3. Starting the simulation

A copy of the executable should be placed in the directory of a particular run. There it will search for the input files with specific names and generate its output.
The input settings are specified in the file `input.txt`. The syntax of the input file is `PARAMETER_NAME=VALUE` at the beginning of a new line without spaces. All other symbols are discarded as comments. The parameters can be specified in any order. Please find the list of parameters in Table 1.
If some of the parameters are not found in `input.txt`, the default values from `NeuroSim-TM-2.1.c` (please see lines 483-571 there) are used. Before the simulation beginning the user is prompted to check if the values imported from the input file are correct and to type ‘Y’ (yes) for confirmation.

After that, the static parameters of the neuronal network are generated (or loaded from the previously saved files; see details below) and written in the output files:
-	`coordinates.txt` with coordinates X and Y of one neuron per row, in units of L,
-	`connections.txt` with the indices of pre- and post-synaptic neurons, and the generated random value from the uniform distribution that should be greater than the connection probability for this pair of neurons (i.e, there are three numbers, two integers and one real, per row),
-	`synaptic_parameters_distribution.txt` with pre- and post-synaptic neuron indices, and the values of synaptic impulse amplitude and time constants A (A>0 if presynaptic neuron is excitatory one, and A<0 if it is inhibitory one), U, tau_I, tau_rec, tau_facil for one synapse per row,
-	`exc_I_distribution.txt` with an excitatory neuron index and the value of its constant background current per row
-	`inh_I_distribution.txt` with an inhibitory neuron index and the value of its constant background current per row (if the neuronal network does not contain inhibitory neurons, the file is empty),
-	`p_sp_distribution.txt` with a neuron index and the value of probability of spontaneous spike generation (per time step dt) by the neuron per row.

These files are then closed at the beginning of the simulation (indicated by the console message ‘Network is created successfully’) and can be used for the analysis of neuronal network realization, for performing other runs with this realization, and for reproducing the simulations.

During the simulation, the program prints in its console window the spike times and indices of spiking neurons, and the number of population spikes detected by crossing the population activity threshold burst_threshold stated in input.txt. At the same time, the dynamical data are periodically written to the output files and can be accessed from there (e.g., by Gnuplot scripts for plotting graphs of network activity and raster of spikes), even if the simulation is still running.

After the simulation finishes, its results are exported to the above-mentioned output files in the simulation folder:
-	`activity.txt` containing (per row) the timestamp and normalized network activity averaged over `AVG_TIME` time chunk stated in `input.txt`,
-	`raster.txt` containing (per row) spiking neuron index, i.e. neuron ID, and the corresponding spike time,
-	`x.txt`, `y.txt`, `z.txt`, `u.txt` in subfolder `\syn_resources_dynamics` contain (per row) the timestamp and the corresponding fractions x, y, z, u of synaptic resources in the TUM model of short-term synaptic plasticity [2], averaged over all synapses of the network and over the time chunk `AVG_TIME`.
Note: As the output data are event-driven (the events are neuronal spikes), the simulation may run without any messages about spikes in the console window. That happens (e.g., during some time after the beginning of simulation) if all neurons are in sub-threshold regime, i.e. for some time they do not receive enough input excitation to generate a spike. In this case, the simulation is still running, and the program window should not be closed until the message `'Releasing memory... Done'` indicating the end of the simulation.


### 4. Reproduction of simulation results

To use previously generated components of the neuronal network realization, please carry out the following sequence:

1. Copy all existing files of neuronal network parameters described in the previous section, and input.txt, from the folder with the accomplished simulation to the folder of a new ‘reproducing’ simulation, where the executable file NeuroSim-TM-2.1 is placed.

2. Append prefix saved_ to the corresponding files described in the previous section. In other words, rename coordinates.txt into saved_coordinates.txt, connections.txt into saved_connections.txt, synaptic_parameters_distribution.txt into saved_synaptic_parameters_distribution.txt etc. However, the file input.txt should not be renamed.

3. In the file input.txt, change values of the flags listed below from 0 to 1:

    3.1. Set use_saved_coordinates=1 to use spatial coordinates of neurons from the file saved_coordinates.txt for the new simulation, instead of generating new coordinates in the case if use_saved_coordinates=0.

    3.2. Set use_saved_topology=1 to use network connectome graph from the file saved_connections.txt. If spatial neuronal network is simulated, then one must previously set use_saved_coordinates=1 for consistency (there will be error message otherwise). If use_saved_topology=0, new network connectome graph will be generated based on the loaded or generated coordinates of neurons.

    3.3. Set use_saved_synaptic_parameters=1 to use the synaptic parameters from the file saved_synaptic_parameters_distribution.txt. Again, one must previously set use_saved_topology=1 for consistency. If use_saved_synaptic_parameters=0, a new set of synaptic parameters will be generated for the loaded or newly-generated network connectome graph.

    3.4. Set use_saved_stimulation_data=1 to use either the values of constant background currents of neurons loaded from the pair of files saved_exc_I_distribution.txt and saved_inh_I_distribution.txt or the corresponding values of the probability of spontaneous spike generation loaded from the file saved_p_sp_distribution.txt. The number of neurons must be consistent with that stated in the input.txt and implied in the data loaded from all relevant files. If use_saved_stimulation_data=0, the corresponding values will be generated anew.

Notes: 1) Due to deterministic nature of the neuron and synapse models, simulations can be reproduced with exact spike times given the same input files (except for the modes of spontaneous neuron spiking). 2) If any inconsistencies between the input files are detected, such as mismatch in overall number of neurons, excitatory/inhibitory types of neurons etc., the simulation will not be started and the corresponding error message will be shown.
 

