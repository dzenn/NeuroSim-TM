#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define MIN(a,b) (((a)>(b))?(b):(a))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define SGN(x) ((x < 0) ? -1 : 1)

/* This program simulates spiking activity of N neurons and creates
 * two main *.txt output files with IDs of neurons that emit spikes and their spiking times (Raster data)
 * and spiking activity.
 *
 * All neuron information is stored in an array of structures for each neuron.
 * Each structure contains parameters of a neuron,
 * a list of output synapses and a list of pointers on input synapses.
 *
 * Connecntions between neurons are created randomly with probability P_CON for BINOMIAL topology,
 * or with exponentially decreasing probability in case of spatially-dependent topology
 *
 * Structure of the program:
 * 1. Reading input data
 * 2. Initializing an array of neurons, setting all initial values, allocating memory
 * 3. Linking neurons, creating synapses, writing (and reading, if set so) all supplementary data to files (e.g. connections, background current distributions etc)
 * 4. Simulating & writing data to files
 * 5. Releasing memory
 *
 * Each simulation step goes with two passes through the neuron array, one for decreasing synaptic resources, another for
 * membrane potential step and spike detection
 * */


#define PARAMETER_N 79   //number of parameters expected to be in the input file
#define BUFFER_SIZE 200
#define TYPE_EXCITATORY 1
#define TYPE_INHIBITORY 0

#define STIM_TYPE_I_BG_GAUSSIAN 1
#define STIM_TYPE_I_BG_TWO_VALUES 2
#define STIM_TYPE_P_SP_GAUSSIAN 3
#define STIM_TYPE_P_SP_TWO_VALUES 4

#define GENERATE_TOPOLOGY 0
#define READ_TOPOLOGY_FROM_FILE 1
#define GENERATE_STIMULATION_DATA 0
#define READ_STIMULATION_DATA_FROM_FILE 1
#define READ_SYNAPTIC_DATA_FROM_FILE 1
#define GENERATE_SYNAPTIC_DATA 0
#define READ_COORDINATES_FROM_FILE 1
#define GENERATE_COORDINATES 0
#define STDP_IS_OFF 0
#define STDP_IS_MULTIPLICATIVE 1
#define STDP_IS_ADDITIVE 2
#define HOMEOSTASIS_IS_OFF 0
#define HOMEOSTASIS_IS_ON 1
#define MAXIMUM_NUMBER_OF_NEURONS 1000000

#define LAYER_TYPE_BINOMIAL 0
#define LAYER_TYPE_UNIFORM 1
#define LAYER_TYPE_SQUARE_LATTICE 2
#define LAYER_TYPE_STRIPED 3
#define LAYER_TYPE_BELL 4
#define LAYER_TYPE_RAMP 5
#define LAYER_TYPE_DOUBLE_RAMP 6
#define LAYER_TYPE_BARBELL 7

#define NEURON_MODEL_PERFECT_INTEGRATE_AND_FIRE 0
#define NEURON_MODEL_LEAKY_INTEGRATE_AND_FIRE 1

//****Incoming spikes
struct delay_timer{

    double time_left;

    struct delay_timer *next;

};

//****Each neuron stores a list of output synapses
struct synapse{

    int pre_id; //id of a presynaptic neuron
    int post_id; //id of a postsynaptic neuron
    double I; //synaptic current

    double w; //synaptic weight

    double l; //length of a connection
    double tau_delay; //axonal delay

    double A; //synaptic current magnitude
    double U; //supplementary synaptic magnitude
    double u; //use of synaptic resources
    double x; //recovered
    double y; //active           states of synaptic resources
    double z; //inactive

    double tau_I; //decay constant
    double tau_rec; //syn depression recovery time
    double tau_facil; //u time constant
    struct synapse *next; //pointer on the next synapse in list
    struct delay_timer *timers;

};


//****Each neuron stores a list of input synapses using pointers on output synapses of other neurons
struct inc_synapse{ // item of a list of input synapses of a neuron

    struct synapse *syn_pointer;

    struct inc_synapse *next;

};

//****Structure of one memory unit that describes a neuron
struct neuron{

    int type; //inhibitory or excitatory

    int id;

    double x; //coordinates of a neuron, set to (0,0) for BINOMIAL topology, otherwise described by LAYER_TYPE
    double y;

    int num_of_outcoming_connections; //amount of outcoming connections of a neuron, stored to be easily accessed without re-calculating

    int num_spikes;

    double V; //mV, membrane potential of a neuron
    double V_rest; //mV, resting value of V
    double V_reset; //mV, resetting value of V, V is set to this value after emitting a spike
    double I_b; //pA, background current value
    double I_b_init; //pA, initital value of background current, as the I_b itself can be decreased by the homeostasis mechanism, if the latter is turned on
    double p_sp; //probability of emitting a spike on each step of simulation

    double ref_time_left; //0 until neuron emits spike, then 3 or 2 ms of excitatory and inhibitory neurons respectively

    double last_spiked_at; //time of the last emitted spike, used for STDP rules

    struct synapse *out_conn; //list for output synapses
    struct inc_synapse *in_conn; //list for input synapses

    struct neuron *next;

};

//****Structure used for reading parameters from input file
struct parameter {
    const char *name;
    const void *value;
};


//****Function for reading parameters from input file
void readfile (FILE *file, struct parameter parameter_array[PARAMETER_N]);


//****Function for cutting off unneccessary digits after a point to make sure program uses same values that is written to output files
double my_round (double input){

    char str[100];
    double output;

    str[0]='\0';
    sprintf(str,"%g",input);
    sscanf(str,"%lf",&output);

    return output;

}

//****Function for generating a gaussian distribution
double gauss (double mean, double sd, double min, double max){

    int i; //number of iterations in procedure

    double s;

    while (1){
        s = 0;
        for (i=0; i< 12; i++){
            s += ((double) rand() / (RAND_MAX));
        }
        s = s - 6;
        if (((mean + s*sd) > min)&&((mean + s*sd) < max))
            return my_round(mean + s*sd);
    }

}

double my_theta (double X){
    if (X > 0) return 1;
    else return 0;
}

char *my_itoa(int num, char *str)
{
    if(str == NULL)
    {
            return NULL;
    }
    sprintf(str, "%d", num);
    return str;
}

//****STDP rules, MULTIPLICATIVE or ADDITIVE
double dw(double w, double delta_t, int STDP_status, double A_plus, double A_minus, double tau_plus_corr, double tau_minus_corr){
    double result = 0;
    if(STDP_status == STDP_IS_MULTIPLICATIVE) {
        if(delta_t > 0)
            result = A_plus * (1 - w) * exp(-(delta_t)/tau_plus_corr);
        if(delta_t < 0)
            result = -A_minus * w * exp((delta_t)/tau_minus_corr);
    }
    else if (STDP_status == STDP_IS_ADDITIVE){
        if(delta_t > 0)
            result = A_plus * exp(-(delta_t)/tau_plus_corr);
        if(delta_t < 0)
            result = -A_minus * exp((delta_t)/tau_minus_corr);
    }
    if(w + result < 0.000001)   //boundary correction
        result = -w;
    else if (w + result > 1)
        result = 1 - w;
    return result;
}


//****T_M model synapse step
void stepSynapse(struct synapse *syn, double dt){

    double old_y, old_z, old_u;

    old_y = syn->y;
    old_z = syn->z;
    old_u = syn->u;


    syn->x += dt*old_z/syn->tau_rec;
    syn->y -= dt*old_y/syn->tau_I;
    syn->z += dt*(old_y/syn->tau_I - old_z/syn->tau_rec);
    if (syn->tau_facil == 0) syn->u = syn->w;
    else syn->u -= dt*old_u/syn->tau_facil;

    syn->I = syn->A*syn->y;

    if (syn->x >= 1) syn->x = 0.9999999; if (syn->x <= 0) syn->x = 0.0000001;//correction
    if (syn->y >= 1) syn->y = 0.9999999; if (syn->y <= 0) syn->y = 0.0000001;
    if (syn->z >= 1) syn->z = 0.9999999; if (syn->z <= 0) syn->z = 0.0000001;
    if (syn->u >= 1) syn->u = 0.9999999; if (syn->u <= 0) syn->u = 0.0000001;

}

//****LIF model step
double stepV_LIF(struct neuron *n, double dt, double R_IN, double TAU_M){
    double I_syn = 0;
    struct inc_synapse *tmp_input_synapse = n->in_conn;

    while (tmp_input_synapse->syn_pointer != NULL){ //a sum of all incoming synaptic currents
        I_syn += tmp_input_synapse->syn_pointer->I;
        tmp_input_synapse = tmp_input_synapse->next;
    }

    return (n->V + dt*(-(n->V - n->V_rest) + (I_syn + n->I_b)*R_IN)/TAU_M);

}

//****PIF model step
double stepV_PIF(struct neuron *n, double dt, double C_M){
    double I_syn = 0;
    struct inc_synapse *tmp_input_synapse = n->in_conn;

    while (tmp_input_synapse->syn_pointer != NULL){ //a sum of all incoming synaptic currents
        I_syn += tmp_input_synapse->syn_pointer->I;
        tmp_input_synapse = tmp_input_synapse->next;
    }

    return (n->V + dt*(I_syn + n->I_b)/C_M);

}

//****Homeostasis model step
double stepM(double M,double M_max,double TAU_M,double dt){

    return dt*(-(M - M_max )/TAU_M);

}

//****Function for generating coordinates for every neuron according to LAYER_TYPE
void setCoordinates(struct neuron **neurons, int layer_type, int N){

    double u,v;
    double *x, *y;
    int i,j, random_i, list_i=0, N_rest=N;
    int *tmp_list = malloc(MAXIMUM_NUMBER_OF_NEURONS*sizeof(int));

    for (i = 0; i < N; i++)
        tmp_list[i] = i;


    for (i = 0; i < N; i++){

        x = &neurons[i]->x;
        y = &neurons[i]->y;

        switch (layer_type){

        case LAYER_TYPE_BINOMIAL:
            *x = 0;
            *y = 0;
            break;

        case LAYER_TYPE_UNIFORM:
            *x = my_round((double)(rand())/RAND_MAX);
            *y = my_round((double)(rand())/RAND_MAX);
            break;

        case LAYER_TYPE_SQUARE_LATTICE:
            list_i = (int)floor((double)N_rest*((double)(rand())/RAND_MAX));
            random_i = tmp_list[list_i];
            *y = my_round((1/sqrt((double)N))*floor((double)random_i/(sqrt((double)N))));
            *x = my_round((random_i - *y * (double)N)/sqrt((double)N));
            if (*x < 0.0000001) *x = 0;

            break;

        case LAYER_TYPE_STRIPED:
            v = (double)(rand())/RAND_MAX;
            if (v<0.1) {
                u = 0.4*(double)(rand())/RAND_MAX;
                if (u<0.2) *x = u+0.2;
                else *x = u+0.4;
            } else {
                u = 0.6*(double)(rand())/RAND_MAX;
                if (u < 0.2) *x = u;
                else if (u < 0.4) *x = u+0.2;
                else *x = u+0.4;
            }
            *x = my_round(*x);
            *y = my_round((double)(rand())/RAND_MAX);
            break;

        case LAYER_TYPE_BELL:
            *x = gauss(0.5,0.25,0,1);
            *y = gauss(0.5,0.25,0,1);
            break;

        case LAYER_TYPE_RAMP:
            while(1){
                v = (double)(rand())/RAND_MAX;
                u = (double)(rand())/RAND_MAX;
                if (u <= v) break;
            }
            *x = my_round(v);
            *y = my_round((double)(rand())/RAND_MAX);
            break;

        case LAYER_TYPE_DOUBLE_RAMP:
            while(1){
                v = (double)(rand())/RAND_MAX;
                u = (double)(rand())/RAND_MAX;
                if (u <= v) break;
            }
            *y = my_round((double)(rand())/RAND_MAX);
            if (*y < 0.5) *x = v;
            else *x = 1 - v;

            *x = my_round(*x);
            break;

        case LAYER_TYPE_BARBELL:

            if (i<0.4*N){
                while(1){
                    u = (double)(rand())/RAND_MAX;
                    v = (double)(rand())/RAND_MAX;
                    if (sqrt(pow(u-0.2,2) + pow (v-0.2,2)) < 0.2 && u+v<0.664){
                        *x=u;
                        *y=v;
                        break;
                    }
                }
            } else if (i<0.8*N){
                while(1){
                    u = (double)(rand())/RAND_MAX;
                    v = (double)(rand())/RAND_MAX;
                    if (sqrt(pow(u-0.8,2) + pow (v-0.8,2)) < 0.2 && u+v>1.336){
                        *x=u;
                        *y=v;
                        break;
                    }
                }
            } else {
                    u = (double)(rand())/RAND_MAX;
                    v = (double)(rand())/RAND_MAX;
                    *x=0.282+0.1*(3.36*u+v);
                    *y=0.382+0.1*(3.36*u-v);
            }
            *x = my_round(*x);
            *y = my_round(*y);
            break;
        }
        for (j=list_i; j<N_rest-1; j++) tmp_list[j]=tmp_list[j+1];
        N_rest--;
    }

    //printf("ERROR: Coordinates were not set for some reason! (neuron %d)\n",i);
    //getchar();
    //exit(1);
    free(tmp_list);


}


//*************Main body***********************************************************************

int main(int argc, char *argv[]){

    //declaration of all needed variables

    struct neuron **neurons = malloc(MAXIMUM_NUMBER_OF_NEURONS*sizeof(struct neuron*));
    int i,j, CONN_PER_NEURON_COUNT, layer_type=1, COORDINATES_LOAD_TYPE=0, SYNAPTIC_DATA_LOAD_TYPE=0, STIMULATION_DATA_LOAD_TYPE=0, W_OUTPUT_COUNTER = 1, TOPOLOGY_LOAD_TYPE = 0, N_SP, N_ACTIVE_CONN, N_con=0, BG_CURRENT_NOISE_MODE=0;
    int burst_flag, burst_counter, eof_flag, eof_flag2, INH_COUNT, EXC_COUNT, N_COUNT, FORCE_BINOMIAL_TOPOLOGY=0, N_PM = 0;
    struct synapse *tmp_synapse, *tmp_synapse2;
    struct inc_synapse *tmp_input_synapse, *tmp_input_synapse2;
    struct delay_timer *tmp_timer, *tmp_timer2;
    double simulation_time, last_burst_time;
    double tmp_I, tmp_prob, prob;
    double averaging_timer, burst_detection_timer, synapse_resource_timer, w_output_timer = 1000;
    double old_u, avg_activity, activity;
    int max_spikes = -1, pre_syn_type, post_syn_type, pre_num, post_num, pre_num2, post_num2; //These temporary variables are used to read parameters from files
    double tmp_A, tmp_tau_I, tmp_tau_rec, tmp_tau_facil, tmp_U, SPIKE_SPEED, burst_threshold;              //parameters from files
    double avg_x, avg_y, avg_z, avg_u, x_init = 0.98, y_init = 0.01, z_init = 0.01, spike_resource;
    char c, tmp_str[100], tmp_str_2[100];



    //declaration of input\output files

    FILE *output_raster;

    FILE *output_file_spiking_activity;
    FILE *output_amount_of_active_connections;

    FILE *output_x;
    FILE *output_y;
    FILE *output_z;
    FILE *output_u;
    FILE *output_M;

    FILE *output_w;
    FILE *output_w_initial;

    FILE *exc_I_distribution;
    FILE *inh_I_distribution;
    FILE *p_sp_distribution;
    FILE *synaptic_parameters_distribution;

    FILE *output_connections;
    FILE *output_connections_distribution;
    FILE *output_coordinates;

    FILE *output_info;

    FILE *output_IBI;
    FILE *output_burst_times;
    FILE *output_lifetimes;
    FILE *output_spike_resource;

    FILE *input = fopen("input.txt", "r");
    FILE *input_connections = NULL;
    FILE *input_coordinates = NULL;
    FILE *input_exc_I_distribution = NULL;
    FILE *input_inh_I_distribution = NULL;
    FILE *input_p_sp_distribution = NULL;
    FILE *input_synaptic_parameters_distribution = NULL;


    //****************Reading input data************************
    // setting default values in case we don't find them in a file
    int N = 50000;        //number of neurons
    double dt = 0.1;   //simulation time step
    double INH_NEURONS_FRACTION = 0.0;  //fraction of inhibitory neurons among N

    int NEURON_MODEL = 1; //membrane potential model of a neuron. 1 for Leaky-integrate-and-fire(LIF), 0 for Perfect Integrate-and-Fire(PIF)

    //LIF model default parameters
    double TAU_M = 20;  //ms, membrane potential V relaxation constant
    double R_IN = 1;   //GOhm, membrane resistance
    double V_REST = 0; //mV, resting potential
    double V_TH = 15;   //mV, threshold potential
    double V_RESET = 13.5;   //mV, threshold potential
    double TAU_REF_EXC = 3;  //ms, refractory period of excitatory neurons
    double TAU_REF_INH = 2;  //ms, refractory period of inhibitory neurons

    //PIF model addition
    double C_M = 20; //pF, membrane capacitance


    //T_M model default parameters (mean values of synaptic parameters)
    double Avg_A[2][2] = {{-72, -72}, {54, 54}}; //pA
    double Avg_U[2][2] = {{0.04, 0.04}, {0.5, 0.5}};
    double Avg_tau_rec[2][2] = {{100, 100}, {800, 800}}; //ms
    double Avg_tau_facil[2][2] = {{1000, 1000}, {0, 0}}; //ms
    double TAU_I = 3; //ms

    double SYN_RESOURCES_OUTPUT_PERIOD = 10; //ms

    int STIMULATION_TYPE = STIM_TYPE_I_BG_GAUSSIAN; //stimulation type, 1 for background currents, 2 for spontaneous spiking

    //stimulation default parameters for background currents
    double I_BG_MEAN = 7.7; //pA
    double I_BG_SD = 4.0; //pA
    double I_BG_MIN = 0; //pA
    double I_BG_MAX = 20; //pA

    double I_BG_1 = 2; //pA
    double I_BG_2 = 37; //pA
    double FRACTION_OF_NEURONS_WITH_I_BG_1 = 0.75;

    double I_BG_NOISE_SD = 0;

    //default parameters for the case of stochastic stimulation
    double P_SP_MEAN = 0.0005;
    double P_SP_SD = 0.0;
    double P_SP_MIN = 0.0;
    double P_SP_MAX = 0.001;

    double P_SP_1 = 0.005;
    double P_SP_2 = 0.00001;
    double FRACTION_OF_NEURONS_WITH_P_SP_1 = 0.25;

    //STDP parameters
    int STDP_status = STDP_IS_OFF;
    double A_plus = 0.01;
    double A_minus = 0.01;
    double w_initial = 0.5;
    double tau_plus_corr = 20;
    double tau_minus_corr = 20;
    double W_OUTPUT_PERIOD = 500;


    //topology parameters
    double lambda = 0.01; //characteristic length of a connection within an exponential distribution in spatially-dependent topology
    double tau_delay = 0.2; //base axonal delay
    SPIKE_SPEED = 0.2; //speed of spike propagating along the axon. All connections are considered straight lines.
    double max_conn_length = 1.415; //maximum length of a connection in spatially-dependent topology

    double P_CON = 0.1; //probability of connecting two neurons in case of BINOMIAL topology

    //Homeostasis parameters
    int HOMEOSTASIS_status = 0;
    double b = 1/(double)N;
    double M=1, M_max = 1;

    //stim disabling protocol
    int stim_disable_protocol=0; //(1 for turning stimulation off when overcoming the activity threshold
                                    // 2 for turning stimulation off at a certain time)
    int pb_number=2;            //number of burst to which the stimulation disabling protocol is applied
    double activity_threshold=0.2; //activity threshold of the pb_number-th burst, at which the stimulation disabling protocol is applied
    double stim_off_time=1400; //(ms - moment when to turn off the stimulation in case stim_disable_protocol = 2)

    double inh_off_time=-1; //ms - moment in time when to set membrane potential V to V_rest until the simulation end. If set to -1 - mechanism is turend OFF

    //simulation parameters
    double AVG_TIME = 2; //ms, activity averaging time
    int SIM_TIME = 10000; //ms, simulation duration
    burst_threshold = 0.2; //activity level for detecting a burst


    //a complete list of parameters' names to read from input file
    struct parameter value_array[PARAMETER_N] = {
        {"N", &N},
        {"dt", &dt},

        {"Inh_neurons_fraction", &INH_NEURONS_FRACTION},

        {"lambda", &lambda},
        {"spike_speed", &SPIKE_SPEED},
        {"tau_delay", &tau_delay},
        {"max_conn_length", &max_conn_length},

        {"P_con", &P_CON},

        {"neuron_model", &NEURON_MODEL},

        {"TAU_M", &TAU_M},
        {"R_IN", &R_IN},
        {"C_M", &C_M},

        {"V_TH", &V_TH},
        {"V_RESET", &V_RESET},
        {"V_REST", &V_REST},

        {"TAU_REF_EXC", &TAU_REF_EXC},
        {"TAU_REF_INH", &TAU_REF_INH},

        {"Aee", &Avg_A[1][1]},
        {"Aei", &Avg_A[1][0]},
        {"Aie", &Avg_A[0][1]},
        {"Aii", &Avg_A[0][0]},

        {"Uee", &Avg_U[1][1]},
        {"Uei", &Avg_U[1][0]},
        {"Uie", &Avg_U[0][1]},
        {"Uii", &Avg_U[0][0]},

        {"tau_rec_ee", &Avg_tau_rec[1][1]},
        {"tau_rec_ei", &Avg_tau_rec[1][0]},
        {"tau_rec_ie", &Avg_tau_rec[0][1]},
        {"tau_rec_ii", &Avg_tau_rec[0][0]},

        {"tau_facil_ee", &Avg_tau_facil[1][1]},
        {"tau_facil_ei", &Avg_tau_facil[1][0]},
        {"tau_facil_ie", &Avg_tau_facil[0][1]},
        {"tau_facil_ii", &Avg_tau_facil[0][0]},

        {"tau_I", &TAU_I},

        {"Stimulation_type", &STIMULATION_TYPE},

        {"I_bg_mean", &I_BG_MEAN},
        {"I_bg_sd", &I_BG_SD},
        {"I_bg_min", &I_BG_MIN},
        {"I_bg_max", &I_BG_MAX},

        {"I_bg_1", &I_BG_1},
        {"I_bg_2", &I_BG_2},
        {"Fraction_of_neurons_with_I_bg_1",&FRACTION_OF_NEURONS_WITH_I_BG_1},

        {"P_sp_mean", &P_SP_MEAN},
        {"P_sp_sd", &P_SP_SD},
        {"P_sp_min", &P_SP_MIN},
        {"P_sp_max", &P_SP_MAX},

        {"P_sp_1", &P_SP_1},
        {"P_sp_2", &P_SP_2},
        {"Fraction_of_neurons_with_P_sp_1",&FRACTION_OF_NEURONS_WITH_P_SP_1},

        {"spatial_layer_type", &layer_type},
        {"force_binomial_topology", &FORCE_BINOMIAL_TOPOLOGY},

        {"stim_disable_protocol",&stim_disable_protocol},
        {"pb_number",&pb_number},
        {"activity_threshold",&activity_threshold},
        {"stim_off_time",&stim_off_time},

        {"inh_off_time",&inh_off_time},

        {"x_init", &x_init},
        {"y_init", &y_init},
        {"z_init", &z_init},

        {"max_spikes_per_neuron", &max_spikes},

        {"use_saved_topology",&TOPOLOGY_LOAD_TYPE},
        {"use_saved_stimulation_data",&STIMULATION_DATA_LOAD_TYPE},
        {"use_saved_synaptic_parameters",&SYNAPTIC_DATA_LOAD_TYPE},
        {"use_saved_coordinates",&COORDINATES_LOAD_TYPE},

        {"STDP_status",&STDP_status},
        {"W_OUTPUT_PERIOD",&W_OUTPUT_PERIOD},

        {"I_bg_noise_mode",&BG_CURRENT_NOISE_MODE},
        {"I_bg_noise_sd",&I_BG_NOISE_SD},

        {"A_plus",&A_plus},
        {"A_minus",&A_minus},
        {"w_initial",&w_initial},
        {"tau_plus_corr",&tau_plus_corr},
        {"tau_minus_corr",&tau_minus_corr},

        {"Homeostasis_status",&HOMEOSTASIS_status},
        {"b",&b},
        {"M_max",&M_max},

        {"burst_threshold", &burst_threshold},

        {"AVG_TIME", &AVG_TIME},
        {"SIM_TIME", &SIM_TIME},

    };




    srand(time(0));  //taking value as a seed of a random function

    printf("********************NeuroSimTM-2.1*****************\n");
    printf("version 16.02.2017\n\n");

    printf("Reading input file...\n");


    if (input == NULL){
        printf("****IMPORTANT! Failed to open input file! using default parameters***********\n");
    } else readfile (input, value_array);


    printf("Ready to start simulation with following parameters:\n");
    printf("\n>>>> GENERAL <<<<\n");
    printf("N = %d\n",N);
    printf("INHIBITORY neurons fraction: %g\n",INH_NEURONS_FRACTION);



    printf("\n\n>>>> NETWORK MODEL <<<<\n");
    if (NEURON_MODEL == NEURON_MODEL_PERFECT_INTEGRATE_AND_FIRE){
        printf("Neuron Model: Perfect Integrate-and-Fire (PIF)\n");
        printf("C_m = %g pF\n",C_M);
    } else {
        printf("Neuron Model: Leaky Integrate-and-Fire (LIF)\n");
        printf("TAU_M = %g ms\n",TAU_M);
        printf("R_in = %g GOhm\n",R_IN);
    }

    printf("V_REST = %g mV\n",V_REST);
    printf("V_TH = %g mV\n",V_TH);
    printf("TAU_REF_EXC = %g ms\n", TAU_REF_EXC);
    if (INH_NEURONS_FRACTION > 0)
        printf("TAU_REF_INH = %g ms\n", TAU_REF_INH);

    printf("\nx_init = %g\n", x_init);
    printf("y_init = %g\n", y_init);
    printf("z_init = %g\n", z_init);


    if ( fabs((x_init + y_init + z_init) - 1) >  0.000001){
        printf("\n\nERROR: x + y + z = %g, though it has to be %g\n", x_init + y_init + z_init, 1.0);
        fflush(stdin);
        getchar();
        exit(1);
    }

    printf("\n\n>>>> STIMULATION <<<<");

    if (STIMULATION_TYPE == STIM_TYPE_I_BG_GAUSSIAN || STIMULATION_TYPE == STIM_TYPE_I_BG_TWO_VALUES){
        printf("\nStimulation type: Background currents\n");
        if (STIMULATION_DATA_LOAD_TYPE == 1){
            printf("\nBackground currents will be read from file\n");
        } else {
            printf("Background currents will be generated with following properties\n");
            if (STIMULATION_TYPE == STIM_TYPE_I_BG_GAUSSIAN){
                printf("I_bg_mean = %g pA\n",I_BG_MEAN);
                printf("I_bg_sd = %g pA\n",I_BG_SD);
                printf("I_bg_min = %g pA\n",I_BG_MIN);
                printf("I_bg_max = %g pA\n",I_BG_MAX);
            } else {
                printf("I_bg_1 = %g pA\n",I_BG_1);
                printf("I_bg_2 = %g pA\n", I_BG_2);
                printf("Fraction_of_neurons_with_I_bg_1 = %g\n", FRACTION_OF_NEURONS_WITH_I_BG_1);
            }
        }
    } else {
        printf("\nStimulation type: Spontaneous spiking\n");
        if (STIMULATION_DATA_LOAD_TYPE == 1){
            printf("\nP_sp probabilities will be read from file\n");
        } else {
            if (STIMULATION_TYPE == STIM_TYPE_P_SP_GAUSSIAN){
                printf("\nP_sp probabilities will be generated with following properties\n");
                printf("P_SP_mean = %g\n",P_SP_MEAN);
                printf("P_SP_sd = %g\n",P_SP_SD);
                printf("P_SP_min = %g\n",P_SP_MIN);
                printf("P_SP_max = %g\n",P_SP_MAX);
            } else {
                printf("p_sp_1 = %g\n",P_SP_1);
                printf("p_sp_2 = %g\n", P_SP_2);
                printf("Fraction_of_neurons_with_p_sp_1 = %g\n", FRACTION_OF_NEURONS_WITH_P_SP_1);
            }
        }
    }

    printf("\n\n>>>> TOPOLOGY <<<<\n");
    printf("Layer type: ");
    switch(layer_type){

    case LAYER_TYPE_BINOMIAL:
        printf("BINOMIAL\n");
        break;

    case LAYER_TYPE_UNIFORM:
        printf("UNIFORM\n");
        break;

    case LAYER_TYPE_SQUARE_LATTICE:
        printf("SQUARE_LATTICE\n");
        break;

    case LAYER_TYPE_BELL:
        printf("BELL\n");
        break;

    case LAYER_TYPE_RAMP:
        printf("RAMP\n");
        break;

    case LAYER_TYPE_DOUBLE_RAMP:
        printf("DOUBLE RAMP\n");
        break;

    case LAYER_TYPE_BARBELL:
        printf("BARBELL\n");
        break;

    case LAYER_TYPE_STRIPED:
        printf("STRIPED\n");
        break;
    }

    if (FORCE_BINOMIAL_TOPOLOGY == 1){
        printf("\nFORCED TO BINOMIAL TOPOLOGY\n");
        printf("P_con = %g\n",P_CON);
    }

    if (layer_type == LAYER_TYPE_BINOMIAL)
        printf("P_con = %g\n",P_CON);
    else{
        if (FORCE_BINOMIAL_TOPOLOGY != 1) printf("lambda = %g\n",lambda);
        printf("max_conn_length = %g L\n",max_conn_length);
        printf("spike_speed = %g L/ms\n",SPIKE_SPEED);
        printf("tau_delay = %g ms\n",tau_delay);
    }

    if (TOPOLOGY_LOAD_TYPE == READ_TOPOLOGY_FROM_FILE || SYNAPTIC_DATA_LOAD_TYPE == READ_SYNAPTIC_DATA_FROM_FILE ||
            COORDINATES_LOAD_TYPE == READ_COORDINATES_FROM_FILE || STIMULATION_DATA_LOAD_TYPE == READ_STIMULATION_DATA_FROM_FILE){
        printf("\n\n>>>> USING SAVED DATA <<<<");

        if (TOPOLOGY_LOAD_TYPE == READ_TOPOLOGY_FROM_FILE){
            printf("\nTOPOLOGY is being READ FROM FILE\n");
        } else {
                              // in case we don't read topology input file we generate it
            printf("\nTOPOLOGY will be GENERATED for this run\n");
        }


        if (SYNAPTIC_DATA_LOAD_TYPE == READ_SYNAPTIC_DATA_FROM_FILE){

            if (TOPOLOGY_LOAD_TYPE != READ_TOPOLOGY_FROM_FILE){
                printf("\nERROR: Synaptic parameters cannot be read from file if topology is not read from file!\n");
                getchar();
                exit(1);
            }

            printf("SYNAPTIC PARAMETERS distribution is being READ FROM FILE\n");

        } else {
                              // in case we don't read synaptic parameters from file we generate them
            printf("SYNAPTIC PARAMETERS will be GENERATED for this run\n");
        }

        if (COORDINATES_LOAD_TYPE == READ_COORDINATES_FROM_FILE){
            printf("COORDINATES are being READ FROM FILE expecting ");
        } else { printf("COORDINATES will be GENERATED with ");}
            switch(layer_type){

         case LAYER_TYPE_BINOMIAL:
                printf("layer type BINOMIAL\n");
                break;

            case LAYER_TYPE_UNIFORM:
                printf("layer type UNIFORM\n");
                break;

            case LAYER_TYPE_SQUARE_LATTICE:
                printf("layer type SQUARE_LATTICE\n");
                break;

            case LAYER_TYPE_BELL:
                printf("layer type BELL\n");
                break;

            case LAYER_TYPE_RAMP:
                printf("layer type RAMP\n");
                break;

            case LAYER_TYPE_DOUBLE_RAMP:
                printf("layer type DOUBLE RAMP\n");
                break;

            case LAYER_TYPE_BARBELL:
                printf("layer type BARBELL\n");
                break;

            case LAYER_TYPE_STRIPED:
                printf("layer type STRIPED\n");
                break;
        }

        if (STIMULATION_DATA_LOAD_TYPE == READ_STIMULATION_DATA_FROM_FILE){
            printf("\nSTIMULATION DATA is being READ FROM FILE\n");
        } else {
            printf("\nSTIMULATION DATA will be GENERATED for this run\n");
        }
    }

    printf("\n\n>>>> ADDITIONAL MECHANISMS <<<<");

    if (BG_CURRENT_NOISE_MODE == 0)
        printf("\nBG current noise is OFF\n");
    else{
        printf("\nBG current noise is ON\n");
        printf("I_bg_noise_sd = %g pA\n",I_BG_NOISE_SD);
    }

    switch(STDP_status){
    case STDP_IS_OFF:
        printf("STDP status is OFF\n");
        break;
    case STDP_IS_ADDITIVE:
        printf("STDP status is ADDITIVE\n");
        break;
    case STDP_IS_MULTIPLICATIVE:
        printf("STDP status is MULTIPLICATIVE\n");
        break;
    }

    if (HOMEOSTASIS_status == HOMEOSTASIS_IS_ON)
        printf("HOMEOSTASIS is ON\n");
    else
        printf("HOMEOSTASIS is OFF\n");

    if (stim_disable_protocol != 0)
            printf("STIMULATION DISABLING protocol is ON\n");
    if (max_spikes < 0){
        printf("LIMITED SPIKE RESOURCE mode is OFF\n");
    } else {
        printf("LIMITED SPIKE RESOURCE mode is ON\nmax_spikes_per_neurons = %d",max_spikes);
    }

    if (inh_off_time != -1)
            printf("\nInhibitory neurons will be forced to V_REST on %g ms\n",inh_off_time);



    printf("\n\n>>>> SIMULATION PARAMETERS <<<<");
    printf("\nBurst detection theshold: %g",burst_threshold);
    printf("\nTime step dt = %g ms",dt);
    printf("\nAveraging time: %g ms\n",AVG_TIME);
    printf("Simulation time: %d ms\n",SIM_TIME);

    if (input != NULL) fclose(input);

    //we let user to check if everything is alright before program creates output files and starts the simulation

    printf("\nStart?(Y/N)\n");
    scanf("%c",&c);
    if (c == 'N' || c == 'n'){
        printf("Aborting simulation!\nPress Enter...\n");
        exit(1);                                        //exiting program if user thinks something's gone wrong
    }
    fflush(stdin);


    //*********************Opening output files********************

    if (TOPOLOGY_LOAD_TYPE == READ_TOPOLOGY_FROM_FILE){
        input_connections = fopen("saved_connections.txt","r");

        if (input_connections == NULL){
            printf("ERROR: Failed to open given topology input file!\n");
            getchar();
            exit(1);
        }
    }

    if (SYNAPTIC_DATA_LOAD_TYPE == READ_SYNAPTIC_DATA_FROM_FILE){
        input_synaptic_parameters_distribution = fopen("saved_synaptic_parameters_distribution.txt","r");

        if (input_synaptic_parameters_distribution == NULL){
            printf("ERROR: Failed to open given synaptic parameters distribution input file!\n");
            getchar();
            exit(1);
        }
    }

    output_raster = fopen("raster.txt","w");

    output_file_spiking_activity = fopen("activity.txt","w");

    output_IBI = fopen("IBI.txt","w");
    output_burst_times = fopen("burst_times.txt","w");

    output_amount_of_active_connections = fopen("active_connections.txt","w");

    system("mkdir syn_resources_dynamics");
    output_x = fopen("syn_resources_dynamics/x.txt","w");
    output_y = fopen("syn_resources_dynamics/y.txt","w");
    output_z = fopen("syn_resources_dynamics/z.txt","w");
    output_u = fopen("syn_resources_dynamics/u.txt","w");

    output_w_initial = fopen("w_initital.txt","w");

    output_M = fopen("M.txt","w");

    if (max_spikes != -1){
        output_lifetimes = fopen("lifetimes.txt","w");
        output_spike_resource = fopen("spiking_resource.txt","w");
    }

    if (STDP_status != STDP_IS_OFF){
        system("mkdir stdp_w_dynamics");
        if (layer_type == LAYER_TYPE_BARBELL){
            system("mkdir \stdp_w_dynamics\\R");
            system("mkdir \stdp_w_dynamics\\L");
        }
    }


    if (STIMULATION_TYPE == STIM_TYPE_I_BG_GAUSSIAN || STIMULATION_TYPE == STIM_TYPE_I_BG_TWO_VALUES){
        exc_I_distribution = fopen("exc_I_distribution.txt","w");
        inh_I_distribution = fopen("inh_I_distribution.txt","w");
    } else
        p_sp_distribution = fopen("p_sp_distribution.txt","w");

    output_connections = fopen("connections.txt","w");
    output_connections_distribution = fopen("connections_per_neuron.txt","w");
    synaptic_parameters_distribution = fopen("synaptic_parameters_distribution.txt","w");
    output_coordinates = fopen("coordinates.txt","w");


    if (output_raster == NULL || output_file_spiking_activity == NULL){
        printf("\nERROR: Failed to open raster or activity output file!\n");
        getchar();
        exit(1);
    }

    if (output_coordinates == NULL){

        printf("\nERROR: Failed to open coordinates output file!\n");
        getchar();
        exit(1);
    }

    if (COORDINATES_LOAD_TYPE == READ_COORDINATES_FROM_FILE){

        if ((input_coordinates = fopen("saved_coordinates.txt","r")) == NULL){
            printf("ERROR: Failed to open saved_coordinates.txt!");
            getchar();
            exit(1);
        }
    }


    //****************Neuron array initialization***************************************************************************

    INH_COUNT = 0;
    EXC_COUNT = 0;

    printf("Creating neurons...\n");

    //we allocate memory for every neuron and set initial values

    for (i=0;i<N;i++){
        neurons[i] = malloc(sizeof(struct neuron));
        neurons[i]->V = V_REST;
        neurons[i]->V_rest = V_REST;
        neurons[i]->V_reset = V_RESET;
        neurons[i]->ref_time_left = 0;
        neurons[i]->id = i;
        neurons[i]->next = NULL;
        neurons[i]->num_of_outcoming_connections = 0;
        neurons[i]->num_spikes = 0;

        neurons[i]->last_spiked_at = -100;

        neurons[i]->x = -1;
        neurons[i]->y = -1;

        //setting type of a neuron, excitatory or inhibitory
        //also here we detect mismatches between topology and synaptic_parameters_distribution files
        if (TOPOLOGY_LOAD_TYPE == READ_TOPOLOGY_FROM_FILE){
            while ((eof_flag = fscanf(input_connections,"%d %d %lf",&pre_num,&post_num,&tmp_prob)) != EOF && abs(pre_num) != i);

            if (eof_flag == EOF && i != N-1){
                printf("ERROR: Neuron amount in saved_connections.txt mismatch!\n");
                getchar();
                exit(1);
            }

            if (pre_num >= 0) {
                neurons[i]->type = TYPE_EXCITATORY;
                EXC_COUNT++;
            } else {
                neurons[i]->type = TYPE_INHIBITORY;
                INH_COUNT++;
            }
        } else {
            prob = (double)(rand())/RAND_MAX;
            if (prob >= INH_NEURONS_FRACTION || i == 0) {
                neurons[i]->type = TYPE_EXCITATORY;
                EXC_COUNT++;
            } else {
                neurons[i]->type = TYPE_INHIBITORY;
                INH_COUNT++;
            }
        }

        //setting background current or p_sp for a neuron in case we don't read them from file

        if (STIMULATION_DATA_LOAD_TYPE == GENERATE_STIMULATION_DATA){

            switch(STIMULATION_TYPE){
                case STIM_TYPE_I_BG_GAUSSIAN:
                    neurons[i]->I_b = gauss (I_BG_MEAN, I_BG_SD,I_BG_MIN,I_BG_MAX);
                    if (neurons[i]->I_b > V_TH/R_IN) N_PM++;
                    if (neurons[i]->type == TYPE_EXCITATORY)
                        fprintf(exc_I_distribution,"%d %g\n",i,neurons[i]->I_b);
                    else
                        fprintf(inh_I_distribution,"%d %g\n",-i,neurons[i]->I_b);
                    break;

                case STIM_TYPE_I_BG_TWO_VALUES:
                    if (i < FRACTION_OF_NEURONS_WITH_I_BG_1*(double)N)
                         neurons[i]->I_b = I_BG_1;
                    else neurons[i]->I_b = I_BG_2;
                    if (neurons[i]->I_b > V_TH/R_IN) N_PM++;
                    if (neurons[i]->type == TYPE_EXCITATORY)
                        fprintf(exc_I_distribution,"%d %g\n",i,neurons[i]->I_b);
                    else
                        fprintf(inh_I_distribution,"%d %g\n",-i,neurons[i]->I_b);
                    break;

                case STIM_TYPE_P_SP_GAUSSIAN:
                    neurons[i]->I_b = 0;
                    neurons[i]->p_sp = gauss (P_SP_MEAN, P_SP_SD, P_SP_MIN, P_SP_MAX);
                    fprintf(p_sp_distribution,"%d %g\n",i,neurons[i]->p_sp);
                    break;

                case STIM_TYPE_P_SP_TWO_VALUES:
                    neurons[i]->I_b = 0;
                    if (i < FRACTION_OF_NEURONS_WITH_P_SP_1*(double)N)
                        neurons[i]->p_sp = P_SP_1;
                    else neurons[i]->p_sp = P_SP_2;
                    fprintf(p_sp_distribution,"%d %g\n",i,neurons[i]->p_sp);
                    break;
            }
        }


        //preparing shortcuts of connection lists for generating connections
        neurons[i]->out_conn = malloc(sizeof(struct synapse));
        neurons[i]->out_conn->next = NULL;
        neurons[i]->out_conn->post_id = -1;

        neurons[i]->in_conn = malloc(sizeof(struct inc_synapse));
        neurons[i]->in_conn->syn_pointer = NULL;
        neurons[i]->in_conn->next = NULL;


    }

    // setting coordinates according to the layer type
    if (COORDINATES_LOAD_TYPE == GENERATE_COORDINATES)
        setCoordinates(neurons,layer_type,N);

    for (i=0; i<N; i++){

        if (COORDINATES_LOAD_TYPE == READ_COORDINATES_FROM_FILE) {
            eof_flag = fscanf(input_coordinates,"%lf %lf", &neurons[i]->x,&neurons[i]->y);
            if (eof_flag == EOF && i != N-1){
            printf("ERROR: Neuron amount in saved_coordinates.txt mismatch!\n");
            getchar();
            exit(1);
            }
        }

        fprintf(output_coordinates,"%g %g\n",neurons[i]->x,neurons[i]->y);
    }

    //reopening saved_connections.txt as we've just passed through the file to set exc\inh types, and now we set it
    //to get connections themselves

    if (TOPOLOGY_LOAD_TYPE == READ_TOPOLOGY_FROM_FILE){
        fclose(input_connections);
        input_connections = fopen("saved_connections.txt","r");
    }


    //****After creating neurons we make one more pass to set background current or p_sp distributions from files**************

    if (STIMULATION_DATA_LOAD_TYPE == READ_STIMULATION_DATA_FROM_FILE){
        if (STIMULATION_TYPE == STIM_TYPE_I_BG_GAUSSIAN  || STIMULATION_TYPE == STIM_TYPE_I_BG_TWO_VALUES){

            input_exc_I_distribution = fopen("saved_exc_I_distribution.txt","r");
            input_inh_I_distribution = fopen("saved_inh_I_distribution.txt","r");

            if (input_exc_I_distribution == NULL){
                printf("ERROR: Failed to open current distribution files!\n");
                getchar();
                exit(1);
            }

            N_COUNT = 0;

            while(fscanf(input_exc_I_distribution,"%d %lf",&i,&tmp_I) != EOF){
                i = abs(i);
                if (i>=N){
                    N_COUNT = 0;
                    break;
                }
                N_COUNT++;

                if (TOPOLOGY_LOAD_TYPE == READ_TOPOLOGY_FROM_FILE && neurons[i]->type != TYPE_EXCITATORY){
                    printf("CURRENT TYPE %d\n",neurons[i]->type);
                    printf("ERROR: Input data on exc\\inh type from input_exc_I_distribution.txt for neuron %d mismatch!\n",i);
                    getchar();
                    exit(1);
                }

                neurons[i]->type = TYPE_EXCITATORY;
                neurons[i]->I_b = tmp_I;
                if (tmp_I > V_TH/R_IN) N_PM++;
                fprintf(exc_I_distribution,"%d %g\n",i,tmp_I);
            }

            if (INH_NEURONS_FRACTION > 0 && input_inh_I_distribution != NULL){
                while(fscanf(input_inh_I_distribution,"%d %lf",&i,&tmp_I) != EOF){
                    i = abs(i);
                    if (i>=N){
                     N_COUNT = 0;
                        break;
                    }
                    N_COUNT++;

                    if (TOPOLOGY_LOAD_TYPE == READ_TOPOLOGY_FROM_FILE && neurons[i]->type != TYPE_INHIBITORY){
                        printf("ERROR: Input data on exc\\inh type from input_inh_I_distribution.txt for neuron %d mismatch!\n",i);
                        getchar();
                        exit(1);
                    }

                    neurons[i]->type = TYPE_INHIBITORY;
                    neurons[i]->I_b = tmp_I;
                    if (tmp_I > V_TH/R_IN) N_PM++;
                    fprintf(inh_I_distribution,"%d %g\n",-i,tmp_I);
                }
            } else if (INH_NEURONS_FRACTION > 0){
                printf("ERROR: Failed to open inh_I_distribution.txt though there are inhibitory neurons!\n");
                getchar();
                exit(1);
            }

            if (N_COUNT != N) {
                printf("ERROR: Neuron amount in stimulation data mismatch!\n"); //Exit if mismatch detected
                getchar();
                exit(1);
            }

        } else {

            input_p_sp_distribution = fopen("saved_p_sp_distribution.txt","r");

            if (input_p_sp_distribution == NULL){
                printf("ERROR: Failed to open p_sp distribution file!\n");
                getchar();
                exit(1);
            }

            N_COUNT = 0;

            while(fscanf(input_p_sp_distribution,"%d %lf",&i,&tmp_I) != EOF){
                i = abs(i);
                if (i>=N){
                    N_COUNT = 0;
                    break;
                }
                N_COUNT++;
                neurons[i]->I_b = 0;
                neurons[i]->p_sp = tmp_I;
                if (neurons[i]->type == TYPE_EXCITATORY)
                    fprintf(p_sp_distribution,"%d %g\n",i,tmp_I);
                else
                    fprintf(p_sp_distribution,"%d %g\n",-i,tmp_I);
            }

            if (N_COUNT != N) {
                printf("ERROR: Neuron amount in stimulation data mismatch!\n");
                getchar();
                exit(1);
            }

        }
    }


    //closing files we don't need anymore
    fclose(output_coordinates);
    if (COORDINATES_LOAD_TYPE == READ_COORDINATES_FROM_FILE)
        fclose(input_coordinates);

    if (STIMULATION_DATA_LOAD_TYPE){
        if (STIMULATION_TYPE == STIM_TYPE_I_BG_GAUSSIAN || STIMULATION_TYPE == STIM_TYPE_I_BG_TWO_VALUES){
            fclose(input_exc_I_distribution);
            if (INH_NEURONS_FRACTION > 0)
                fclose(input_inh_I_distribution);
        } else
            fclose(input_p_sp_distribution);
    }

    if (STIMULATION_TYPE == STIM_TYPE_I_BG_GAUSSIAN || STIMULATION_TYPE == STIM_TYPE_I_BG_TWO_VALUES){
        fclose(exc_I_distribution);
        fclose(inh_I_distribution);
    } else
        fclose(p_sp_distribution);


    //***************Establishsing connections******************************************************************************
    printf("Establishing neuron connections...\n");


    N_COUNT = 0;

    //in case we have topology input file
    if (TOPOLOGY_LOAD_TYPE == READ_TOPOLOGY_FROM_FILE)
        fscanf(input_connections,"%d %d %lf",&pre_num,&post_num,&tmp_prob);

    for (i=0;i<N;i++){
        tmp_synapse = neurons[i]->out_conn;  //a pointer that runs through the list of outcoming connections

        pre_syn_type = neurons[i]->type;
        CONN_PER_NEURON_COUNT = 0;



        //looking over every other neuron in the network and checking whether we create a connection
        //(if we read topology from file we still check if given probability value tmp_prob satisfies our connection occurrence condidtion)
        for (j=0;j<N;j++){

            if (j != i && ((TOPOLOGY_LOAD_TYPE == READ_TOPOLOGY_FROM_FILE && i == abs(pre_num) && j == abs(post_num)) || (TOPOLOGY_LOAD_TYPE == GENERATE_TOPOLOGY))){

                //reading given probability value from file or generating it
                if (TOPOLOGY_LOAD_TYPE == READ_TOPOLOGY_FROM_FILE){
                    prob = tmp_prob;
                    if (post_num < 0){

                        if (STIMULATION_DATA_LOAD_TYPE == READ_STIMULATION_DATA_FROM_FILE && neurons[j]->type == TYPE_EXCITATORY) { //checking postsynaptic neuron type
                            printf("ERROR: Input data on exc\\inh type for neuron %d mismatch!\n",j);
                            getchar();
                            exit(1);
                        }

                        neurons[j]->type = TYPE_INHIBITORY;
                    } else {
                        if (STIMULATION_DATA_LOAD_TYPE == READ_STIMULATION_DATA_FROM_FILE && neurons[j]->type == TYPE_INHIBITORY) {
                            printf("ERROR: Input data on exc\\inh type for neuron %d mismatch!\n",j);
                            getchar();
                            exit(1);

                        }

                        neurons[j]->type = TYPE_EXCITATORY;
                    }


                } else
                    prob = (double)(rand())/RAND_MAX;

                //****Checking probability whether we should create a connection (synapse)
                if ((layer_type != LAYER_TYPE_BINOMIAL && FORCE_BINOMIAL_TOPOLOGY != 1 && prob < my_theta(max_conn_length - sqrt(pow(neurons[i]->x - neurons[j]->x,2) + pow(neurons[i]->y - neurons[j]->y,2)))*exp(-sqrt(pow(neurons[i]->x - neurons[j]->x,2) + pow(neurons[i]->y - neurons[j]->y,2))/lambda)) ||
                        (layer_type == LAYER_TYPE_BINOMIAL && prob <= P_CON) || (FORCE_BINOMIAL_TOPOLOGY == 1 && prob <= P_CON)){

                    //counting for statistics
                    N_con++;
                    CONN_PER_NEURON_COUNT++;                    

                    neurons[i]->num_of_outcoming_connections++;
                    post_syn_type = neurons[j]->type;

                    //setting all initial values

                    tmp_synapse->pre_id = i;
                    tmp_synapse->post_id = j;
                    tmp_synapse->l = sqrt(pow(neurons[i]->x - neurons[j]->x,2) + pow(neurons[i]->y - neurons[j]->y,2));
                    if (layer_type == LAYER_TYPE_BINOMIAL)
                        tmp_synapse->tau_delay = 0;
                    else
                        tmp_synapse->tau_delay = tau_delay + tmp_synapse->l/SPIKE_SPEED;
                    tmp_synapse->timers = NULL;

                    if (STDP_status){
                        tmp_synapse->w = w_initial;
                        fprintf(output_w_initial,"%g\n",tmp_synapse->w);
                    } else tmp_synapse->w = 0.5;


                    if (SYNAPTIC_DATA_LOAD_TYPE == READ_SYNAPTIC_DATA_FROM_FILE){
                        while((eof_flag2 = fscanf(input_synaptic_parameters_distribution,"%d %d %lf %lf %lf %lf %lf",&pre_num2,&post_num2,&tmp_A,&tmp_U,&tmp_tau_I,&tmp_tau_rec,&tmp_tau_facil)) != EOF && !(abs(pre_num) == abs(pre_num2) && abs(post_num) == abs(post_num2)));
                        if (SGN(pre_num) != SGN(tmp_A) || tmp_U < 0 || tmp_tau_facil < 0 || tmp_tau_I < 0 || tmp_tau_rec < 0){
                            printf("ERROR: Synaptic data from saved_synaptic_parameters_distribution.txt mismatch!\nFailed on neuron %d\n\n (Check if U,tau_I,tau_rec,tau_facil>0, Sign of A is correct)\n",pre_num);
                            getchar();
                            exit(1);
                        }
                        if (eof_flag2 == EOF && i != N-1){
                            printf("ERROR: Neuron parameters in saved_synaptic_parameters_distribution.txt mismatch!\n");
                            getchar();
                            exit(1);
                        }

                        tmp_synapse->A = tmp_A;
                        tmp_synapse->U = tmp_U;
                        tmp_synapse->tau_I = tmp_tau_I;
                        tmp_synapse->tau_rec = tmp_tau_rec;
                        tmp_synapse->tau_facil = tmp_tau_facil;

                    } else {

                        if (Avg_A[pre_syn_type][post_syn_type] > 0)
                            tmp_synapse->A = gauss(Avg_A[pre_syn_type][post_syn_type], 0.5*Avg_A[pre_syn_type][post_syn_type], 0, 4*Avg_A[pre_syn_type][post_syn_type]);
                        else
                            tmp_synapse->A = gauss(Avg_A[pre_syn_type][post_syn_type], -0.5*Avg_A[pre_syn_type][post_syn_type], 4*Avg_A[pre_syn_type][post_syn_type], 0);

                        tmp_synapse->U = gauss(Avg_U[pre_syn_type][post_syn_type], 0.5*Avg_U[pre_syn_type][post_syn_type], 0, MIN(1, 4*Avg_U[pre_syn_type][post_syn_type]));

                        tmp_synapse->tau_I = TAU_I;
                        tmp_synapse->tau_rec = gauss(Avg_tau_rec[pre_syn_type][post_syn_type], 0.5*Avg_tau_rec[pre_syn_type][post_syn_type], dt, 4*Avg_tau_rec[pre_syn_type][post_syn_type]);

                        if (Avg_tau_facil[pre_syn_type][post_syn_type] == 0) tmp_synapse->tau_facil = 0;
                        else tmp_synapse->tau_facil = gauss(Avg_tau_facil[pre_syn_type][post_syn_type], 0.5*Avg_tau_facil[pre_syn_type][post_syn_type], dt, 4*Avg_tau_facil[pre_syn_type][post_syn_type]);
                    }

                    fprintf(synaptic_parameters_distribution,"%d %d %g %g %g %g %g\n",i,j,tmp_synapse->A,tmp_synapse->U,tmp_synapse->tau_I,tmp_synapse->tau_rec,tmp_synapse->tau_facil);
                    fflush(synaptic_parameters_distribution);



                    tmp_synapse->x = x_init;
                    tmp_synapse->y = y_init;
                    tmp_synapse->z = z_init;
                    tmp_synapse->u = tmp_synapse->U;

                    //adding this synapse in a list if incoming connections in a list of the postsynaptic neuron
                    tmp_input_synapse = neurons[j]->in_conn;
                    while (tmp_input_synapse->syn_pointer != NULL)
                        tmp_input_synapse = tmp_input_synapse->next;
                    tmp_input_synapse->syn_pointer = tmp_synapse;
                    tmp_input_synapse->next = malloc(sizeof(struct inc_synapse));
                    tmp_input_synapse = tmp_input_synapse->next;
                    tmp_input_synapse->syn_pointer = NULL;
                    tmp_input_synapse->next = NULL;

                    // moving pointer to the next list item
                    tmp_synapse->next = malloc(sizeof(struct synapse));
                    tmp_synapse = tmp_synapse->next;
                    tmp_synapse->post_id = -1;
                    tmp_synapse->next = NULL;

                    // finally writing the created connection down to the file
                    if (pre_syn_type == TYPE_INHIBITORY)
                        fprintf(output_connections,"%d ",-i);
                    else
                        fprintf(output_connections,"%d ",i);

                    if (post_syn_type == TYPE_INHIBITORY)
                        fprintf(output_connections,"%d ",-j);
                    else
                        fprintf(output_connections,"%d ",j);

                    fprintf(output_connections,"%g\n",prob);

                }
                if (TOPOLOGY_LOAD_TYPE == READ_TOPOLOGY_FROM_FILE)
                    if (fscanf(input_connections,"%d %d %lf",&pre_num,&post_num,&tmp_prob) == EOF && i != N-1){
                        printf("ERROR: Saved topology file ended, neuron number mismatch on n=%d\n",i);
                        getchar();
                        exit(1);
                    }

            }


        }

        //writing the amount of outcoming connections per neuron down to the file
       fprintf(output_connections_distribution,"%d %d\n",i,CONN_PER_NEURON_COUNT);
       if ((i+1) % 100 == 0)
        printf("\rConnections done for %d neurons...",i+1);
    }

    //closing files we don't need anymore
    fclose(output_connections);
    fclose(output_connections_distribution);
    fclose(synaptic_parameters_distribution);

    if (TOPOLOGY_LOAD_TYPE == READ_TOPOLOGY_FROM_FILE){
        fclose(input_connections);
    }

    if (SYNAPTIC_DATA_LOAD_TYPE == READ_SYNAPTIC_DATA_FROM_FILE){
        fclose(input_synaptic_parameters_distribution);
    }

    printf("\nNetwork is created successfully.\nReady to start simulation.\n");

    //****************************************Simulation**********************************************************************************************

    printf("Simulating...\n");
    output_info = fopen("info.txt","w");

    //Initializing all necessary variables
    averaging_timer = AVG_TIME;
    synapse_resource_timer = dt;
    N_SP = 0;
    N_ACTIVE_CONN = 0;
    avg_x = 0;
    avg_y = 0;
    avg_z = 0;
    avg_u = 0;
    spike_resource = (N - N_PM)*max_spikes;

    last_burst_time=0;
    burst_flag = 0;
    burst_counter = 0;
    avg_activity = 0;
    burst_detection_timer = 0;

    if (HOMEOSTASIS_status == HOMEOSTASIS_IS_ON)
        M=M_max;
    else M = 1;

    for (i=0;i<N;i++) neurons[i]->I_b_init = neurons[i]->I_b;

    //the main cycle of simulation
    for (simulation_time = 0; simulation_time < SIM_TIME; simulation_time += dt){

        // PASS 1: step for all synapses

        for (i=0;i<N;i++){

            // if homeostasis is turned ON, 0 < M < 1, M = 1 otherwise.
            neurons[i]->I_b = neurons[i]->I_b_init*M;
            // if background current noise is turned ON, the deviation from
            if (BG_CURRENT_NOISE_MODE == 1) neurons[i]->I_b += gauss(0,I_BG_NOISE_SD,-neurons[i]->I_b,100*neurons[i]->I_b);

            //For each incoming synapse of a neuron i we
            //                  a)make a step for timers in a queue of incoming spikes and process them
            //                  b)make a step for synaptic resources
            tmp_input_synapse = neurons[i]->in_conn;
            while (tmp_input_synapse->syn_pointer != NULL){
                stepSynapse(tmp_input_synapse->syn_pointer, dt);

                if (tmp_input_synapse->syn_pointer->timers != NULL){

                    tmp_timer = tmp_input_synapse->syn_pointer->timers;

                    while (tmp_timer->next != NULL){
                        tmp_timer2 = tmp_timer;
                        tmp_timer->time_left -= dt;
                        tmp_timer = tmp_timer->next;
                    }

                    tmp_timer->time_left -= dt;
                    //    Processing a spike
                    if (tmp_timer->time_left <= 0){

                        old_u = tmp_input_synapse->syn_pointer->u;
                        if (neurons[tmp_input_synapse->syn_pointer->pre_id]->type == TYPE_EXCITATORY)
                            tmp_input_synapse->syn_pointer->u = tmp_input_synapse->syn_pointer->w;
                        else tmp_input_synapse->syn_pointer->u += tmp_input_synapse->syn_pointer->U*(1 - old_u);
                        tmp_input_synapse->syn_pointer->y += old_u*tmp_input_synapse->syn_pointer->x;
                        tmp_input_synapse->syn_pointer->x -= old_u*tmp_input_synapse->syn_pointer->x;

                        if (STDP_status != 0)
                            tmp_input_synapse->syn_pointer->w += dw(tmp_input_synapse->syn_pointer->w,neurons[i]->last_spiked_at - (neurons[tmp_input_synapse->syn_pointer->pre_id]->last_spiked_at + tmp_input_synapse->syn_pointer->tau_delay),STDP_status,A_plus,A_minus,tau_plus_corr,tau_minus_corr);

                        if (tmp_input_synapse->syn_pointer->timers->next == NULL)
                            tmp_input_synapse->syn_pointer->timers = NULL;
                        else
                            tmp_timer2->next = NULL;

                        free(tmp_timer);

                    }

                }

                tmp_input_synapse = tmp_input_synapse->next;
            }
        }



        //PASS 2: step for membrane potentials V, detecting and proccessing spikes


        for (i=0;i<N;i++){
            //decreasing refractory state timers
            if (neurons[i]->ref_time_left > 0){
                neurons[i]->ref_time_left -= dt;
                continue;
            }

            //membrane potential step
            switch (NEURON_MODEL){
                case NEURON_MODEL_LEAKY_INTEGRATE_AND_FIRE:
                    neurons[i]->V = stepV_LIF(neurons[i], dt, R_IN, TAU_M);
                    break;
                case NEURON_MODEL_PERFECT_INTEGRATE_AND_FIRE:
                    neurons[i]->V = stepV_PIF(neurons[i],dt,C_M);
                    break;
            }

            if (simulation_time > inh_off_time && inh_off_time != -1 && neurons[i]->type==TYPE_INHIBITORY)
                neurons[i]->V = V_REST;

            if (max_spikes != -1 && neurons[i]->I_b_init < V_TH/R_IN) {
                if (max_spikes - neurons[i]->num_spikes <= 0){
                    neurons[i]->V = neurons[i]->V_rest;
                }
            }

            //check if a neuron emits a spontaneous spike if the stimulation type is corresponding
            if (STIMULATION_TYPE == STIM_TYPE_P_SP_GAUSSIAN  || STIMULATION_TYPE == STIM_TYPE_P_SP_TWO_VALUES){

                prob = (double)(rand())/RAND_MAX;

                if (prob < neurons[i]->p_sp && neurons[i]->ref_time_left <= 0){ //*****Spontaneous spiking
                    neurons[i]->V = V_TH;
                    printf("Stimulating %d!\n",i);

                }
            }


            if (neurons[i]->V >= V_TH) { //In case of exceeding threshold value neuron emits a spike

                N_SP++; //we count this spike to get the activity later
                N_ACTIVE_CONN += neurons[i]->num_of_outcoming_connections;

                neurons[i]->num_spikes++;
                spike_resource--;

                if (max_spikes != -1 && neurons[i]->I_b_init < V_TH/R_IN)
                    if (max_spikes - neurons[i]->num_spikes <= 0){
                        neurons[i]->I_b_init = 0;
                        fprintf(output_lifetimes,"%d %g %g %g\n",i,simulation_time,neurons[i]->x, neurons[i]->y);
                        fflush(output_lifetimes);
                    }

                // we give a message to user


                //we write the spike down to the raster file
                if (neurons[i]->I_b >= V_TH/R_IN){
                    printf("Neuron %d \t1__PM emitted spike at %.1lf ms \t|| Bursts registered: %d\n",i,simulation_time,burst_counter);
                    fprintf(output_raster,"%g %d pm\n", simulation_time, i);
                } else {
                    printf("Neuron %d \t0_reg emitted spike at %.1lf ms \t|| Bursts registered: %d\n",i,simulation_time,burst_counter);
                    fprintf(output_raster,"%g %d\n", simulation_time, i);
                }
                fflush(output_raster);

                neurons[i]->last_spiked_at = simulation_time;

                //We send a spike to all outcoming synapses by adding a new item in a queue of timers of incoming spikes
                tmp_synapse = neurons[i]->out_conn;
                while (tmp_synapse->post_id != -1){

                    tmp_timer = tmp_synapse->timers;
                    tmp_synapse->timers = malloc(sizeof(struct delay_timer));
                    tmp_synapse->timers->time_left = tmp_synapse->tau_delay;
                    tmp_synapse->timers->next = tmp_timer;

                    tmp_synapse = tmp_synapse->next;
                }

                //Applying STDP rules, if turned ON
                if (STDP_status != 0){
                    tmp_input_synapse = neurons[i]->in_conn;
                    while (tmp_input_synapse->syn_pointer != NULL){

                        tmp_input_synapse->syn_pointer->w += dw(tmp_input_synapse->syn_pointer->w,neurons[i]->last_spiked_at - (neurons[tmp_input_synapse->syn_pointer->pre_id]->last_spiked_at + tmp_input_synapse->syn_pointer->tau_delay),STDP_status,A_plus,A_minus,tau_plus_corr,tau_minus_corr);
                        tmp_input_synapse = tmp_input_synapse->next;
                    }
                }

                // setting refractory period
                if (neurons[i]->type == TYPE_EXCITATORY) neurons[i]->ref_time_left = TAU_REF_EXC;
                else neurons[i]->ref_time_left = TAU_REF_INH;

                //setting resetting potential
                neurons[i]->V = neurons[i]->V_reset;

                if (HOMEOSTASIS_status == HOMEOSTASIS_IS_ON)
                    M -= M*b; //Homeostasis mechanism

            }

        }

    //****Homeostasis step
        if (HOMEOSTASIS_status == HOMEOSTASIS_IS_ON)
            M += stepM(M,M_max,TAU_M,dt);

//********Last part of an iteration is decreasing all existing timers and processing events if some of them hit 0

        //****Getting activity and counting bursts
        averaging_timer -= dt;
        if (burst_detection_timer > 0) burst_detection_timer -= dt;
        if (averaging_timer <= 0){

            //printf("\nSimulation time: %.0f ms, Bursts registered: %d\n",simulation_time,burst_counter);

            activity = ((double)N_SP)/N;
            fprintf(output_file_spiking_activity,"%g %g\n", simulation_time, activity);
            fflush(output_file_spiking_activity);

            if (max_spikes != -1){
                fprintf(output_spike_resource,"%g %g\n",simulation_time,spike_resource/(double)((N-N_PM)*max_spikes));
                fflush(output_spike_resource);
            }

            fprintf(output_amount_of_active_connections,"%g %g\n", simulation_time, (double)N_ACTIVE_CONN/(double)N_con);
            fflush(output_amount_of_active_connections);

            avg_activity += activity;
                if (activity > burst_threshold && burst_flag == 0 && burst_detection_timer <= 0){
                    burst_flag = 1;
                    burst_detection_timer = 100;                    
                    fprintf(output_burst_times,"%g\n",simulation_time);
                    fflush(output_burst_times);

                }
                if (activity < burst_threshold && burst_flag == 1){
                    burst_flag = 0;
                    if (last_burst_time != 0){
                        fprintf(output_IBI,"%g\n",simulation_time - last_burst_time);
                        fflush(output_IBI);
                    }
                    burst_counter++;
                    last_burst_time = simulation_time;
                }
            averaging_timer = AVG_TIME;
            N_SP = 0;
            N_ACTIVE_CONN = 0;

            //****Disabling stimulation (e.g. background currents or spontaneous spiking)
            if (stim_disable_protocol != 0){
                        if (stim_disable_protocol == 1 && burst_counter >= pb_number - 1 && ((burst_detection_timer <= 0) || (burst_flag == 1 && burst_detection_timer > 0))&& activity >= activity_threshold ){
                            for (i=0;i<N;i++){
                                neurons[i]->I_b_init = 0;
                                neurons[i]->p_sp = 0;

                            }
                            stim_disable_protocol = 0;
                            printf("\n\n\n********\n\nSTIMULATION DISABLING\n\n********\n\n");
                            fprintf(output_info,"EVENTS: Stimulation was disabled at %g ms\n", simulation_time);
                        }
                        if (stim_disable_protocol == 2 && simulation_time >= stim_off_time){
                            for (i=0;i<N;i++){
                                neurons[i]->I_b_init = 0;
                                neurons[i]->p_sp = 0;
                            }
                            stim_disable_protocol = 0;
                            printf("\n\n\n********\n\nSTIMULATION DISABLING\n\n********\n\n");
                            fprintf(output_info,"EVENTS: Stimulation was disabled at %g ms\n", simulation_time);
                        }
                        if (stim_disable_protocol == 3 && simulation_time >= stim_off_time && burst_counter >= pb_number - 1 && activity >= activity_threshold){
                            for (i=0;i<N;i++){
                                neurons[i]->I_b_init = 0;
                                neurons[i]->p_sp = 0;
                            }
                            stim_disable_protocol = 0;
                            printf("\n\n\n********\n\nSTIMULATION DISABLING\n\n********\n\n");
                            fprintf(output_info,"EVENTS: Stimulation was disabled at %g ms\n", simulation_time);
                        }
                    }
        }

        //****Exporting synaptic resource values
        synapse_resource_timer -= dt;
        if(synapse_resource_timer <= 0){

            for(i=0;i<N;i++){

                tmp_synapse = neurons[i]->out_conn;

                while (tmp_synapse->post_id != -1){
                    avg_x += tmp_synapse->x;
                    avg_y += tmp_synapse->y;
                    avg_z += tmp_synapse->z;
                    avg_u += tmp_synapse->u;
                    tmp_synapse = tmp_synapse->next;
                }

            }

            fprintf(output_x,"%g %g\n", simulation_time, avg_x/((double)N_con));
            fprintf(output_y,"%g %g\n", simulation_time, avg_y/((double)N_con));
            fprintf(output_z,"%g %g\n", simulation_time, avg_z/((double)N_con));
            fprintf(output_u,"%g %g\n", simulation_time, avg_u/((double)N_con));
            fprintf(output_M,"%g %g\n", simulation_time, M);
            fflush(output_M);

            synapse_resource_timer = SYN_RESOURCES_OUTPUT_PERIOD;

            avg_x = 0;
            avg_y = 0;
            avg_z = 0;
            avg_u = 0;

        }

        //****Exporting connection weights periodically if STDP is ON
        if (STDP_status != STDP_IS_OFF){
            w_output_timer -= dt;
            if (w_output_timer <= 0){

                if (layer_type == LAYER_TYPE_BARBELL){ //writing weights for "right" and "left" connections separately in case of BARBELL topology
                    tmp_str[0]='\0';
                    strcat(tmp_str,"stdp_w_dynamics/R/r_w");
                    strcat(tmp_str,my_itoa(W_OUTPUT_COUNTER,tmp_str_2));
                    strcat(tmp_str,".txt");
                    output_w = fopen(tmp_str,"w");
                    for (i=0;i<N;i++){
                        tmp_synapse = neurons[i]->out_conn;
                        if (neurons[i]->type == TYPE_EXCITATORY)
                            while (tmp_synapse->post_id != -1){
                                if (neurons[tmp_synapse->post_id]->x + neurons[tmp_synapse->post_id]->y > neurons[tmp_synapse->pre_id]->x + neurons[tmp_synapse->pre_id]->y)
                                    fprintf(output_w,"%d %d %g\n",tmp_synapse->pre_id, tmp_synapse->post_id, tmp_synapse->w);
                                tmp_synapse = tmp_synapse->next;
                            }
                    }
                    fclose(output_w);

                    tmp_str[0]='\0';
                    strcat(tmp_str,"stdp_w_dynamics/L/l_w");
                    strcat(tmp_str,my_itoa(W_OUTPUT_COUNTER,tmp_str_2));
                    strcat(tmp_str,".txt");
                    output_w = fopen(tmp_str,"w");
                    for (i=0;i<N;i++){
                        tmp_synapse = neurons[i]->out_conn;
                        if (neurons[i]->type == TYPE_EXCITATORY)
                            while (tmp_synapse->post_id != -1){
                                if (neurons[tmp_synapse->post_id]->x + neurons[tmp_synapse->post_id]->y < neurons[tmp_synapse->pre_id]->x + neurons[tmp_synapse->pre_id]->y)
                                    fprintf(output_w,"%d %d %g\n",tmp_synapse->pre_id, tmp_synapse->post_id, tmp_synapse->w);
                                tmp_synapse = tmp_synapse->next;
                            }
                    }
                    fclose(output_w);

                } else { //proceeding with usual synaptic weight output

                    tmp_str[0]='\0';
                    strcat(tmp_str,"stdp_w_dynamics/w");
                    strcat(tmp_str,my_itoa(W_OUTPUT_COUNTER,tmp_str_2));
                    strcat(tmp_str,".txt");
                    output_w = fopen(tmp_str,"w");
                    for (i=0;i<N;i++){
                        tmp_synapse = neurons[i]->out_conn;
                        if (neurons[i]->type == TYPE_EXCITATORY)
                        while (tmp_synapse->post_id != -1){
                            fprintf(output_w,"%d %d %g\n",tmp_synapse->pre_id, tmp_synapse->post_id, tmp_synapse->w);
                            tmp_synapse = tmp_synapse->next;
                        }
                    }
                    fclose(output_w);
                }
                w_output_timer = W_OUTPUT_PERIOD;
                W_OUTPUT_COUNTER++;
            }
        }




    }

    //*******************After Simulation********************
    //Writing simulation parameters to the info.txt file


    fprintf(output_info,"\n>>>> GENERAL <<<<\n");
    fprintf(output_info,"N = %d\n",N);
    fprintf(output_info,"INHIBITORY neurons fraction: %g\n",INH_NEURONS_FRACTION);

    fprintf(output_info,"\n\n>>>> NEURON MODEL <<<<\n");
    if (NEURON_MODEL == NEURON_MODEL_PERFECT_INTEGRATE_AND_FIRE){
        fprintf(output_info,"Neuron Model: Perfect Integrate-and-Fire (PIF)\n");
        fprintf(output_info,"C_m = %g pF\n",C_M);
    } else {
        fprintf(output_info,"Neuron Model: Leaky Integrate-and-Fire (LIF)\n");
        fprintf(output_info,"TAU_M = %g ms\n",TAU_M);
        fprintf(output_info,"R_in = %g GOhm\n",R_IN);
    }

    fprintf(output_info,"V_REST = %g mV\n",V_REST);
    fprintf(output_info,"V_TH = %g mV\n",V_TH);
    fprintf(output_info,"TAU_REF_EXC = %g ms\n",TAU_REF_EXC);
    fprintf(output_info,"TAU_REF_INH = %g ms\n",TAU_REF_INH);

    fprintf(output_info,"\n\n>>>> STIMULATION <<<<");

    if (STIMULATION_TYPE == STIM_TYPE_I_BG_GAUSSIAN || STIMULATION_TYPE == STIM_TYPE_I_BG_TWO_VALUES){
        fprintf(output_info,"\nStimulation type: Background currents\n");
        if (STIMULATION_DATA_LOAD_TYPE == 1){
            fprintf(output_info,"\nBackground currents will be read from file\n");
        } else {
            fprintf(output_info,"Background currents will be generated with following properties\n");
            if (STIMULATION_TYPE == STIM_TYPE_I_BG_GAUSSIAN){
                fprintf(output_info,"I_bg_mean = %g pA\n",I_BG_MEAN);
                fprintf(output_info,"I_bg_sd = %g pA\n",I_BG_SD);
                fprintf(output_info,"I_bg_min = %g pA\n",I_BG_MIN);
                fprintf(output_info,"I_bg_max = %g pA\n",I_BG_MAX);
            } else {
                fprintf(output_info,"I_bg_1 = %g pA\n",I_BG_1);
                fprintf(output_info,"I_bg_2 = %g pA\n", I_BG_2);
                fprintf(output_info,"Fraction_of_neurons_with_I_bg_1 = %g\n", FRACTION_OF_NEURONS_WITH_I_BG_1);
            }
        }
    } else {
        fprintf(output_info,"Stimulation type: Spontaneous spiking\n");
        if (STIMULATION_DATA_LOAD_TYPE == 1){
            fprintf(output_info,"\nP_sp probabilities will be read from file\n");
        } else {
            if (STIMULATION_TYPE == STIM_TYPE_P_SP_GAUSSIAN){
                fprintf(output_info,"\nP_sp probabilities will be generated with following properties\n");
                fprintf(output_info,"P_SP_mean = %g\n",P_SP_MEAN);
                fprintf(output_info,"P_SP_sd = %g\n",P_SP_SD);
                fprintf(output_info,"P_SP_min = %g\n",P_SP_MIN);
                fprintf(output_info,"P_SP_max = %g\n",P_SP_MAX);
            } else {
                fprintf(output_info,"p_sp_1 = %g\n",P_SP_1);
                fprintf(output_info,"p_sp_2 = %g\n", P_SP_2);
                fprintf(output_info,"Fraction_of_neurons_with_p_sp_1 = %g\n", FRACTION_OF_NEURONS_WITH_P_SP_1);
            }
        }
    }

    fprintf(output_info,"\n\n>>>> TOPOLOGY <<<<\n");
    fprintf(output_info,"Layer type: ");
    switch(layer_type){
    case LAYER_TYPE_BINOMIAL:
        fprintf(output_info,"BINOMIAL\n");
        break;
    case LAYER_TYPE_UNIFORM:
        fprintf(output_info,"UNIFORM\n");
        break;
    case LAYER_TYPE_SQUARE_LATTICE:
        fprintf(output_info,"SQUARE_LATTICE\n");
        break;
    case LAYER_TYPE_BELL:
        fprintf(output_info,"BELL\n");
        break;
    case LAYER_TYPE_RAMP:
        fprintf(output_info,"RAMP\n");
        break;
    case LAYER_TYPE_DOUBLE_RAMP:
        fprintf(output_info,"DOUBLE RAMP\n");
        break;
    case LAYER_TYPE_BARBELL:
        fprintf(output_info,"BARBELL\n");
        break;
    case LAYER_TYPE_STRIPED:
        fprintf(output_info,"STRIPED\n");
        break;
    }

    if (layer_type == LAYER_TYPE_BINOMIAL)
        fprintf(output_info,"P_con = %g\n",P_CON);
    else{
        fprintf(output_info,"lambda = %g\n",lambda);
        fprintf(output_info,"max_conn_length = %g L\n",max_conn_length);
        fprintf(output_info,"spike_speed = %g L/ms\n",SPIKE_SPEED);
        fprintf(output_info,"tau_delay = %g ms\n",tau_delay);
    }

    if (TOPOLOGY_LOAD_TYPE == READ_TOPOLOGY_FROM_FILE || SYNAPTIC_DATA_LOAD_TYPE == READ_SYNAPTIC_DATA_FROM_FILE ||
            COORDINATES_LOAD_TYPE == READ_COORDINATES_FROM_FILE || STIMULATION_DATA_LOAD_TYPE == READ_STIMULATION_DATA_FROM_FILE){
        fprintf(output_info,"\n\n>>>> USING SAVED DATA <<<<");

        if (TOPOLOGY_LOAD_TYPE == READ_TOPOLOGY_FROM_FILE){
            fprintf(output_info,"\nTOPOLOGY is being READ FROM FILE\n");
        } else {
                              // in case we don't read topology input file we generate it
            printf("\nTOPOLOGY will be GENERATED for this run\n");
        }

        if (SYNAPTIC_DATA_LOAD_TYPE == READ_SYNAPTIC_DATA_FROM_FILE){

            if (TOPOLOGY_LOAD_TYPE != READ_TOPOLOGY_FROM_FILE){
                fprintf(output_info,"\nERROR: Synaptic parameters cannot be read from file if topology is not read from file!\n");
                getchar();
                exit(1);
            }

            fprintf(output_info,"SYNAPTIC PARAMETERS distribution is being READ FROM FILE\n");

        } else {
                              // in case we don't read synaptic parameters from file we generate them
            fprintf(output_info,"SYNAPTIC PARAMETERS will be GENERATED for this run\n");
        }
        if (COORDINATES_LOAD_TYPE == READ_COORDINATES_FROM_FILE){
            fprintf(output_info,"COORDINATES are being READ FROM FILE expecting ");
        } else { fprintf(output_info,"COORDINATES will be GENERATED with ");}
            switch(layer_type){

         case LAYER_TYPE_BINOMIAL:
                fprintf(output_info,"layer type BINOMIAL\n");
                break;
            case LAYER_TYPE_UNIFORM:
                fprintf(output_info,"layer type UNIFORM\n");
                break;
            case LAYER_TYPE_SQUARE_LATTICE:
                fprintf(output_info,"layer type SQUARE_LATTICE\n");
                break;
            case LAYER_TYPE_BELL:
                fprintf(output_info,"layer type BELL\n");
                break;
            case LAYER_TYPE_RAMP:
                fprintf(output_info,"layer type RAMP\n");
                break;
            case LAYER_TYPE_DOUBLE_RAMP:
                fprintf(output_info,"layer type DOUBLE RAMP\n");
                break;
            case LAYER_TYPE_BARBELL:
                fprintf(output_info,"layer type BARBELL\n");
                break;
            case LAYER_TYPE_STRIPED:
                fprintf(output_info,"layer type STRIPED\n");
                break;
        }
        if (STIMULATION_DATA_LOAD_TYPE == READ_STIMULATION_DATA_FROM_FILE){
            fprintf(output_info,"\nSTIMULATION DATA is being READ FROM FILE\n");
        } else {
            fprintf(output_info,"\nSTIMULATION DATA will be GENERATED for this run\n");
        }
    }
    fprintf(output_info,"\n\n>>>> ADDITIONAL MECHANISMS <<<<");

    if (BG_CURRENT_NOISE_MODE == 0)
        fprintf(output_info,"\nBG current noise is OFF\n");
    else{
        fprintf(output_info,"\nBG current noise is ON\n");
        fprintf(output_info,"I_bg_noise_sd = %g pA\n",I_BG_NOISE_SD);
    }
    switch(STDP_status){
    case STDP_IS_OFF:
        fprintf(output_info,"STDP status is OFF\n");
        break;
    case STDP_IS_ADDITIVE:
        fprintf(output_info,"STDP status is ADDITIVE\n");
        break;
    case STDP_IS_MULTIPLICATIVE:
        fprintf(output_info,"STDP status is MULTIPLICATIVE\n");
        break;
    }
    if (HOMEOSTASIS_status == HOMEOSTASIS_IS_ON)
        fprintf(output_info,"HOMEOSTASIS is ON\n");
    else
        fprintf(output_info,"HOMEOSTASIS is OFF\n");

    if (stim_disable_protocol != 0)
            fprintf(output_info,"STIMULATION DISABLING protocol is ON\n");

    if (inh_off_time != -1)
            fprintf(output_info,"\nInhibitory neurons will be forced to V_REST on %g ms\n",inh_off_time);

    fprintf(output_info,"\n\n>>>> SIMULATION PARAMETERS <<<<");
    fprintf(output_info,"\nBurst detection theshold: %g",burst_threshold);
    fprintf(output_info,"\nTime step dt = %g ms",dt);
        fprintf(output_info,"\nAveraging time: %g ms\n",AVG_TIME);
        fprintf(output_info,"Simulation time: %d ms\n",SIM_TIME);
        fprintf(output_info,"Number of bursts in simulation: %d\n",burst_counter);
        fclose(output_info);

//*************Releasing memory************************************************************************
    printf("Releasing memory...\n");

    for (i=0;i<N;i++){
        tmp_synapse = neurons[i]->out_conn;
        tmp_synapse2 = tmp_synapse;
        while (tmp_synapse->post_id != -1){
            tmp_synapse = tmp_synapse->next;
            free(tmp_synapse2);
            tmp_synapse2 = tmp_synapse;
        }
        free(tmp_synapse2);

        tmp_input_synapse = neurons[i]->in_conn;
        tmp_input_synapse2 = tmp_input_synapse;
        while (tmp_input_synapse->syn_pointer != NULL){
            tmp_input_synapse = tmp_input_synapse->next;
            free(tmp_input_synapse2);
            tmp_input_synapse2 = tmp_input_synapse;
        }
        free(tmp_input_synapse2);


        free(neurons[i]);
    }

    free(neurons);
    printf("Done\n");
    fflush(stdin);
    getchar();

    //***************** Closing rest of the files *************************************************************************
    fclose(output_raster);

    fclose(output_burst_times);
    fclose(output_IBI);


    fclose(output_file_spiking_activity);
    fclose(output_amount_of_active_connections);

    fclose(output_x);
    fclose(output_y);
    fclose(output_z);
    fclose(output_u);

    if (max_spikes != -1){
        fclose(output_lifetimes);
        fclose(output_spike_resource);
    }

    fclose(output_M);

    return 0;
}

//****Input file-reading function
void readfile (FILE *file, struct parameter value_array[PARAMETER_N]){

    int i;
    char str[BUFFER_SIZE];
    char *value;

    while (fgets(str, BUFFER_SIZE, file) != NULL){

        value = str;
        while (*value != '=' && *value != '\0') value++;

        if (*value == '\0') continue;

        *value = '\0';
        value++;

        for (i=0; i<PARAMETER_N; i++){
            if (strcmp(str, value_array[i].name) == 0){

                if (strcmp(str,"N") == 0 || strcmp(str,"SIM_TIME") == 0 || strcmp(str,"pb_number") == 0 || strcmp(str,"neuron_model") == 0 || strcmp(str,"stim_disable_protocol") == 0 || strcmp(str,"Homeostasis_status") == 0  || strcmp(str,"I_bg_noise_mode") == 0 || strcmp(str,"STDP_status") == 0 || strcmp(str,"force_binomial_topology") == 0 || strcmp(str,"use_saved_coordinates") == 0 || strcmp(str,"spatial_layer_type") == 0 || strcmp(str,"max_spikes_per_neuron") == 0 || strcmp(str,"Stimulation_type") == 0 || strcmp(str,"use_saved_topology") == 0 || strcmp(str,"use_saved_stimulation_data") == 0 || strcmp(str,"use_saved_synaptic_parameters") == 0) //these parameters are integer, so they need special care
                    *(int *)value_array[i].value = atoi(value);

                else *(double *)value_array[i].value = atof(value);


            }
        }
    }
}

