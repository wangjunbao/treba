/**************************************************************************/
/*   treba - probabilistic finite-state automaton training and decoding   */
/*   Copyright © 2012 Mans Hulden                                         */

/*   This file is part of treba.                                          */

/*   Treba is free software: you can redistribute it and/or modify        */
/*   it under the terms of the GNU General Public License version 2 as    */
/*   published by the Free Software Foundation.                           */

/*   Treba is distributed in the hope that it will be useful,             */
/*   but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/*   GNU General Public License for more details.                         */

/*   You should have received a copy of the GNU General Public License    */
/*   along with treba.  If not, see <http://www.gnu.org/licenses/>.       */
/**************************************************************************/

#define LIKELIHOOD_FORWARD      1
#define LIKELIHOOD_BACKWARD     2
#define LIKELIHOOD_VITERBI      3
#define DECODE_FORWARD          4
#define DECODE_BACKWARD         5
#define DECODE_VITERBI          6
#define DECODE_FORWARD_PROB     7
#define DECODE_BACKWARD_PROB    8
#define DECODE_VITERBI_PROB     9
#define TRAIN_BAUM_WELCH        10
#define TRAIN_DA_BAUM_WELCH     11
#define TRAIN_VITERBI           12
#define TRAIN_VITERBI_BW        13
#define TRAIN_VARIATIONAL_BAYES 14
#define TRAIN_GIBBS_SAMPLING    15
#define GENERATE_WORDS          16

#define GENERATE_NONDETERMINISTIC 1
#define GENERATE_DETERMINISTIC    2
#define GENERATE_UNIFORM          3
#define GENERATE_BAKIS            4

#define SMRONE_REAL  1

#define FORMAT_REAL     0
#define FORMAT_LOG10    1
#define FORMAT_LOG2     2
#define FORMAT_LN       3
#define FORMAT_NLOG10   4
#define FORMAT_NLOG2    5
#define FORMAT_NLN      6

#ifdef MATH_FLOAT
 #define LOG(X)        (log2f((X)))
 #define EXP(X)        (exp2f((X)))
 #define SMRZERO_LOG  -FLT_MAX
 #define LOGZERO       FLT_MAX
 #define ABS(X)        (fabsf((X)))
 typedef float PROB;
#else
 #define LOG(X)        (log2((X)))
 #define EXP(X)        (exp2((X)))
 #define SMRZERO_LOG  -DBL_MAX
 #define LOGZERO       DBL_MAX
 #define ABS(X)        (fabs((X)))
 typedef double PROB;
#endif /* MATH_FLOAT */


/* Auxiliary macros to access trellis, FSMs, and FSM counts */
#define FINALPROB(FSM, STATE) ((FSM)->final_table + (STATE))
#define TRELLIS_CELL(STATE, TIME) ((trellis) + (fsm->num_states) * (TIME) + (STATE))
#define TRANSITION(FSM, SOURCE_STATE, SYMBOL, TARGET_STATE) ((FSM)->state_table + (FSM)->num_states * (FSM)->alphabet_size * (SOURCE_STATE) + (SYMBOL) * (FSM)->num_states + (TARGET_STATE))
#define FSM_COUNTS(FSMC, SOURCE_STATE, SYMBOL, TARGET_STATE) ((FSMC) + (fsm->num_states * fsm->alphabet_size * (SOURCE_STATE) + (SYMBOL) * fsm->num_states + (TARGET_STATE)))

PROB smrzero = SMRZERO_LOG;
PROB smrone  = 0;

PROB *fsm_counts, *fsm_totalcounts, *fsm_finalcounts;
_Bool *fsm_counts_spin;

pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
PROB g_loglikelihood = 0;

struct thread_args {
    struct trellis *trellis;
    struct observations **obsarray;
    int minobs;
    int maxobs;
    struct wfsa *fsm;
    PROB beta;
};

struct observations *g_obsarray;

struct wfsa {
    int num_states;
    int alphabet_size;
    PROB *state_table;
    PROB *final_table;
};

struct observations {
    int size;
    int *data;
    int occurrences;
    struct observations *next;
};

struct trellis {
    PROB fp;
    PROB bp;
    int backstate;
};

/*******************/
/* Default options */
/*******************/

int g_numcpus = 1;
int g_maxiterations = 100000;
PROB g_maxdelta = 0.1;
int g_input_format = FORMAT_REAL;
int g_output_format = FORMAT_REAL;
int g_viterbi_pseudocount = 1;
int g_random_restarts = 0;
int g_random_restart_iterations = 3;
int g_generate_type = 0;
int g_gen_max_length = 100000;
int g_num_states = 5;
int g_alphabet_size = -1;
int g_initialize_uniform = 0;
int g_train_da_bw = 0;
int g_generate_words = 0;
/* Thread variables */
int g_num_threads = 1;
/* Deterministic annealing default parameters */
PROB g_betamin = 0.02;
PROB g_betamax = 1;
PROB g_alpha = 1.01;
/* Variational Bayes Dirichlet parameter */
PROB g_vb_alpha = 0.02;
/* Collapsed Gibbs sampling concentration, lag, burnin parameters */
PROB g_gibbs_beta = 0.02;
int g_gibbs_lag = 100;
int g_gibbs_burnin = 10000;
/* Flag whether to adjust counts in BW by VB strategy */
int g_bw_vb = 0;

struct wfsa *g_lastwfsa = NULL; /* Stores pointer to last WFSA to spit out in case of SIGINT */

static char *versionstring = "treba v0.1";
static char *usagestring = 
"treba [options] observationfile\n";
static char *helpstring = 
"\n"
"Train FSMs, calculate probabilities and best paths of sequences.\n\n"

"Main options:\n"
"-T bw|dabw|vit|vitbw\tTrain FSM with Baum-Welch, B-W + deterministic annealing, Viterbi training, Viterbi+B-W.\n"
"-D f|b|vit(,p)\t\tDecode (find best path through automaton) with forward, backward, or Viterbi.\n"
"-L f|vit\t\tCalculate probability of sequences; forward probability or best path (Viterbi).\n"
"-i format\t\tSet input format (real, log10, ln, log2, nlog10, nln, nlog2). Default is real.\n"
"-o format\t\tSet output format (real, log10, ln, log2, nlog10, nln, nlog2). Default is real.\n"
"-f fsmfile\t\tSpecify FSM file.\n"
"-g b|d|n#states(,#syms)\tSpecify type of initial random FSM. b=Bakis,d=deterministic,n=ergodic.\n"
"-u\t\t\tSet uniform probabilities on initially generated FSM.\n\n"
"Training options:\n"
"-d max-delta\t\tMaximum change in log likelihood for convergence. Default is 0.1.\n"
"-x maxiterations\tMaximum number of iterations for training. Default is 100000.\n"
"-p pseudocounts\t\tPseudocounts to use for Viterbi training. Default is 1.\n"
"-r restarts,iterations-per-restart\tNumber of restarts and iterations per restart before beginning B-W.\n"
"-a betamin,betamax,alpha\tParameters for deterministic annealing.\n"
"-t num-threads\t\tNumber of threads to launch in Baum-Welch training.\n";

/* Input and output conversion */
PROB output_convert(PROB x);
PROB input_convert(PROB x);

void interrupt_sigproc(void);

inline void spinlock_lock(_Bool *ptr);
inline void spinlock_unlock(_Bool *ptr);

/* Functions to handle file and text input */
char *file_to_mem(char *name);
int char_in_array(char c, char *array);
int line_count_elements(char **ptr);
char *line_to_int_array(char *ptr, int **line, int *size);

/* WFSA functions */
struct wfsa *wfsa_read_file(char *filename);
void wfsa_print(struct wfsa *fsm);
void wfsa_randomize_deterministic(struct wfsa *fsm, int uniform);
void wfsa_randomize_nondeterministic(struct wfsa *fsm, int bakis, int uniform);
struct wfsa *wfsa_init(int num_states, int alphabet_size);
struct wfsa *wfsa_copy(struct wfsa *fsm);
void wfsa_destroy(struct wfsa *fsm);
void wfsa_to_log2(struct wfsa *fsm);
PROB wfsa_sum_prob(struct wfsa *fsm, int state);
int wfsa_random_transition(struct wfsa *fsm, int state, int *symbol, PROB *prob);

/* Generation functions */
void generate_words(struct wfsa *fsm, int numwords);

/* Observation file/array functions */
int obssortcmp(struct observations **a, struct observations **b);
int observations_alphabet_size(struct observations *ohead);
struct observations **observations_to_array(struct observations *ohead, int *numobs);
struct observations *observations_uniq(struct observations *ohead);
struct observations *observations_sort(struct observations *ohead);
void observations_destroy(struct observations *ohead);
struct observations *observations_read(char *filename);

/* Trellis functions */
PROB trellis_backward(struct trellis *trellis, int *obs, int length, struct wfsa *fsm);
PROB trellis_viterbi(struct trellis *trellis, int *obs, int length, struct wfsa *fsm);
PROB trellis_forward(struct trellis *trellis, int *obs, int length, struct wfsa *fsm);
struct trellis *trellis_init(struct observations *o, struct wfsa *fsm);
void trellis_print(struct trellis *trellis, struct wfsa *fsm, int obs_len);

/* Trellis path printing functions */
void forward_print_path(struct trellis *trellis, struct wfsa *fsm, int obs_len);
void backward_print_path(struct trellis *trellis, struct wfsa *fsm, int obs_len);
void viterbi_print_path(struct trellis *trellis, struct wfsa *fsm, int obs_len);

/* Main decoding and likelihood calculations */
void viterbi(struct wfsa *fsm, struct observations *o, int algorithm);
void forward(struct wfsa *fsm, struct observations *o, int algorithm);
void backward(struct wfsa *fsm, struct observations *o, int algorithm);

/* Main training functions */
PROB train_viterbi(struct wfsa *fsm, struct observations *o, int maxiterations, PROB maxdelta);
PROB train_baum_welch(struct wfsa *fsm, struct observations *o, int maxiterations, PROB maxdelta, int vb);
PROB train_bw(struct wfsa *fsm, struct observations *o, int maxiterations, PROB maxdelta);
PROB train_viterbi_bw(struct wfsa *fsm, struct observations *o);
void *trellis_fill_bw(void *threadargs);

int main(int argc, char **argv);
