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

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <pthread.h>
#include <signal.h>

#include "treba.h"
#include "fastlogexp.h"

#ifdef _WIN32
#include <windows.h>
#define srandom srand
#define random rand
#endif

PROB myexp(PROB x) {
    return(x <= smrzero ? 0 : EXP(x));

}
PROB output_convert(PROB x) {
    /* Internal format is log2 */
    switch (g_output_format) {
    case FORMAT_LOG2:   return(x);
    case FORMAT_LOG10:  return(0.30102999566398119521 * x);
    case FORMAT_LN:     return(0.69314718055994530942 * x);
    case FORMAT_REAL:   return(x <= smrzero ? 0 : EXP(x));
    case FORMAT_NLOG2:  return(x == 0 ? 0 : -x);
    case FORMAT_NLOG10: return(x == 0 ? 0 : 0.30102999566398119521 * -x);
    case FORMAT_NLN:    return(x == 0 ? 0 : 0.69314718055994530942 * -x);
    }
    return(x);
}

PROB input_convert(PROB x) {
    /* Internal format is log2 */
    switch (g_input_format) {
    case FORMAT_LOG2:   return(x);
    case FORMAT_LOG10:  return(3.32192809488736234787 * x);
    case FORMAT_LN:     return(1.44269504088896340736 * x);
    case FORMAT_REAL:   return(x == 0 ? smrzero : log2(x));
    case FORMAT_NLOG2:  return(-x);
    case FORMAT_NLOG10: return(x == 0 ? 0 : 3.32192809488736234787 * -x);
    case FORMAT_NLN:    return(x == 0 ? 0 : 1.44269504088896340736 * -x);
    }
    return(x);
}

char *file_to_mem(char *name) {
    FILE    *infile;
    size_t    numbytes;
    char *buffer;
    infile = fopen(name, "r");
    if(infile == NULL) {
        printf("Error opening file '%s'\n",name);
        return NULL;
    }
    fseek(infile, 0L, SEEK_END);
    numbytes = ftell(infile);
    fseek(infile, 0L, SEEK_SET);
    buffer = (char*)malloc((numbytes+1) * sizeof(char));
    if(buffer == NULL) {
        perror("Error reading file");
        return NULL;
    }
    if (fread(buffer, sizeof(char), numbytes, infile) != numbytes) {
        perror("Error reading file");
        return NULL;
    }
    fclose(infile);
    *(buffer+numbytes)='\0';
    return(buffer);
}

static PROB rand_double() {
    return (drand48());
}

static PROB rand_int(int n) {
    long rnd;
    rnd = rand();
    return (int)(rnd % n);
}

int rand_int_range(int from, int to) {
    double ffrom = (double)from;
    double fto = (double)to;
    return from + (int)(fto * random() / (RAND_MAX + ffrom));
}

void wfsa_print(struct wfsa *fsm) {
    int i,j,k;
    PROB thisprob;
    for (i=0; i < fsm->num_states; i++) {
	for (j = 0; j < fsm->alphabet_size; j++) {
	    for (k = 0; k < fsm->num_states; k++) {
		thisprob = *(TRANSITION(fsm,i,j,k));
		if (thisprob > smrzero) {
		    if (g_output_format != FORMAT_REAL || output_convert(thisprob) > 0) {
			printf("%i %i %i %.17g\n", i,k,j, output_convert(thisprob));
		    }
		}
	    }
	}
    }
    for (i=0; i < fsm->num_states; i++) {
	thisprob = *(fsm->final_table + i);
	if (thisprob > smrzero) {
	    if (g_output_format != FORMAT_REAL || output_convert(thisprob) > 0) {
		printf("%i %.17g\n",i,output_convert(thisprob));
	    }
	}
    }
}

void wfsa_randomize_deterministic(struct wfsa *fsm, int uniform) {
    int i, j, k, tablesize;
    PROB randnum, randsum, *tableptr;
    for (i=0; i < fsm->num_states; i++) {
        randsum = 0;
        tablesize = fsm->num_states * fsm->alphabet_size;
        for (j = 0; j < fsm->alphabet_size; j++) {
            k = rand_int(fsm->num_states);
            tableptr = TRANSITION(fsm,i,j,k);
            randnum = (uniform == 0) ? rand_double() : 1;
            *tableptr = randnum;
            randsum += randnum;
        }
        randnum = (uniform == 0) ? rand_double() : 1;
        *(fsm->final_table + i) = randnum;
        randsum += randnum;
        for (tableptr = TRANSITION(fsm,i,0,0), j = 0; j < tablesize; j++) {
            *(tableptr+j) /= randsum;
        }
        *(fsm->final_table + i) /= randsum;
    }
}

void wfsa_randomize_nondeterministic(struct wfsa *fsm, int bakis, int uniform) {
    int i, j, k, tablesize;
    PROB randnum, randsum, *tableptr;
    for (i = 0; i < fsm->num_states; i++) {
	randsum = 0;
	tablesize = fsm->num_states * fsm->alphabet_size;
	for (j = (bakis == 1 ? i : 0); j < fsm->num_states; j++) {
	    for (k = 0 ; k < fsm->alphabet_size; k++) {
		randnum = (uniform == 0) ? rand_double() : 1;
		*TRANSITION(fsm,i,k,j) = randnum;
		randsum += randnum;
	    }
	}
	randnum = (uniform == 0) ? rand_double() : 1;
	*(fsm->final_table + i) = randnum;
	randsum += randnum;
	for (tableptr = TRANSITION(fsm,i,0,0), j = 0; j < tablesize; j++) {
	    *(tableptr+j) /= randsum;
	}
	*(fsm->final_table + i) /= randsum;
    }
}

struct wfsa *wfsa_init(int num_states, int alphabet_size) {
    struct wfsa *fsm;
    fsm = malloc(sizeof(struct wfsa));
    fsm->num_states = num_states;
    fsm->alphabet_size = alphabet_size;
    fsm->state_table = calloc(num_states * num_states * alphabet_size, sizeof(PROB));
    fsm->final_table = calloc(num_states, sizeof(PROB));
    return(fsm);
}

struct wfsa *wfsa_copy(struct wfsa *fsm) {
    struct wfsa *newfsm;
    newfsm = malloc(sizeof(struct wfsa));
    newfsm->num_states = fsm->num_states;
    newfsm->alphabet_size = fsm->alphabet_size;
    newfsm->state_table = malloc(fsm->num_states * fsm->num_states * fsm->alphabet_size * sizeof(PROB));
    newfsm->final_table = malloc(fsm->num_states * sizeof(PROB));
    memcpy(newfsm->state_table, fsm->state_table, fsm->num_states * fsm->num_states * fsm->alphabet_size * sizeof(PROB));
    memcpy(newfsm->final_table, fsm->final_table, fsm->num_states * sizeof(PROB));
    return(newfsm);
}

void wfsa_destroy(struct wfsa *fsm) {
    free(fsm->state_table);
    free(fsm->final_table);
    free(fsm);
}

void wfsa_to_log2(struct wfsa *fsm) {
    int i,j,k;
    for (i = 0; i < fsm->num_states; i++) {
	*FINALPROB(fsm,i) = input_convert(*FINALPROB(fsm,i));
    }
    for (i = 0; i < fsm->num_states; i++) {
	for (j = 0; j < fsm->alphabet_size; j++) {
	    for (k = 0; k < fsm->num_states; k++) {
		*TRANSITION(fsm,i,j,k) = input_convert(*TRANSITION(fsm,i,j,k));
	    }
	}
    }
}

int char_in_array(char c, char *array) {
    int i;
    for (i = 0; *(array+i) != '\0'; i++) {
	if (c == *(array+i)) {
	    return 1;
	}
    }
    return 0;
}

int line_count_elements(char **ptr) {
    int i, elements;
    char seps[] = {'0','1','2','3','4','5','6','7','8','9','.','-','e','E','\0'};
    for (i = 0, elements = 0; *(*ptr+i) != '\n' && *(*ptr+i) != '\0'; i++) {
	if (!char_in_array(*(*ptr+i), seps)) {
	    continue;
	}
	if (!char_in_array(*(*ptr+i+1), seps)) {
	    elements++;
	}
    }
    if (*(*ptr+i) == '\0') {
	*ptr = *ptr+i;
    } else {
	*ptr = *ptr+i+1;
    }
    return(elements);
}

struct wfsa *wfsa_read_file(char *filename) {
    char *wfsa_char_data, *w, *lastline;
    int elements, source, target, symbol, finalstate, maxstate, maxsymbol;
    PROB prob;
    struct wfsa *fsm;
    if ((wfsa_char_data = file_to_mem(filename)) == NULL) {
	exit(1);
    }
    /* Figure out alphabet size and number of states */
    for (w = wfsa_char_data, maxstate = 0, maxsymbol = 0; ; ) {
	lastline = w;
	elements = line_count_elements(&w);
	if (elements == 0) {
	    break;
	}
	switch (elements) {
	case 1:
	    sscanf(lastline, "%i", &finalstate);
	    maxstate = maxstate > finalstate ? maxstate : finalstate;
	    break;
	case 2:
#ifdef MATH_FLOAT
	    sscanf(lastline, "%i %g", &finalstate, &prob);
#else
	    sscanf(lastline, "%i %lg", &finalstate, &prob);
#endif
	    maxstate = maxstate > finalstate ? maxstate : finalstate;
	    break;
	case 3:
	    sscanf(lastline, "%i %i %i", &source, &target, &symbol);
	    maxstate = maxstate > source ? maxstate : source;
	    maxstate = maxstate > target ? maxstate : target;
	    maxsymbol = maxsymbol > symbol ? maxsymbol : symbol;
	    break;
	case 4:
#ifdef MATH_FLOAT
	    sscanf(lastline, "%i %i %i %g", &source, &target, &symbol, &prob);
#else
	    sscanf(lastline, "%i %i %i %lg", &source, &target, &symbol, &prob);
#endif
	    maxstate = maxstate > source ? maxstate : source;
	    maxstate = maxstate > target ? maxstate : target;
	    maxsymbol = maxsymbol > symbol ? maxsymbol : symbol;
	    break;
	default:
	    perror("WFSA file format error");
	    free(wfsa_char_data);
	    exit(1);
	}
    }
    fsm = wfsa_init(maxstate+1, maxsymbol+1);
    for (w = wfsa_char_data, maxstate = 0, maxsymbol = 0; ; ) {
	lastline = w;
	elements = line_count_elements(&w);
	if (elements == 0) {
	    break;
	}
	switch (elements) {
	case 1:
	    sscanf(lastline, "%i", &finalstate);
	    *FINALPROB(fsm, finalstate) = SMRONE_REAL;
	    break;
	case 2:
#ifdef MATH_FLOAT
	    sscanf(lastline, "%i %g", &finalstate, &prob);
#else
	    sscanf(lastline, "%i %lg", &finalstate, &prob);
#endif
	    *FINALPROB(fsm, finalstate) = prob;
	    break;
	case 3:
	    sscanf(lastline, "%i %i %i", &source, &target, &symbol);
	    *TRANSITION(fsm, source, symbol, target) = SMRONE_REAL;
	    break;
	case 4:
#ifdef MATH_FLOAT
	    sscanf(lastline, "%i %i %i %g", &source, &target, &symbol, &prob);
#else
	    sscanf(lastline, "%i %i %i %lg", &source, &target, &symbol, &prob);
#endif
	    *TRANSITION(fsm, source, symbol, target) = prob;
	    break;
	default:
	    perror("WFSA file format error");
	    free(wfsa_char_data);
	    exit(1);
	}
    }
    free(wfsa_char_data);
    return(fsm);
}

char *line_to_int_array(char *ptr, int **line, int *size) {
    /* Reads (destructively) a line of integers (separated by non-integers) and returns a malloced array  */
    /* of numbers (ints) with the line in it + a size count, and also a pointer to the next line.         */
    char *nextline, *startptr;
    int i, j, elements, lastnumber;
    if (*ptr == '\n') {
	*line = NULL;
	*size = 0;
	return ptr+1;
    }
    for (i = 0, elements = 0; ptr[i] != '\n' && ptr[i] != '\0'; i++) {
	/* If ptr+i is a digit and ptr+i+1 is not, we've seen a number */
	if (!isdigit(ptr[i])) {
	    continue;
	}
	if (!isdigit(ptr[i+1])) {
	    elements++;
	}
    }
    if (ptr[i] == '\0') {
	nextline = ptr+i;
    } else {
	nextline = ptr+i+1;
    }
    *size = elements;
    if (!elements) {
	return NULL;
    }
    *line = malloc(sizeof(int) * elements);
    for (i = 0, j = 0, lastnumber = 0; ; i++, j++) {
	/* Find next digit */
	startptr = ptr + i;
	while (!isdigit(ptr[i]) && ptr[i] != '\n' && ptr[i] != '\0') {
	    i++;
	    startptr = ptr+i;
	}
	if (ptr[i] == '\n' || ptr[i] == '\0') {
	    break;
	}
	/* Find number end */
	while (isdigit(ptr[i])) {
	    i++;
	}
	if (ptr[i] == '\n' || ptr[i] == '\0') {
	    lastnumber = 1;
	}
	ptr[i] = '\0';
	*(*line+j) = atoi(startptr);
	if (lastnumber)
	    break;
    }
    return(nextline);
}

int obssortcmp(struct observations **a, struct observations **b) {
    int *one, *two, i;
    one = (*a)->data; two = (*b)->data;
    for (i = 0; i < (*a)->size && i < (*b)->size; i++) { 
	if (*(one+i) < *(two+i)) {
	    return -1;
	}
	if (*(one+i) > *(two+i)) {
	    return 1;
	}
    }
    if ((*a)->size > (*b)->size) {
	return -1;
    }
    if ((*a)->size < (*b)->size) {
	return 1;
    }
    return 0;
}

int observations_alphabet_size(struct observations *ohead) {
    int i, maxsigma;
    for (maxsigma = -1; ohead != NULL ; ohead = ohead->next) {
	for (i = 0; i < ohead->size; i++) {
	    maxsigma = *(ohead->data+i) > maxsigma ? *(ohead->data+i) : maxsigma;		
	}
    }
    return(maxsigma+1);
}

struct observations **observations_to_array(struct observations *ohead, int *numobs) {
    int i;
    struct observations **obsarray, *o;
    for (o = ohead, i = 0; o != NULL; o = o->next, i++) { }
    *numobs = i;
    obsarray = malloc(sizeof(struct observations *) * i) ;
    for (o = ohead, i = 0; o != NULL; o = o->next, i++) {
	*(obsarray+i) = o;
    }
    return obsarray;
}

struct observations *observations_uniq(struct observations *ohead) {
    struct observations *o, *onext;
    for (o = ohead; o != NULL && o->next != NULL; ) {
	onext = o->next;
	if (onext != NULL) {
	    if (obssortcmp(&o,&onext) == 0) {
		o->next = onext->next;
		free(onext->data);
		free(onext);
		o->occurrences = o->occurrences + 1;
		continue;
	    }
	}
	o = o->next;
    }
    return (ohead);
}

struct observations *observations_sort(struct observations *ohead) {
    struct observations *o, **sptr, *curro, *lasto;
    int i, numobs;
    int (*sorter)() = obssortcmp;
    for (o = ohead, numobs = 0; o != NULL; o = o->next) {
	numobs++;
    }
    if (numobs < 1) {
	return ohead;
    }
    sptr = malloc(sizeof(struct observations *) * numobs);
    for (o = ohead, i = 0; o != NULL; o = o->next, i++) {
	*(sptr+i) = o;
    }
    qsort(sptr, numobs, sizeof(struct observations *), sorter);
    for (i = 1, ohead = lasto = *sptr, lasto->next = NULL; i < numobs; i++) {	
	curro = *(sptr+i);
	curro->next = NULL;
	lasto->next = curro;
	lasto = curro;
    }
    return(ohead);
}

void interrupt_sigproc() {
    fprintf(stderr, "Received SIGINT. Exiting.\n");
    if (g_lastwfsa != NULL) {
	wfsa_print(g_lastwfsa);
    }
    exit(1);
}

void observations_destroy(struct observations *ohead) {
    struct observations *o, *olast;
    for (o = ohead; o != NULL; ) {
	olast = o;
	o = o->next;
	free(olast->data);
	free(olast);
    }
}

struct observations *observations_read(char *filename) {
    char *obs_char_data, *optr;
    struct observations *ohead, *o, *olast;
    int *line, size;
    if ((obs_char_data = file_to_mem(filename)) == NULL) {
	return(NULL);
    }
    ohead = olast = NULL;
    for (optr = obs_char_data, o = ohead ; ;) {
	if (*optr == '\0') {
	    break;
	}
	optr = line_to_int_array(optr, &line, &size);
	if (optr == NULL) {
	    break;
	}
	if (olast == NULL) {
	    ohead = olast = malloc(sizeof(struct observations));
	    o = olast;
	} else {
	    o = malloc(sizeof(struct observations));
	    olast->next = o;
	    olast = o;
	}
	o->size = size;
	o->data = line;	
	o->occurrences = 1;	
	o->next = NULL;	
    }
    free(obs_char_data);
    return(ohead);
}

inline PROB log_add(PROB x, PROB y) {
    PROB temp, negdiff;
    PROB result, result2;
    if (x == LOGZERO) return (y);
    if (y == LOGZERO) return (x);
    if (y > x) {
	temp = x;
	x = y;
	y = temp;
    }

    negdiff = y - x;
    result2 = (double) log2l(exp2l((long double)negdiff)+1);
    if (negdiff <= -61) { return x; }
#ifdef LOG_LUT
    result = log1plus_table_interp(negdiff);
#elif LOG_LIB
    result = log2(exp2(negdiff)+1);
#else
    //result = log1plus_taylor(negdiff);
    result = log1plus_minimax(negdiff);
    //result = log1plus_table_interp(negdiff);
#endif /* LOG_LUT */
    //printf("Asked %.17g REAL: %.17g APPROX: %.17g DIFF: %.17g\n",negdiff,result2,result,fabs(fabs(result)-(fabs(result2))));
    return(result+x);
}

PROB trellis_backward(struct trellis *trellis, int *obs, int length, struct wfsa *fsm) {
    int i, sourcestate, targetstate, symbol;
    PROB target_prob;
    for (i = 0; i <= length + 1; i++)
	for (sourcestate = 0; sourcestate < fsm->num_states; sourcestate++)
	    TRELLIS_CELL(sourcestate,i)->bp = LOGZERO;
    
    /* Fill last and penultimate column */
    for (targetstate = 0; targetstate < fsm->num_states; targetstate++) {
	TRELLIS_CELL(targetstate,length+1)->bp = smrone;
	TRELLIS_CELL(targetstate,length)->bp = smrone + *FINALPROB(fsm, targetstate);
    }
    /* Fill rest */
    for (i = length-1; i >= 0 ; i--) {
	symbol = obs[i];
	for (sourcestate = 0; sourcestate < fsm->num_states; sourcestate++) {
	    TRELLIS_CELL(sourcestate,i)->bp = LOGZERO;
	    for (targetstate = 0; targetstate < fsm->num_states; targetstate++) {
		target_prob = *TRANSITION(fsm, sourcestate, symbol, targetstate);
		if (target_prob <= smrzero) { continue; }
		TRELLIS_CELL(sourcestate,i)->bp = log_add(TRELLIS_CELL(sourcestate,i)->bp, TRELLIS_CELL(targetstate,i+1)->bp + target_prob);
	    }
	}
    }
    return(TRELLIS_CELL(0,0)->bp);
}

PROB trellis_viterbi(struct trellis *trellis, int *obs, int length, struct wfsa *fsm) {
    int i, sourcestate, targetstate, symbol, final_state;
    PROB target_prob, final_prob;
    
    for (i = 0; i <= length + 1; i++)
	for (sourcestate = 0; sourcestate < fsm->num_states; sourcestate++)
	    TRELLIS_CELL(sourcestate,i)->fp = LOGZERO;
    
    /* Calculate first transition */
    TRELLIS_CELL(0,0)->fp = smrone;
    for (i = 0; i < 1 && i < length; i++) {
	symbol = obs[i];
	for (targetstate = 0; targetstate < fsm->num_states; targetstate++) {
	    target_prob = *TRANSITION(fsm, 0, symbol, targetstate);
	    if (target_prob > smrzero) {
		TRELLIS_CELL(targetstate,1)->fp = target_prob;
		TRELLIS_CELL(targetstate,1)->backstate = 0;
	    }
	}
    }
    /* Calculate remaining transitions */
    for (i = 1; i < length; i++) {
	symbol = obs[i];
	for (sourcestate = 0; sourcestate < fsm->num_states; sourcestate++) {
	    if (TRELLIS_CELL(sourcestate,i)->fp == LOGZERO) { continue; }
	    for (targetstate = 0; targetstate < fsm->num_states; targetstate++) {
		target_prob = *TRANSITION(fsm, sourcestate, symbol, targetstate);	       
		if (target_prob <= smrzero) { continue; }
		if (TRELLIS_CELL(targetstate,(i+1))->fp == LOGZERO) {
		    TRELLIS_CELL(targetstate,(i+1))->fp = TRELLIS_CELL(sourcestate,i)->fp + target_prob;
		    TRELLIS_CELL(targetstate,(i+1))->backstate = sourcestate;
		} else {		    
		    if (TRELLIS_CELL(targetstate,(i+1))->fp < TRELLIS_CELL(sourcestate,i)->fp + target_prob) {
			TRELLIS_CELL(targetstate,(i+1))->fp = TRELLIS_CELL(sourcestate,i)->fp + target_prob;
			TRELLIS_CELL(targetstate,(i+1))->backstate = sourcestate;
		    }
		}
	    }
	}
    }
    
    /* Calculate final state probabilities */
    i = length;
    final_state = -1;
    for (targetstate = 0, final_prob = smrzero; targetstate < fsm->num_states; targetstate++) {
	if (TRELLIS_CELL(targetstate,i)->fp == LOGZERO) {
	    TRELLIS_CELL(targetstate,(i+1))->backstate = -1;
	    continue;
	}
	if (*FINALPROB(fsm, targetstate) > smrzero && TRELLIS_CELL(targetstate,i)->fp > smrzero) { 
	    TRELLIS_CELL(targetstate,(i+1))->fp = TRELLIS_CELL(targetstate,i)->fp + *FINALPROB(fsm, targetstate);
	} else {
	    continue;
	}
	if (TRELLIS_CELL(targetstate,(i+1))->fp > final_prob) {
	    final_prob = TRELLIS_CELL(targetstate,(i+1))->fp;
	    final_state = targetstate;
	}
    }
    for (targetstate = 0; targetstate < fsm->num_states; targetstate++) {
	if (targetstate != final_state) {
	    TRELLIS_CELL(targetstate,(i+1))->backstate = -1;
	} else {
	    TRELLIS_CELL(targetstate,(i+1))->backstate = targetstate;
	}
    }
    return(final_prob);
}

PROB trellis_forward(struct trellis *trellis, int *obs, int length, struct wfsa *fsm) {
    int i, sourcestate, targetstate, symbol;
    PROB target_prob, final_prob;

    for (i = 0; i <= length + 1; i++)
	for (sourcestate = 0; sourcestate < fsm->num_states; sourcestate++)
	    TRELLIS_CELL(sourcestate,i)->fp = LOGZERO;
    
    /* Calculate first transition */
    TRELLIS_CELL(0,0)->fp = smrone;
    for (i = 0; i < 1 && i < length; i++) {
	symbol = obs[i];
	for (targetstate = 0; targetstate < fsm->num_states; targetstate++) {
	    target_prob = *TRANSITION(fsm, 0, symbol, targetstate);
	    if (target_prob > smrzero) {
		TRELLIS_CELL(targetstate,1)->fp = target_prob;
		TRELLIS_CELL(targetstate,1)->backstate = 0;
	    }
	}
    }
    /* Calculate remaining transitions */
    for (i = 1; i < length; i++) {
	symbol = obs[i];
	for (sourcestate = 0; sourcestate < fsm->num_states; sourcestate++) {
	    if (TRELLIS_CELL(sourcestate,i)->fp == LOGZERO) { continue; }
	    for (targetstate = 0; targetstate < fsm->num_states; targetstate++) {
		target_prob = *TRANSITION(fsm, sourcestate, symbol, targetstate);
		if (target_prob <= smrzero) { continue; }
		TRELLIS_CELL(targetstate,(i+1))->fp = log_add(TRELLIS_CELL(sourcestate,i)->fp + target_prob, TRELLIS_CELL(targetstate,(i+1))->fp);
	    }
	}
    }
    
    /* Calculate final state probabilities */
    i = length;
    for (targetstate = 0, final_prob = smrzero; targetstate < fsm->num_states; targetstate++) {
	if (TRELLIS_CELL(targetstate,i)->fp == LOGZERO) { continue; }
	if (*FINALPROB(fsm, targetstate) > smrzero && TRELLIS_CELL(targetstate,i)->fp > smrzero) {
	    TRELLIS_CELL(targetstate,(i+1))->fp = TRELLIS_CELL(targetstate,i)->fp + *FINALPROB(fsm, targetstate);
	} else {
	    continue;
	}
	final_prob = log_add(final_prob, TRELLIS_CELL(targetstate,(i+1))->fp);
    }
    return(final_prob);
}

struct trellis *trellis_init(struct observations *o, struct wfsa *fsm) {
    int olenmax;
    struct trellis *trellis;
    for (olenmax = 0 ; o != NULL; o = o->next) {
	olenmax = olenmax < o->size ? o->size : olenmax;
    }
    trellis = calloc((olenmax + 2) * fsm->num_states, sizeof(struct trellis));
    return(trellis);
}

void trellis_print(struct trellis *trellis, struct wfsa *fsm, int obs_len) {
    int i, j;
    for (j = fsm->num_states-1; j >= 0 ; j--) {
	for (i = 0; i <= obs_len + 1 ; i++) {
	    printf("FP:%4.4f BP:%4.4f B:%i\t",  EXP(TRELLIS_CELL(j,i)->fp), EXP(TRELLIS_CELL(j,i)->bp), TRELLIS_CELL(j,i)->backstate);
	}
	printf("\n");
    }
}

void forward_print_path (struct trellis *trellis, struct wfsa *fsm, int obs_len) {
    int i, j, beststate;
    PROB bestprob;
    for (i = 0; i <= obs_len+1; i++) {
	if (i == obs_len)
	    continue;
	bestprob = smrzero;
	beststate = -1;
	for (j = 0; j < fsm->num_states; j++) {
	    if (TRELLIS_CELL(j,i)->fp == LOGZERO) {
		continue;
	    }
	    if (TRELLIS_CELL(j,i)->fp > bestprob) {
		bestprob = TRELLIS_CELL(j,i)->fp;
		beststate = j;
	    }
	}
	printf("%i", beststate);
	if (i < obs_len) printf(" ");
    }
    printf("\n");
}

void backward_print_path (struct trellis *trellis, struct wfsa *fsm, int obs_len) {
    int i, j, beststate;
    PROB bestprob;
    for (i = 0; i <= obs_len; i++) {
	bestprob = smrzero;
	beststate = -1;
	for (j = 0; j < fsm->num_states; j++) {
	    if (TRELLIS_CELL(j,i)->bp == LOGZERO) {
		continue;
	    }
	    if (TRELLIS_CELL(j,i)->bp > bestprob) {
		bestprob = TRELLIS_CELL(j,i)->fp;
		beststate = j;
	    }
	}
	printf("%i", beststate);
	if (i < obs_len) printf(" ");
    }
    printf("\n");
}

void viterbi_print_path (struct trellis *trellis, struct wfsa *fsm, int obs_len) {
    int i, laststate, *path;
    path = malloc(sizeof(int) * (obs_len+1));
    for (i = 0; i < fsm->num_states; i++) {
	if (TRELLIS_CELL(i,obs_len+1)->backstate != -1) {
	    laststate = i;
	}
    }
    *(path+obs_len) = laststate;
    for (i = obs_len; i > 0; i--) {
	*(path+i-1) = TRELLIS_CELL(laststate,i)->backstate;
	laststate = TRELLIS_CELL(laststate,i)->backstate;
    }
    for (i = 0 ; i <= obs_len; i++) {
	printf("%i",path[i]);
	if (i < obs_len) {
	    printf(" ");
	}
    }
    printf("\n");
    free(path);
}

void viterbi(struct wfsa *fsm, struct observations *o, int algorithm) {
    struct observations *obs;
    struct trellis *trellis;
    PROB viterbi_prob;
    trellis = trellis_init(o,fsm);
    for (obs = o; obs != NULL; obs = obs->next) {
	viterbi_prob = trellis_viterbi(trellis, obs->data, obs->size, fsm);
	if (algorithm == DECODE_VITERBI_PROB)
	    printf("%.17g\t", output_convert(viterbi_prob));
	if (algorithm == LIKELIHOOD_VITERBI) 
	    printf("%.17g\n", output_convert(viterbi_prob));
	if (algorithm == DECODE_VITERBI_PROB || algorithm == DECODE_VITERBI) {
	    if (viterbi_prob > smrzero) {
		viterbi_print_path(trellis, fsm, obs->size);
	    } else {
		printf("\n");
	    }
	}
    }
    free(trellis);
}

PROB wfsa_sum_prob(struct wfsa *fsm, int state) {
    PROB sum;
    int i, s;
    /* Get sum of probabilities for a state (in reals) */
    sum = *FINALPROB(fsm,state) >= smrzero ? EXP(*FINALPROB(fsm,state)) : 0;
    for (i = 0; i < fsm->num_states; i++) {
	for (s = 0 ; s < fsm->alphabet_size; s++) {
	    if (*TRANSITION(fsm, state, s, i) >= smrzero) {
		sum += EXP(*TRANSITION(fsm, state, s, i));
	    }
	}
    }
    return(sum);
}

int wfsa_random_transition(struct wfsa *fsm, int state, int *symbol, PROB *prob) {
    /* Choose a random arc from state state           */
    /* Return target state, and put symbol in *symbol */
    /* If stop: symbol = -1                           */
    PROB thissum, r;
    int i, s;
    r = (PROB) random() / RAND_MAX;
    r = r * wfsa_sum_prob(fsm, state);
    thissum = 0;
    
    if (*FINALPROB(fsm,state) >= smrzero) {
	thissum += EXP(*FINALPROB(fsm, state));
	if (thissum >= r) {
	    *symbol = -1;
	    *prob = *FINALPROB(fsm,state);
	    return state;
	}
    }
    for (i = 0; i < fsm->num_states; i++) {
	for (s = 0 ; s < fsm->alphabet_size; s++) {
	    if ((*prob = *TRANSITION(fsm, state, s, i)) >= smrzero) {		
		thissum += EXP(*prob);
		if (thissum >= r) {
		    *symbol = s;
		    return(i);
		}
	    }
	}
    }
    perror("Inconsistent probabilities in FSM");
    exit(1);
}

void generate_words(struct wfsa *fsm, int numwords) {
    int i, j, k, state, symbol, *output = NULL, *stateseq = NULL;
    PROB prob, totalprob;

    output = malloc(sizeof(int) * g_gen_max_length);
    stateseq = malloc(sizeof(int) * g_gen_max_length);
    
    for (i = 0; i < numwords; i++) {
	stateseq[0] = 0;
	totalprob = LOGZERO;
	for (state = 0, j = 0; j < g_gen_max_length ; j++) {
	    state = wfsa_random_transition(fsm, state, &symbol, &prob);
	    totalprob = (totalprob == LOGZERO) ? prob : totalprob + prob;
	    stateseq[j+1] = state;
	    if (symbol >= 0) {
		output[j] = symbol;
	    } else {
		break;
	    }
	}
	printf("%.17g\t",output_convert(totalprob));
	if (j < g_gen_max_length) {
	    for (k = 0; k < j; k++) {
		printf("%i", output[k]);
		if (k < j-1) {
		    printf(" ");
		}
	    }
	    printf("\t");
	    for (k = 0; k <= j; k++) {
		printf("%i", stateseq[k]);
		if (k < j) {
		    printf(" ");
		}
	    }
	    printf("\n");

	} else {
	    i--;
	}
    }
    free(output);
    free(stateseq);
}

void forward(struct wfsa *fsm, struct observations *o, int algorithm) {
    struct observations *obs;
    struct trellis *trellis;
    PROB forward_prob;
    trellis = trellis_init(o,fsm);
    for (obs = o; obs != NULL; obs = obs->next) {
	forward_prob = trellis_forward(trellis, obs->data, obs->size, fsm);
	if (algorithm == DECODE_FORWARD_PROB)
	    printf("%.17g\t", output_convert(forward_prob));
	if (algorithm == LIKELIHOOD_FORWARD)
	    printf("%.17g\n", output_convert(forward_prob));
	if (algorithm == DECODE_FORWARD_PROB || algorithm == DECODE_FORWARD) {
	    if (forward_prob > smrzero) {
		forward_print_path(trellis, fsm, obs->size);
	    } else {
		printf("\n");
	    }	
	}
    }
    free(trellis);
}

void backward(struct wfsa *fsm, struct observations *o, int algorithm) {
    struct observations *obs;
    struct trellis *trellis;
    PROB backward_prob;
    trellis = trellis_init(o,fsm);
    for (obs = o; obs != NULL; obs = obs->next) {
	backward_prob = trellis_backward(trellis, obs->data, obs->size, fsm);
	if (algorithm == DECODE_BACKWARD_PROB)
	    printf("%.17g\t", output_convert(backward_prob));
	if (algorithm == LIKELIHOOD_BACKWARD)
	    printf("%.17g\n", output_convert(backward_prob));
	if (algorithm == DECODE_BACKWARD_PROB || algorithm == DECODE_BACKWARD) {
	    if (backward_prob > smrzero) {
		backward_print_path(trellis, fsm, obs->size);
	    } else {
		printf("\n");
	    }
	}
    }
    free(trellis);
}

PROB train_viterbi(struct wfsa *fsm, struct observations *o, int maxiterations, PROB maxdelta) {
    struct observations *obs;
    struct trellis *trellis;
    int i,j,k,iter, source, target, laststate, symbol, occurrences, *fsm_vit_counts, *fsm_vit_totalcounts, *fsm_vit_finalcounts;
    PROB viterbi_prob, loglikelihood, prevloglikelihood, newprob;
    trellis = trellis_init(o,fsm);
    fsm_vit_counts = malloc(sizeof(int) * fsm->num_states * fsm->num_states * fsm->alphabet_size);
    fsm_vit_totalcounts = malloc(sizeof(int) * fsm->num_states);
    fsm_vit_finalcounts = malloc(sizeof(int) * fsm->num_states);
    
    prevloglikelihood = 0;
    for (iter = 0 ; iter < maxiterations; iter++) {
        /* Clear counts */
        for (i = 0; i < fsm->num_states; i++) {
            fsm_vit_totalcounts[i] = fsm_vit_finalcounts[i] = 0;
        }
        for (i = 0; i < fsm->num_states * fsm->num_states * fsm->alphabet_size; i++) {
            fsm_vit_counts[i] = 0;
        }
        loglikelihood = 0;
        for (obs = o; obs != NULL; obs = obs->next) {
	    occurrences = obs->occurrences;
            viterbi_prob = trellis_viterbi(trellis, obs->data, obs->size, fsm);
            if (viterbi_prob <= smrzero) {
                continue;
            } else {
                loglikelihood += viterbi_prob * occurrences;
                /* Update final counts */
                for (i = 0, laststate = -1; i < fsm->num_states; i++) {
                    if (TRELLIS_CELL(i,(obs->size+1))->backstate != -1) {
                        laststate = i;
                        break;
                    }
                }
                if (laststate == -1) {
                    printf("Could not find last state\n");
                    continue;
                }
                fsm_vit_finalcounts[laststate] += occurrences;
                fsm_vit_totalcounts[laststate] += occurrences;
                /* Update arc counts */
                for (i = obs->size; i > 0; i--) {
                    target = laststate;
                    laststate = TRELLIS_CELL(laststate,i)->backstate;
                    source = laststate;
                    symbol = *((obs->data)+i-1);
                    *FSM_COUNTS(fsm_vit_counts, source, symbol, target) += occurrences;
                    fsm_vit_totalcounts[source]+= occurrences;
                }
            }
        }
	fprintf(stderr, "iteration %i loglikelihood=%.17g delta: %.17g\n", iter+1, loglikelihood, ABS(prevloglikelihood - loglikelihood));
	if (ABS(prevloglikelihood - loglikelihood) < maxdelta)  {
	    break;
	}
        /* Update WFSA */
        for (i=0; i < fsm->num_states; i++) {
	    fsm_vit_totalcounts[i] += g_viterbi_pseudocount * fsm->alphabet_size * fsm->num_states + g_viterbi_pseudocount;
            for (j = 0; j < fsm->alphabet_size; j++) {
                for (k = 0; k < fsm->num_states; k++) {		    
                    if (fsm_vit_totalcounts[i] > 0) {
                        newprob = LOG(*FSM_COUNTS(fsm_vit_counts,i,j,k) + g_viterbi_pseudocount) - LOG(fsm_vit_totalcounts[i]);
                    } else {
                        newprob = smrzero;
                    }
                    *(TRANSITION(fsm,i,j,k)) = newprob;
                }
            }
        }
        for (i=0; i < fsm->num_states; i++) {
            if (fsm_vit_totalcounts[i] > 0) {
                newprob = LOG(fsm_vit_finalcounts[i] + g_viterbi_pseudocount) - LOG(fsm_vit_totalcounts[i]);
            } else {
                newprob = smrzero;
            }
            *(fsm->final_table + i) = newprob;
        }
        prevloglikelihood = loglikelihood;
    }
    free(fsm_vit_counts);
    free(fsm_vit_totalcounts);
    free(fsm_vit_finalcounts);
    return(loglikelihood);
}

inline void spinlock_lock(_Bool *ptr) {
    if (g_num_threads > 1)
    	while (__sync_lock_test_and_set(ptr,1)) { }
}

inline void spinlock_unlock(_Bool *ptr) {
    if (g_num_threads > 1)
     	__sync_lock_release(ptr);
}

void *trellis_fill_bw(void *threadargs) {
    struct trellis *trellis;
    struct observations **obsarray, *obs;
    struct wfsa *fsm;
    PROB backward_prob, forward_prob, thisxi, beta;
    int i, t, symbol, source, target, minobs, maxobs, occurrences;

    trellis = ((struct thread_args *)threadargs)->trellis;
    obsarray = ((struct thread_args *)threadargs)->obsarray;
    minobs = ((struct thread_args *)threadargs)->minobs;
    maxobs = ((struct thread_args *)threadargs)->maxobs;
    fsm = ((struct thread_args *)threadargs)->fsm;
    beta = ((struct thread_args *)threadargs)->beta;
    
    for (i = minobs; i <= maxobs; i++) {
	obs = *(obsarray+i);
	occurrences = obs->occurrences;
	/* E-step */
	backward_prob = trellis_backward(trellis, obs->data, obs->size, fsm);
	forward_prob = trellis_forward(trellis, obs->data, obs->size, fsm);
	pthread_mutex_lock(&mutex1);
	g_loglikelihood += backward_prob * occurrences;
	pthread_mutex_unlock(&mutex1);
	/* Traverse trellis and add */
	for (t = 0; t < obs->size; t++) {
	    symbol = obs->data[t];
	    for (source = 0; source < fsm->num_states; source++) {
		if (TRELLIS_CELL(source,t)->fp == LOGZERO) { continue; }
		for (target = 0; target < fsm->num_states; target++) {
		    if (TRELLIS_CELL(target,t+1)->bp == LOGZERO) { continue; }
		    if (*TRANSITION(fsm,source,symbol,target) <= smrzero) { continue; }
		    thisxi = TRELLIS_CELL(source,t)->fp + *TRANSITION(fsm,source,symbol,target) + TRELLIS_CELL(target,t+1)->bp;
		    thisxi = thisxi - backward_prob;
		    thisxi = g_train_da_bw == 0 ? thisxi : thisxi * beta;
		    thisxi += LOG(occurrences);

		    spinlock_lock(&fsm_counts_spin[source]);
		    fsm_totalcounts[source] = log_add(fsm_totalcounts[source], thisxi);
		    *FSM_COUNTS(fsm_counts,source,symbol,target) = log_add(*FSM_COUNTS(fsm_counts,source,symbol,target), thisxi);
		    spinlock_unlock(&fsm_counts_spin[source]);
		}
	    }
	}
	/* Final states */
	for (source = 0; source < fsm->num_states; source++) {
	    target = source;
	    if (TRELLIS_CELL(source,t)->fp == LOGZERO)   { continue; }
	    if (TRELLIS_CELL(target,t+1)->bp == LOGZERO) { continue; }
	    thisxi = TRELLIS_CELL(source,t)->fp + *FINALPROB(fsm, source);
	    thisxi = thisxi - backward_prob ;
	    thisxi = g_train_da_bw == 0 ? thisxi : thisxi * beta;
	    thisxi += LOG(occurrences);

	    spinlock_lock(&fsm_counts_spin[source]);
	    fsm_totalcounts[source] = log_add(fsm_totalcounts[source], thisxi);
	    fsm_finalcounts[source] = log_add(fsm_finalcounts[source], thisxi);
	    spinlock_unlock(&fsm_counts_spin[source]);
	}
    }
    return(NULL);
}

PROB train_baum_welch(struct wfsa *fsm, struct observations *o, int maxiterations, PROB maxdelta, int vb) {
    struct trellis *trellis, *trellisarray[32];
    struct thread_args *threadargs[32];
    struct observations **obsarray;
    int i, source, target, symbol, iter, numobs, obsperthread;
    PROB newprob, prevloglikelihood, da_beta = 1.0, numstatetrans;
    pthread_t threadids[32];
    
    if (g_train_da_bw) { da_beta = g_betamin; }
    obsarray = observations_to_array(o, &numobs);
    
    /* Each thread gets its own trellis */
    for (i = 0; i < g_num_threads; i++) {
	trellis = trellis_init(o,fsm);
	trellisarray[i] = trellis;
	threadargs[i] = malloc(sizeof(struct thread_args));
    }

    /* These are accessed by all threads through a spinlock */
    fsm_counts = malloc(fsm->num_states * fsm->num_states * fsm->alphabet_size * sizeof(PROB));
    fsm_totalcounts = malloc(fsm->num_states * sizeof(PROB));
    fsm_finalcounts = malloc(fsm->num_states * sizeof(PROB));
    fsm_counts_spin = calloc(fsm->num_states, sizeof(_Bool));
    
    prevloglikelihood = 0;

    obsperthread = (int) ((double) numobs/(double) g_num_threads);
    /* Helper thread parameters (divide observations equally among all) */
    for (i = 1; i < g_num_threads; i++) {
	threadargs[i]->minobs = (i-1) * obsperthread;
	threadargs[i]->maxobs = i * obsperthread - 1;
	threadargs[i]->trellis = trellisarray[i];
	threadargs[i]->obsarray = obsarray;
	threadargs[i]->fsm = fsm;
    }
    /* Main thread parameters */
    threadargs[0]->minobs = (i-1) * obsperthread;
    threadargs[0]->maxobs = numobs - 1;
    threadargs[0]->trellis = trellisarray[0];
    threadargs[0]->obsarray = obsarray;
    threadargs[0]->fsm = fsm;

   
    for (iter = 0 ; iter < maxiterations ; iter++) {
	g_loglikelihood = 0;
	for (i = 0; i < fsm->num_states * fsm->num_states * fsm->alphabet_size ; i++) { fsm_counts[i] = LOGZERO; }
	for (i = 0; i < fsm->num_states ; i++) { fsm_totalcounts[i] =  fsm_finalcounts[i] = LOGZERO; }
	for (i = 1; i < g_num_threads; i++) {
	    /* Launch threads */
	    threadargs[i]->beta = da_beta;
	    pthread_create(&threadids[i], NULL, &trellis_fill_bw, threadargs[i]);
	}

	/* Run main thread Baum-Welch */
	threadargs[0]->beta = da_beta;
	trellis_fill_bw(threadargs[0]);
	/* Wait for all to finish */
	for (i = 1; i < g_num_threads; i++) {
	    pthread_join(threadids[i],NULL);
	}

	if (!g_train_da_bw)
	    fprintf(stderr, "iteration %i loglikelihood=%.17g delta: %.17g\n", iter+1, g_loglikelihood, ABS(prevloglikelihood - g_loglikelihood));
	else 
	    fprintf(stderr, "iteration %i loglikelihood=%.17g delta: %.17g beta: %.17g\n", iter+1, g_loglikelihood, ABS(prevloglikelihood - g_loglikelihood), da_beta);
	    
	if (ABS(prevloglikelihood - g_loglikelihood) < maxdelta)  {
	    if (g_train_da_bw == 1 && da_beta < g_betamax) {
		da_beta *= g_alpha;
		if (da_beta > g_betamax) {
		    da_beta = g_betamax;
		}
	    } else {
		break;
	    }
	}

	/* Modify WFSA (M-step) */
	signal(SIGINT, SIG_IGN); /* Disable interrupts to prevent corrupted WFSA in case of SIGINT while updating */	    
        
	/* if (vb == 1) { */
	/*     for (source = 0; source < fsm->num_states; source++) { */
	/* 	fsm_totalcounts[source] = LOGZERO; */
	/* 	for (symbol = 0; symbol < fsm->alphabet_size; symbol++) {		 */
	/* 	    for (target = 0; target < fsm->num_states; target++) { */
	/* 		newprob = *FSM_COUNTS(fsm_counts,source,symbol,target); */
	/* 		if (newprob == LOGZERO) { */
	/* 		    continue; */
	/* 		} */
	/* 		newprob = digamma(myexp(log_add(newprob, LOG(g_vb_alpha))))/M_LN2; */
	/* 		*FSM_COUNTS(fsm_counts,source,symbol,target) = newprob; */
	/* 		fsm_totalcounts[source] = log_add(fsm_totalcounts[source], newprob); */
	/* 	    } */
	/* 	} */
	/*     } */
	/*     for (source = 0; source < fsm->num_states; source++) { */
	/* 	newprob = fsm_finalcounts[source]; */
	/* 	if (newprob == LOGZERO) { */
	/* 	    continue; */
	/* 	}		 */
	/* 	newprob = digamma(myexp(log_add(newprob, LOG(g_vb_alpha))))/M_LN2; */
	/* 	fsm_totalcounts[source] = log_add(fsm_totalcounts[source], newprob); */
	/* 	fsm_finalcounts[source] = newprob; */
	/*     } */
	/* } */
	numstatetrans = fsm->num_states * fsm->alphabet_size + 1;

	for (source = 0; source < fsm->num_states; source++) {
	    for (symbol = 0; symbol < fsm->alphabet_size; symbol++) {
		for (target = 0; target < fsm->num_states; target++) {
		    newprob = *FSM_COUNTS(fsm_counts,source,symbol,target);
		    if (newprob == LOGZERO) { newprob = smrzero; }
		    /* Variational Bayes */
		    if (vb) {
			*TRANSITION(fsm,source,symbol,target) = (digamma(EXP(newprob) + g_vb_alpha) - digamma(EXP(fsm_totalcounts[source]) + numstatetrans * g_vb_alpha))/M_LN2;
		    } else {
			*TRANSITION(fsm,source,symbol,target) = newprob - fsm_totalcounts[source];
		    }
		}
	    }
	}
	for (source = 0; source < fsm->num_states; source++) {
	    newprob = fsm_finalcounts[source];
	    if (newprob == LOGZERO) { newprob = smrzero; }
	    /* Variational Bayes */
	    if (vb) {
		*FINALPROB(fsm,source) = (digamma(EXP(newprob) + g_vb_alpha) - digamma(EXP(fsm_totalcounts[source]) + numstatetrans * g_vb_alpha))/M_LN2;
	    } else {
		*FINALPROB(fsm,source) = newprob - fsm_totalcounts[source];
	    }
	}
	g_lastwfsa = fsm;                          /* Put fsm into global var to recover in case of SIGINT */
	signal(SIGINT, (void *)interrupt_sigproc); /* Re-enable interrupt */
	prevloglikelihood = g_loglikelihood;
    }

    for (i = 0; i < g_num_threads; i++) {
	free(trellisarray[i]);
	free(threadargs[i]);
    }
    free(fsm_counts);
    free(fsm_totalcounts);
    free(fsm_finalcounts);
    free(fsm_counts_spin);
    free(obsarray);
    return(g_loglikelihood);
}

PROB train_bw(struct wfsa *fsm, struct observations *o, int maxiterations, PROB maxdelta) {
    struct wfsa *bestfsm;
    int i;
    PROB thisll, minll;
    if (g_random_restarts == 0)
	return(train_baum_welch(fsm,o,maxiterations,maxdelta,g_bw_vb));

    /* Random restarts */
    minll = -DBL_MAX;
    bestfsm = NULL;
    for (i = 0; i < g_random_restarts; i++) {
	fprintf(stderr, "===Running BW restart #%i===\n", i+1);
	thisll = train_baum_welch(fsm, o, g_random_restart_iterations, maxdelta,g_bw_vb);
	if (thisll > minll) {
	    minll = thisll;
	    if (bestfsm != NULL)
		wfsa_destroy(bestfsm);
	    bestfsm = wfsa_copy(fsm);
	}
	if (i < g_random_restarts + 1) {
	    fsm = wfsa_init(g_num_states, g_alphabet_size);
	    if (g_generate_type == GENERATE_NONDETERMINISTIC)
		wfsa_randomize_nondeterministic(fsm,0,0);
	    if (g_generate_type == GENERATE_DETERMINISTIC)
		wfsa_randomize_deterministic(fsm,0);
	    if (g_generate_type == GENERATE_BAKIS)
		wfsa_randomize_nondeterministic(fsm,1,0);
	    wfsa_to_log2(fsm);
	}
    }
    fprintf(stderr, "===Running final BW===\n");
    return(train_baum_welch(bestfsm,o,maxiterations,maxdelta,g_bw_vb));
}

PROB train_viterbi_bw(struct wfsa *fsm, struct observations *o) {
    PROB ll;
    train_viterbi(fsm,o,g_maxiterations,g_maxdelta);
    ll = train_baum_welch(fsm,o,g_maxiterations,g_maxdelta,g_bw_vb);
    return(ll);
}

struct gibbs_state_chain {
    int state;
    int sym;
};

struct gibbs_state_chain *gibbs_init(struct observations *o, int num_states, int alphabet_size, int *obslen) {
    int i,j,*data;
    struct gibbs_state_chain *go;
    struct observations *tempo;

    /* Initialize the sequence of observations for Gibbs sampling */
    /* It's a simple array of two ints, like so:                  */

    /*        0     ...  obslen-1                                 */ 
    /*     --------     --------                                  */
    /*     |state |     |state |                                  */
    /*     |sym   |     |sym   |                                  */
    /*     --------     --------                                  */
    
    /* We chain together multiple observations into one big       */
    /* observation by adding the symbol # at the end of each o    */
    /* which transitions to state 0. State 0 after # is never     */
    /* resampled. Also, we augment the alphabet by one to         */
    /* accommodate for #. After sampling, # represents the        */
    /* halting probability at that state, and is removed when     */
    /* constructing an automaton based on the samples taken.      */

    for (tempo = o, (*obslen) = 0; tempo != NULL; tempo = tempo->next, (*obslen)++) {
	(*obslen) += tempo->size;
    }
    (*obslen)++;

    go = malloc(sizeof(struct gibbs_state_chain) * (*obslen));

    for (j = 0; o != NULL; o = o->next) {
	data = o->data;
	/* First state is always 0 */
	(go+j)->state = 0;
	for (i = 0; i < o->size; i++) {
	    if (i > 0) {
		(go+j)->state = rand_int_range(0,num_states);
	    }
	    /* Add transition */
	    (go+j)->sym = *(data+i);
	    j++;
	}
	(go+j)->state = rand_int_range(0,num_states);
	(go+j)->sym = alphabet_size;
	j++;
    }
    (go+j)->state = 0;
    (go+j)->sym = -1; /* sentinel */
    return(go);
}

int gibbs_weighted_select(PROB *weight_list, PROB weight_sum, int num_elements) {

    /* Do a weighted selection of num_elements, each with their own */
    /* weight specified in the array prob_list                      */
    /* The sum of all weights should be given in weight_sum         */

    int i;
    PROB select;
    select = rand_double() * weight_sum;
    for (i = 0; i < num_elements; i++) {
	select -= weight_list[i];
	if (select < 0)
	    return i;
    }
    perror("Something is rotten with gibbs_weighted_select(). Shouldn't have reached here.");
    return i;
}

PROB gibbs_sampler(struct wfsa *fsm, struct observations *o, double beta, int num_states, int maxiter, int burnin, int lag) {

  /* Macros to access counts easily */
  #define Ccurr(SOURCE_STATE, SYMBOL, TARGET_STATE) (*((gibbs_counts) + (num_states * alphabet_size * (SOURCE_STATE) + (SYMBOL) * num_states + (TARGET_STATE))))

  #define Ctimestamp(SOURCE_STATE, SYMBOL, TARGET_STATE) (*((gibbs_timestamps) + (num_states * alphabet_size * (SOURCE_STATE) + (SYMBOL) * num_states + (TARGET_STATE))))

  #define Csampled(SOURCE_STATE, SYMBOL, TARGET_STATE) (*((gibbs_sampled_counts) + (num_states * alphabet_size * (SOURCE_STATE) + (SYMBOL) * num_states + (TARGET_STATE))))

    int i, j, k, l, obslen, alphabet_size, z, zprev, znext, a, aprev, anext, indicator, newstate, samplecount;
    unsigned int steps, *gibbs_counts, *gibbs_counts_states, *gibbs_sampled_counts, *gibbs_counts_sampled_states, *gibbs_timestamps;
    PROB ANbeta, g_k, g_sum, *current_prob;
    struct gibbs_state_chain *chain;
	
    alphabet_size = g_alphabet_size + 1; /* Use extra symbol for end-of-word (#) */
    /* Build initial array of states */
    chain = gibbs_init(o, num_states, g_alphabet_size, &obslen);
    gibbs_counts = calloc(num_states * num_states * alphabet_size, sizeof(int));
    gibbs_sampled_counts = calloc(num_states * num_states * alphabet_size, sizeof(int));
    gibbs_timestamps = calloc(num_states * num_states * alphabet_size, sizeof(int));
    gibbs_counts_states = calloc(num_states, sizeof(int));
    gibbs_counts_sampled_states = calloc(num_states, sizeof(int));
    current_prob = calloc(num_states, sizeof(PROB));

    /* Initialize timestamps to 0 */
    for (i = 0 ; i < num_states * num_states * alphabet_size; i++)
	gibbs_timestamps[i] = 0;

    /* Accumulate initial counts from initial random sequence */
    for (i = 0; i < obslen-1; i++) {
	Ccurr( (chain+i)->state , (chain+i)->sym, (chain+i+1)->state)++;
	gibbs_counts_states[(chain+i)->state]++;
    }

    ANbeta = alphabet_size * num_states * beta;

    for (i = 0, samplecount = 1, steps = 0; i < maxiter; i++) {
	if (i > 0 && i % 100 == 0) {
	    fprintf(stderr, "Iteration: %i  Samples collected: %i\n", i, samplecount-1);
	}
	for (j = 0; j < obslen; j++) {
	    if (j == 0 || (chain+j-1)->sym == g_alphabet_size) { /* Don't resample "initial" states */
		continue;
	    }
	    a = (chain+j)->sym;          /* Current symbol  */
	    aprev = (chain+j-1)->sym;    /* Previous symbol */
	    z = (chain+j)->state;        /* Current state   */
	    zprev = (chain+j-1)->state;  /* Previous state  */
	    znext = (chain+j+1)->state;  /* Next state      */

	    /* Settle up with total counts depending on timestamp */
	    if (Ctimestamp(zprev,aprev,z) != samplecount) {
		Csampled(zprev,aprev,z) += ( samplecount - Ctimestamp(zprev,aprev,z) ) * Ccurr(zprev,aprev,z);
		Ctimestamp(zprev,aprev,z) = samplecount;
	    }
	    if (Ctimestamp(z,a,znext) != samplecount) {
		Csampled(z,a,znext) += ( samplecount - Ctimestamp(z,a,znext) ) * Ccurr(z,a,znext);
		Ctimestamp(z,a,znext) = samplecount;
	    }

	    /* Adjust counts */
	    Ccurr(zprev,aprev,z)--;
	    Ccurr(z,a,znext)--;
	    Csampled(zprev,aprev,z)--;
	    Csampled(z,a,znext)--;

	    gibbs_counts_states[z]--;

	    for (k = 0, g_sum = 0; k < num_states; k++) {
		/* Sample */
		indicator = (k == zprev && aprev == a && znext == k) ? 1 : 0;
		g_k = (((double)Ccurr(k,a,znext)) + beta + indicator) * (((double)Ccurr(zprev,aprev,k)) + beta) / ((double)gibbs_counts_states[k] + ANbeta);
		if (g_k < 0) {
		    exit(0);
		}
		current_prob[k] = g_k;
		g_sum += g_k;
	    }
	    /* Do #num_states-way weighted coin toss to select new state */
	    newstate = gibbs_weighted_select(current_prob, g_sum, num_states);

	    /* Settle up with total counts depending on timestamp */
	    if (Ctimestamp(zprev,aprev,newstate) != samplecount) {
		Csampled(zprev,aprev,newstate) += ( samplecount - Ctimestamp(zprev,aprev,newstate) ) * Ccurr(zprev,aprev,newstate);
		Ctimestamp(zprev,aprev,newstate) = samplecount;
	    }
	    if (Ctimestamp(newstate,a,znext) != samplecount) {
		Csampled(newstate,a,znext) += ( samplecount - Ctimestamp(newstate,a,znext) ) * Ccurr(newstate,a,znext);
		Ctimestamp(newstate,a,znext) = samplecount;
	    }
	    Ccurr(zprev,aprev,newstate)++;
	    Ccurr(newstate,a,znext)++;
	    Csampled(zprev,aprev,newstate)++;
	    Csampled(newstate,a,znext)++;

	    gibbs_counts_states[newstate]++;
	    (chain+j)->state = newstate;
	    steps++;
	    if (i > burnin && steps % lag == 0) {
		/* Update samplecount: sample count table entries will be updated */
		/* as we go along based on the timestamp                          */
		samplecount++;
	    }
	}
    }
    /* Settle all remaining debts Ccurr has with Csampled depending on timestamps */
    for (l = 0; l < num_states * num_states * alphabet_size; l++) {
	gibbs_sampled_counts[l] += ( samplecount - gibbs_timestamps[l] ) * gibbs_counts[l];
    }

    /* Construct WFSA */
    for (i = 0 ; i < num_states; i++) {
	for (j = 0; j < num_states; j++) {
	    for (k = 0; k < alphabet_size; k++) {
		gibbs_counts_sampled_states[i] += Csampled(i,k,j);
	    }
	}
    }
    for (i = 0 ; i < num_states; i++) {
	for (j = 0; j < num_states; j++) {
	    for (k = 0; k < alphabet_size; k++) {
		if (k == alphabet_size - 1) {
		    *FINALPROB(fsm,i) += ((double)Csampled(i,k,j) + beta) / ((double)gibbs_counts_sampled_states[i] + ANbeta);
		} else {		    
		    *(TRANSITION(fsm,i,k,j)) = log2(((double)Csampled(i,k,j) + beta) / ((double)gibbs_counts_sampled_states[i] + ANbeta));
		}
	    }
	}
    }
    for (i = 0; i < num_states; i++)
	*FINALPROB(fsm,i) = log2(*FINALPROB(fsm,i));

    free(chain);
    free(gibbs_counts);
    free(gibbs_sampled_counts);
    free(gibbs_counts_states);
    free(gibbs_counts_sampled_states);
    free(current_prob);
    return(0.0); /* log likelihood */
}

int main(int argc, char **argv) {
    int opt, algorithm = 0, numelem, obs_alphabet_size;
    char *fsmfile = NULL, optionchar;
    struct wfsa *fsm;
    struct observations *o = NULL;
	
    log1plus_taylor_init();
#ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    g_numcpus = sysinfo.dwNumberOfProcessors;
#else
    g_numcpus = (int) sysconf(_SC_NPROCESSORS_ONLN);
#endif

    srandom((unsigned int)time((time_t *)NULL));
    srand48((unsigned int)time((time_t *)NULL));

    while ((opt = getopt(argc, argv, "uhva:b:d:f:l:i:o:p:r:g:t:x:G:T:L:D:")) != -1) {
	switch(opt) {
	case 'v':
	    printf("This is %s\n", versionstring);
	    exit(0);
	case 'h':
	    printf("%s", helpstring);
	    exit(0);
	case 'u':
	    g_initialize_uniform = 1;
	    break;
	case 'a':
#ifdef MATH_FLOAT
	    numelem = sscanf(optarg,"%g,%g,%g",&g_betamin,&g_betamax,&g_alpha);
#else
	    numelem = sscanf(optarg,"%lg,%lg,%lg",&g_betamin,&g_betamax,&g_alpha);
#endif
	    if (numelem < 3) {
		printf("-a option requires betamin,betamax,alpha\n"); 
		exit(1);
	    }
	    break;
	case 'r':
	    numelem = sscanf(optarg,"%i,%i",&g_random_restarts,&g_random_restart_iterations);
	    if (numelem < 1) {
		printf("-r option requires #num restarts[,#num iterations per restart]\n"); 
		exit(1);
	    }
	    break;
	case 'g':
	    numelem = sscanf(optarg,"%c%i,%i",&optionchar,&g_num_states,&g_alphabet_size);
	    switch(optionchar) {
	    case 'b': g_generate_type = GENERATE_BAKIS; break;
	    case 'd': g_generate_type = GENERATE_DETERMINISTIC; break;
	    case 'n': g_generate_type = GENERATE_NONDETERMINISTIC; break;
	    default: printf("-g option requires b|d|n\n"); exit(1);
	    }
			
	    if (numelem < 2) { printf("Error in option -g\n"); exit(1); }
	    if (numelem < 3) { g_alphabet_size = -1;}
	    break;
	case 't':
	    if (strncmp(optarg,"c/",2) == 0) {
		g_num_threads = g_numcpus / atoi(optarg+2);
	    } else if (strncmp(optarg,"c",1) == 0) {
		g_num_threads = g_numcpus - atoi(optarg+1);
	    } else {
		g_num_threads = atoi(optarg);
	    }
	    if (g_num_threads <= 0) { g_num_threads = 1; }
	    break;
	case 'd':
	    g_maxdelta = strtod(optarg,NULL);
	    break;
	case 'i':
	    if (strcmp(optarg,"real") == 0)   {  g_input_format = FORMAT_REAL;   }
	    if (strcmp(optarg,"log10") == 0)  {  g_input_format = FORMAT_LOG10;  }
	    if (strcmp(optarg,"ln") == 0)     {  g_input_format = FORMAT_LN;     }
	    if (strcmp(optarg,"log2") == 0)   {  g_input_format = FORMAT_LOG2;   }
	    if (strcmp(optarg,"nlog10") == 0) {  g_input_format = FORMAT_NLOG10; }
	    if (strcmp(optarg,"nln") == 0)    {  g_input_format = FORMAT_NLN;    }
	    if (strcmp(optarg,"nlog2") == 0)  {  g_input_format = FORMAT_NLOG2;  }
	    break;
	case 'o':
	    if (strcmp(optarg,"real") == 0)   {  g_output_format = FORMAT_REAL;   }
	    if (strcmp(optarg,"log10") == 0)  {  g_output_format = FORMAT_LOG10;  }
	    if (strcmp(optarg,"ln") == 0)     {  g_output_format = FORMAT_LN;     }
	    if (strcmp(optarg,"log2") == 0)   {  g_output_format = FORMAT_LOG2;   }
	    if (strcmp(optarg,"nlog10") == 0) {  g_output_format = FORMAT_NLOG10; }
	    if (strcmp(optarg,"nln") == 0)    {  g_output_format = FORMAT_NLN;    }
	    if (strcmp(optarg,"nlog2") == 0)  {  g_output_format = FORMAT_NLOG2;  }
	    break;
	case 'f':
	    fsmfile = strdup(optarg);
	    break;
	case 'p':
	    g_viterbi_pseudocount = atoi(optarg);
	    g_vb_alpha = strtod(optarg,NULL);
	    g_gibbs_beta = strtod(optarg,NULL);			
	    break;
	case 'b':
	    g_gibbs_burnin = atoi(optarg);
	    break;
	case 'l':
	    g_gibbs_lag = atoi(optarg);
	    break;
	case 'x':
	    g_maxiterations = atoi(optarg);
	    break;
	case 'G':
	    g_generate_words = atoi(optarg);
	    algorithm = GENERATE_WORDS;
	    break;
	case 'T':
	    if (strcmp(optarg,"vit") == 0)
		algorithm = TRAIN_VITERBI;
	    if (strcmp(optarg,"bw") == 0)
		algorithm = TRAIN_BAUM_WELCH;
	    if (strcmp(optarg,"gs") == 0)
		algorithm = TRAIN_GIBBS_SAMPLING;
	    if (strcmp(optarg,"vb") == 0) {
		algorithm = TRAIN_VARIATIONAL_BAYES;
		g_bw_vb = 1;
	    }
	    if (strcmp(optarg,"vitbw") == 0)
		algorithm = TRAIN_VITERBI_BW;
	    if (strcmp(optarg,"dabw") == 0) {
		algorithm = TRAIN_DA_BAUM_WELCH;
		g_train_da_bw = 1;
	    }
	    break;
	case 'L':
	    if (strcmp(optarg,"vit") == 0)
		algorithm = LIKELIHOOD_VITERBI;
	    if (strcmp(optarg,"f") == 0)
		algorithm = LIKELIHOOD_FORWARD;
	    if (strcmp(optarg,"b") == 0)
		algorithm = LIKELIHOOD_BACKWARD;
	    break;
	case 'D':
	    if (strcmp(optarg,"vit") == 0)
		algorithm = DECODE_VITERBI;
	    if (strcmp(optarg,"vit,p") == 0)
		algorithm = DECODE_VITERBI_PROB;
	    if (strcmp(optarg,"f") == 0)
		algorithm = DECODE_FORWARD;
	    if (strcmp(optarg,"f,p") == 0)
		algorithm = DECODE_FORWARD_PROB;
	    if (strcmp(optarg,"b") == 0)
		algorithm = DECODE_BACKWARD;
	    if (strcmp(optarg,"b,p") == 0)
		algorithm = DECODE_BACKWARD_PROB;
	    break;
	}
    }
    argc -= optind;
    argv += optind;

    if (argc < 1 && ((algorithm && algorithm != GENERATE_WORDS) || (g_alphabet_size < 0 && fsmfile == NULL))) {
	printf("Missing observation filename\n");
	printf("Usage: %s",usagestring);
	exit(1);
    }
    if (argc > 0) {
	if ((o = observations_read(argv[0])) == NULL) {
	    perror("Error reading observations file");	    
	    exit(1);
	}
	obs_alphabet_size = observations_alphabet_size(o);
    }
    if (fsmfile == NULL && g_generate_type == 0) {
	printf("You must either specify a FSM file with -f, or initialize a random FSM with -g\n");
	exit(1);
    }
    if (g_alphabet_size < 0 && o != NULL) {
	g_alphabet_size = obs_alphabet_size;
    }
    if (g_generate_type > 0) {
	fsm = wfsa_init(g_num_states, g_alphabet_size);
	if (algorithm != TRAIN_GIBBS_SAMPLING) {
	    if (g_generate_type == GENERATE_NONDETERMINISTIC)
		wfsa_randomize_nondeterministic(fsm,0,g_initialize_uniform);
	    if (g_generate_type == GENERATE_DETERMINISTIC)
		wfsa_randomize_deterministic(fsm,g_initialize_uniform);
	    if (g_generate_type == GENERATE_BAKIS)
		wfsa_randomize_nondeterministic(fsm,1,g_initialize_uniform);
	    wfsa_to_log2(fsm);
	}
    }

    if (fsmfile != NULL) {
	fsm = wfsa_read_file(fsmfile);
	if (g_input_format != FORMAT_LOG2) {
	    wfsa_to_log2(fsm);
	}
    }
    if (o != NULL && fsm->alphabet_size < obs_alphabet_size) {
	printf("Error: the observations file has symbols outside the FSA alphabet.\n");
	printf("Observations alphabet size %i, FSA alphabet size: %i.\n", fsm->alphabet_size, obs_alphabet_size);
	exit(1);
    }
	
    log1plus_init();
	
    switch (algorithm) {
    case GENERATE_WORDS:
	generate_words(fsm,g_generate_words);
	break;
    case DECODE_VITERBI:
    case DECODE_VITERBI_PROB:
    case LIKELIHOOD_VITERBI:
	viterbi(fsm,o,algorithm);
	break;
    case DECODE_FORWARD:
    case DECODE_FORWARD_PROB:
    case LIKELIHOOD_FORWARD:
	forward(fsm,o,algorithm);
	break;
    case DECODE_BACKWARD:
    case DECODE_BACKWARD_PROB:
    case LIKELIHOOD_BACKWARD:
	backward(fsm,o,algorithm);
	break;
    case TRAIN_VARIATIONAL_BAYES:
    case TRAIN_BAUM_WELCH:
    case TRAIN_DA_BAUM_WELCH:
	o = observations_sort(o);
	o = observations_uniq(o);
	train_bw(fsm,o,g_maxiterations,g_maxdelta);
	wfsa_print(fsm);
	break;
    case TRAIN_VITERBI:
	o = observations_sort(o);
	o = observations_uniq(o);
	train_viterbi(fsm,o,g_maxiterations,g_maxdelta);
	wfsa_print(fsm);
	break;
    case TRAIN_VITERBI_BW:
	o = observations_sort(o);
	o = observations_uniq(o);
	train_viterbi_bw(fsm,o);
	wfsa_print(fsm);
	break;
    case TRAIN_GIBBS_SAMPLING:
	o = observations_sort(o);
	gibbs_sampler(fsm, o, g_gibbs_beta, g_num_states, g_maxiterations, g_gibbs_burnin, g_gibbs_lag);
	wfsa_print(fsm);
	break;
    default:
	wfsa_print(fsm);
    }
    log1plus_free();
    if (fsm != NULL)
	wfsa_destroy(fsm);
    if (o != NULL)
	observations_destroy(o);
    exit(0);
}
