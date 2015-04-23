# News #

  * Treba 1.01 released (11/27/2013):
  * Now supports **Baum-Welch** training for HMMs (with parallel threads for multiple CPUs). Also adds **Viterbi training** and **Variational Bayes** for HMMs.

  * Treba 1.0 released (10/30/2013):
  * Adds **Gibbs sampler**/likelihood calculation support for for HMMs
  * Features new algorithms for Gibbs sampling (PFSA/HMMs)
  * Contains **NVIDIA CUDA** GPU implementation of Gibbs sampling for PFSA/HMMs, yielding up to 100x speedups over the CPU-based C implementation with NVIDIA cards
  * New algorithms for deterministic PFSA induction: ALERGIA/MDI and variations thereof
  * Now supports **Variational Bayes**

# Summary #

Treba is a basic command-line tool for training, decoding, and calculating with Hidden Markov Models (HMMs) and weighted (probabilistic) finite state automata (PFSA). It

  * trains PFSA/HMMs, decodes, calculates sequence probabilities and generates sequences of PFSA.
  * includes a parallel implementation of HMM/PFSA training with Baum-Welch (optionally with deterministic annealing) for multiple CPUs/cores.
  * includes a Gibbs sampler for estimating probabilities of PFSA/HMMs.
  * supports inference of deterministic PFSA using state-merging algorithms with various knobs to tweak (ALERGIA/MDI)
  * uses a simple text format compatible with AT&T finite-state tools or OpenFST.
  * is written in C.
  * is very fast: the CUDA-enabled Gibbs sampler produces state-of-the-art inference accuracy in a few seconds on large problems in the [PAutomaC](http://ai.cs.umbc.edu/icgi2012/challenge/Pautomac/) dataset.

# Installation #

  * Download the latest tarball
  * Compile source and install with `make; sudo make install`
  * If you have the [CUDA](https://developer.nvidia.com/cuda-downloads) developer toolkit installed, you can compile with CUDA support by `make CUDA=1`; `sudo make install`

or

  * Run one of the pre-compiled binaries for OSX, Win32 and Linux (i386 and x86\_64) found in the tarball in the bin/ directory

# Overview #

Treba has four main modes of operation illustrated in the following examples:

  * (1) **Training**: `treba --train=ALGORITHM [options] OBSERVATION-FILE`
  * Train a PFSA with ALGORITHM using the observations in OBSERVATION-FILE. Or, train a HMM, by adding the flag `--hmm`

  * (2) **Likelihood**: `treba --likelihood=f --file fsmfile.fsm OBSERVATION-FILE`
  * Calculate the probability of each sequence (in this case the forward probability) in `OBSERVATION-FILE` with respect to the PFSA in `fsmfile.fsm`. You can also do the same for HMMs by adding the flag `--hmm`

  * (3) **Decoding**: `treba --decode=vit --file=fsmfile.fsm OBSERVATION-FILE`
  * Calculate the most likely state-path (the Viterbi path) for sequences in OBSERVATION-FILE according to the PFSA `fsmfile.fsm` (or a HMM file)

  * (4) **Generation**: `treba --generate=NUM --file=fsmfile.fsm`
  * Generate `NUM` sequences at random from the PFSA in `fsmfile.fsm` (or a HMM file), weighted by its transition probabilities

(1)--(3) correspond to the three classical problems regarding HMMs (Rabiner, 1989): adjusting the model parameters given known observation sequences, calculating the probability of a sequence, and calculating the optimal state sequence for an observation.


---


# Brief tutorial #

The following tutorial explains and illustrates common usage patterns and the file formats of `treba` - this includes training PFSA and HMMs with EM (Baum-Welch), a Gibbs sampler, as well as likelihood calculations, decoding, and inferring deterministic PFSA using state-merging algorithms. For details on allowed options, see the [man page](TrebaManPage.md).

## Calculating the probability of sequences with PFSA and HMMs ##

## File Formats ##

### PFSA file format ###

PFSA are represented in a simple text format.  For example, the following PFSA


<img src='http://treba.googlecode.com/svn/wiki/automatonexample.png' alt='Logo' width='50%' />


is represented by the file

```
# myfsm.fsm
0 0 0 0.15
0 0 1 0.06
0 0 2 0.09
0 1 0 0.35
0 1 1 0.14
0 1 2 0.21
1 1 2 0.18
1 1 3 0.63
1 1 4 0.09
1 2 2 0.02
1 2 3 0.07
1 2 4 0.01
2 2 3 0.04
2 2 5 0.2
2 2 6 0.16
2 3 3 0.06
2 3 5 0.30
2 3 6 0.24
3 1.0
```

That is, each line is either of the format:

```
SOURCE-STATE TARGET-STATE SYMBOL PROBABILITY
```

or

```
STATE HALTING-PROBABILITY
```

Note that the probabilities for each state (outgoing transitions + the state's halting probability) are normalized to sum to one.  This is advisable, though not strictly necessary.  The only possible initial state is state 0.

Also note that "symbols" are always integers.  If you need to use alphabetic symbols, you need to do the mapping to integers yourself before using `treba`.

### HMM file format ###

The HMM file format is similar, and it is assumed that state 0 is a designated **start** state, from which emissions never occur. Likewise, the highest numbered state is assumed to be the **end** state, with no emissions or transitions to other states.  This means, for instance, that emitting the empty string only is possible, and its corresponding path will be **start**-**end**.

The below 2-symbol, 4-state HMM, for instance

<img src='http://treba.googlecode.com/svn/wiki/hmmmodels.png' alt='Logo' width='50%' />

is represented by the file:

```
# myhmm.hmm
0 > 1 0.44
0 > 2 0.03
0 > 3 0.53
1 > 1 0.11
1 > 2 0.08
1 > 3 0.81
2 > 1 0.48
2 > 2 0.25
2 > 3 0.27
1 0 0.21
1 1 0.79
2 0 0.92
2 1 0.08
```

Here, each line is either of the format:

```
SOURCE-STATE > TARGET-STATE TRANSITION-PROBABILITY
```

or

```
STATE SYMBOL EMISSION-PROBABILITY
```

Note that you may use weights, or logarithms of probabilities in file specifications (which treba uses internally for calculations) instead of real numbers. However, this needs to be specified with the `--input-format` option (see below).


### Observation file format ###

To train or to calculate the probability of observation sequences, we need to store them in a file, one line for each observation, which are whitespace-separated numbers representing symbols:

```
# myobs.obs
0 2 2 3 3
0 1 2 3
```

### Likelihoods (PFSA) ###

Now, we can calculate the probability of the two above sequences with respect to the above automaton:

```
treba --likelihood=f --file myfsm.fsm myobs.obs
```

which calculates the forward probability for each observation (the sum of the probabilities of all paths through the automaton), producing

```
8.7885015411769804e-05
2.5199999999999972e-05
```


### Likelihoods (HMMs) ###

The observation file format is the same regardless of whether you're using PFSA or HMMs.  However, when doing calculations, you must specify that you're assuming a HMM as input as PFSA is the default. For example, to calculate the likelihood of the two sequences

```
# myobs2.obs
0 1 1 1 1 0
0 1 1
```

with respect to the above HMM, you can do:

```
treba --hmm --likelihood=f --file=myhmm.hmm myobs2.obs
```

yielding

```
1.2647161259981213e-06
0.0016911986494436745
```

### The most likely path (PFSA) ###

We can also calculate the Viterbi probability (the probability of the most likely path) by

```
treba --likelihood=vit --file myfsm.fsm myobs.obs
```

producing

```
4.7627999999999965e-05
2.5199999999999972e-05
```

Note that the Viterbi probability for the second observation (**0123**) is identical to the forward probability since the sequence **0123** only has one possible path through the automaton.

### The most likely path (HMMs) ###

As with PFSA, you can also find the probability of the most likely path through a HMM:

```
treba --hmm --likelihood=vit --file myhmm.hmm myobs2.obs
```

producing

```
1.3648292411249901e-07
0.00073668564287999983
```


## Decoding (PFSA) ##

If we also want to know the actual most likely paths for each observation sequence, we can run:

```
treba --decode=vit --file myfsm.fsm myobs.obs
```

producing
```
0 1 1 1 2 3
0 0 1 2 3
```

This will print out the Viterbi (most likely) paths.  If we want to see the probability as well as the path (TAB-separated), we can issue:

```
treba --decode=vit,p --file myfsm.fsm myobs.obs
```

producing

```
4.7627999999999965e-05	0 1 1 1 2 3
2.5199999999999972e-05	0 0 1 2 3
```

## Decoding (HMMs) ##

Decoding for HMMs works similarly. For example,

```
treba --hmm --decode=vit,p --file myhmm.hmm myobs2.obs
```

will produce

```
1.3648292411249901e-07	0 2 1 1 1 1 2 3
0.00073668564287999983	0 2 1 1 3
```

### Forward path/posterior decoding (PFSA) ###

There is also the possibility of printing the "forward path".  This is the state sequence where for each step in time, we choose the most likely state to be in at that time.  Note that for automata that are not fully connected, this will often not correspond to an actually valid path through the automaton, as in the below:

```
treba --decode=f --file myfsm.fsm myobs.obs
```

and the output is

```
0 1 1 1 1 3
0 1 1 1 3
```

### Forward path/posterior decoding (HMMs) ###

This is similar for HMMs, for example:

```
treba --hmm --decode=f,p --file myhmm.hmm myobs2.obs
```

will print out the best "forward path", along with the probability (since we used the `--decode=f,p` option):

```
1.2647161259981213e-06	0 1 1 1 1 1 1 3
0.0016911986494436745	0 1 1 1 3
```


## Training a PFSA/WFSA ##

### EM/Baum-Welch ###

Given that we have access to observation sequences, and possibly the structure of a PFSA, we can estimate its parameters (transition probabilities) with the EM algorithm (aka Baum-Welch or forward-backward).  To train a PFSA with Baum-Welch, we can either

  * provide an initial PFSA in a file
  * initialize one randomly with various options

### Using an existing PFSA ###

If we can provide an initial PFSA and a set of observations we can re-estimate its transition probabilities by

```
treba --train=bw --file myfsm.fsm myobs.obs > learnedpfsa.fsm
```

which will run Baum-Welch until the difference in log likelihood between iterations reaches the default threshold delta of 0.1 (tunable with the `--max-delta` option) and output the results to the standard output, in this example redirected to the file `learnedpfsa.fsm`.

### Random initialization ###

We can also generate the initial PFSA randomly for training using the `--initialize` flag and specifying the topology of the automaton and its size. For example:

  * `--initialize=4`   would generate a fully connected (ergodic), and hence nondeterministic automaton with 4 states (see image below)
  * `--initialize=b4`  would generate a Bakis (left-to-right), automaton with 4 states (see image below)
  * `--initialize=d4`  would generate an initial deterministic automaton with 4 states (mostly used for testing purposes)

The automaton alphabet is inferred from the observation sequence. That is, if a nondeterministic automaton is specified, it will contain a transition from each state to each other state with every symbol in the alphabet (with some initial random probability).

<img src='http://treba.googlecode.com/svn/wiki/automatamodels.png' alt='Logo' width='400' />

Note that the alphabet size is always determined automatically from the observation sequences.  If another alphabet size is required (perhaps larger than what the largest occurring symbol in the observations would indicate), it can be specified after the number of states: e.g. `--initialize=4,10` would force an alphabet size of 10.

### Parallel threads in Baum-Welch ###

If we have access to multiple CPUs, we can use the `--threads` flag to specify the number of threads to launch in Baum-Welch training.  For example:

```
treba --train=bw --threads=8 --initialize=n40 myobs.obs > learnedpfsa.fsm
```
would launch 8 threads for calculating Baum-Welch.

### Termination ###

Baum-Welch can sometimes take a long time to converge, especially if the delta (`--max-delta`) parameter is set low. To prevent losing valuable calculations, hitting **CTRL-C** in the middle of Baum-Welch will always output the results of the last iteration, and then exit.

## Baum-Welch for HMMs ##

Baum-Welch for HMMs works the same way, except the machine output is an hmm.  This is specified, again, with the `--hmm` flag.  For example:

```
treba --hmm --train=bw --initialize=10 myobs.obs > learnedhmm.hmm
```

will train a 10-state fully connected HMM from the observations in `myobs.obs`.


## Viterbi training ##

For very large HMMs/PFSAs and a very large number of observation sequences, Viterbi training (aka hard EM) may be a useful alternative to Baum-Welch.  It converges much faster than Baum-Welch since it only uses the most probable path (hence the name Viterbi) to estimate the transition probabilities.  It has an additional optional parameter (`--prior`) used to add pseudocounts of transitions between states to prevent producing zero-probability paths in the result. The default value of `--prior` is 1.0.  To specify a transition and emission pseudocount separately, one can issue a comma-separated value for the prior parameter, e.g. `--prior=0.5,1.0`.

To train with Viterbi instead of Baum-Welch, one would issue, for example:

```
treba --hmm --train=vit --initialize=40 myobs.obs > learnedhmm.hmm
```

to run Viterbi training with an initial random HMM of 40 states using the observations in `myobs.obs`.


## Variational Bayes ##

Variational Bayes training may also be specified with the option `--train=vb`.  The algorithm is the same core algorithm as for EM/Baum-Welch, however, a prior may be specified with the `--prior` option.  For example:

```
treba --train=vb --prior=0.5 --initialize=n40 myobs.obs > learnedpfsa.fsm
```

would train a 40-state randomly initialized PFSA with Variational Bayes.  Variational Bayes tends to converge faster than EM; however, the PFSA it returns are not normalized for each state, i.e. the transitions may sum to less than 1.

Again, the prior for emissions and transitions may be specified separately in the HMM case, by `--prior=t,e` where `t` is the transition prior and `e` the emission prior.

# Gibbs sampling #

Another alternative to estimating the probabilities of a PFSA/HMM through Baum-Welch, is to use a Gibbs sampler to sample a posterior distribution, and then reconstruct a PFSA/HMM from the samples.  Treba includes a collapsed Gibbs sampler that does this.  Inferring a PFSA as in the above Baum-Welch example, but using a Gibbs sampler instead, can be done as:

```
treba --train=gs --initialize=n40 myobs.obs > learnedpfsa.fsm
```

The main parameters to be tune for the Gibbs sampler are `--burnin`, the number of iterations to run before starting to collect samples (the burn-in), `--max-iterations`, the total number of iterations to run, and `--lag`, how many iterations to skip during sample collection to avoid correlated samples.  Also, the Dirichlet prior may be set by `--prior`. The default lag is 10, the default burn-in 100, and the default prior is 0.02.

For example:

```
treba --train=gs --burnin=1000 --lag=100 --max-iterations=20000 --prior=0.05 --initialize=n40 myobs.obs > learnedpfsa.fsm
```

Would run a Gibbs sampler for 20000 iterations, only collecting samples from every 100th iteration after the first 1000 iterations, using a prior of 0.05.

## Gibbs sampling for HMMs ##

HMM (instead of PFSA) inference is also supported by the Gibbs sampler. To infer the probabilities of a 40-state HMM using the same data as in the above example could be done by adding the `--hmm`-flag, i.e:

```
treba --hmm --train=gs --burnin=1000 --lag=100 --max-iterations=4000 --prior=0.05,0.01 --initialize=n40 myobs.obs > learnedhmm.hmm
```

When training HMMs, we assume two priors: one prior for the state-to-state transitions (0.05 in the example), and one for the symbols emissions (0.01 in the example), which are specified comma-separated.  Defaults are 0.02 (transition) and 0.01 (emission).


## Gibbs sampling with CUDA ##

To run the Gibbs sampler for HMMs/PFSA using a GPU instead of CPU with NVIDIA CUDA, use the `--cuda` flag.  Note that this requires that you have compiled support for CUDA (or downloaded a precompiled version with CUDA support and have the CUDA runtime), and have a supported NVIDIA graphics card.

For example:

```
treba --cuda --hmm --train=gs --initialize=n40 myobs.obs > learnedhmm.hmm
```


# State-merging inference algorithms #

State-merging inference algorithms produce deterministic probabilistic finite automata as output given a sequence of observations as input.  Their sizes are not predefined; rather, the algorithms attempt to produce the smallest possible PFSA by merging states. The algorithms initialize a trie-shaped structure (a deterministic finite frequency automaton), after which states are merged whenever they are deemed compatible.  The compatibility test and its parameters are variable, and altering them yields different generalizations.

## ALERGIA ##

The standard state-merging algorithm is `ALERGIA`.  It is run as default if state-merging algorithms are called:

```
treba --train=merge myobs.obs > learnedpfsa.fsm
```

All state-merging algorithms rely on a parameter (α) for determining when to merge states, which is set through `--alpha` .  For `ALERGIA`, normal values are between 0.01 and 1, the default is 0.05.  Also, `ALERGIA` can be supplied a `--t0`-parameter, indicating a minimum number of attested strings passing a state in the prefix tree to consider merging.  The default value is 3. A high value will prevent many potential merges.  For example:

```
treba --train=merge --alpha=0.5 --t0=10 myobs.obs > learnedpfsa.fsm
```

will run `ALERGIA` with α=0.5 and `t0` being 10.

Several alternative tests can be used instead of the default `ALERGIA` (Hoeffding bound) test for merging.  A simple test that sometimes works well is the _likelihood ratio_ test for statistical significance `--merge-test=lr`:

```
treba --train=merge --merge-test=lr myobs.obs > learnedpfsa.fsm
```

Here, the α-value corresponds to the minimum desired p-value for statistical significance of a test concerning the distribution of the arcs in a state before and after merging.  With this, and other available statistical tests, the t0-value is ignored since the statistical test should be equally representative independently of the number of attested times a transition has been traversed in the prefix automaton.

## MDI ##

`MDI` is another state-merging algorithm, though it is not grouped with the other state-merging algorithms for technical reasons.  The main parameter to set is again `--alpha`, e.g.:

```
treba --train=mdi --alpha=0.005 myobs.obs > learnedpfsa.fsm
```

## Smoothing in state-merging ##

Sometimes the PFSA resulting from state merging operations does not contain transitions from certain states with certain symbols.  This may be a problem for some data sets over the same alphabet, as some strings will be deemed infinitely improbable. There is a option to smooth the transitions and finality of states when creating the PFSA after induction; this is controlled with the `--prior` flag, e.g.:

```
treba --train=merge --prior=0.1 myobs.obs > learnedpfsa.fsm
```

would add a pseudocount of 0.1 to each missing transition and each final state before reconstructing the PFSA from the frequency automaton.

Doing so will guarantee a probability to any observation over the same alphabet as the training set.

This value can be set to `0` for no smoothing, e.g.

```
treba --train=merge --prior=0.0 myobs.obs > learnedpfsa.fsm
```

## Generation ##

Treba can also be used to generate sequences according to the parameters of a HMM/PFSA.  To generate 100 sequences from `myfsm.fsm`, one would launch `treba` with:

```
treba --generate-words=100 --file myfsm.fsm
```

which prints out, in three TAB-separated columns: the **probability**, the **observation**, and the **state sequence** generated at random, e.g.:

```
3.0240000000000012e-05	2 0 2 6 5	0 0 1 2 2 3
3.3339600000000063e-05	0 3 2 3 3 5	0 1 1 1 2 2 3
...
```

Similarly, for a HMM

```
treba --hmm --generate-words=100 --file myhmm.hmm
```

would print out 100 words at random in the above three-column format.


# Miscellaneous options #

## Baum-Welch restarts ##

Since Baum-Welch training is very sensitive to the initial parameters, there exists an option to run the algorithm repeatedly with different random initializations (option `--restarts`) for a specified number of iterations each, and then choosing the one with the highest log likelihood to continue until convergence.  For example

```
treba --hmm --train=bw --initialize=20 --restarts=10,100
```

would run Baum-Welch on a 20-state HMM for 100 iterations 10 times, and then continuing until convergence with the best one found.

## Deterministic annealing ##

Treba also supports augmenting Baum-Welch with deterministic annealing (see references).  In effect, this will weight the re-estimation procedure in such a way that in the early stages, paths are weighted toward more uniform probabilities, while repeatedly running Baum-Welch until convergence.  Each time Baum-Welch converges, the `beta` parameter which weights the observation counts is raised by a factor of `alpha` (usually until it reaches 1.0, which is equivalent to running Baum-Welch in standard fashion).  The idea of the maneuver is to escape being stuck in local optima to which Baum-Welch is susceptible. For example,

```
treba --train=dabw --initialize=20 myobs.obs > learnedpfsa.fsm
```

would run a deterministically annealed Baum-Welch with a random 20-state PFSA using the default parameters (initial beta = 0.02, maximum beta = 1.0, alpha = 1.01).  These parameters can be changed with the `--annealing-params` flag (e.g. `--annealing-params=0.02,1.0,1.2` would change the alpha to 1.2).

Note that deterministic annealing requires a very long time to run compared with Baum-Welch, and there are no real guarantees of it yielding a better model than straight-up Baum-Welch.  See the literature for details.

## Input/output formats ##

By default, `treba` uses real numbers for probabilities (though internally log2 calculations are performed).  The `--input-format` and the `--output-format` flags can be used to control the format of the input and the output.  For example, `--output-format=log10`, would print all output probabilities in log10.  Negative logprobs can be specified with an n-prefix: `--output-format=nln` would, for example, output negative log probabilities with base **e**. Negative logprobs are often used to convert probabilities to weights, and doing so allows easy export of PFSA files to other tools.

## Uniform probabilities ##

The `--uniform-probs` flag can be used in conjunction with the `--initialize` flag for generating initial uniform probability automata (instead of the default random ones).  If this is used in conjunction with a restart scheme, only the first automaton generated will have uniform probabilities; subsequent automata will be random.  Note that initializing HMMs/PFSAs uniformly for Baum-Welch is not advisable, since this is often a local optimum the algorithm cannot escape.

# Bibliography #

Baum, L. E., T. Petrie, G. Soules, and N. Weiss. (1970). A maximization technique occurring in the statistical analysis of probabilistic functions of Markov chains. _Annals of Mathematical Statistics_, vol. 41, no. 1, pages. 164–171.

Carrasco, R. C., & Oncina, J. (1994). Learning stochastic regular grammars by means of a state merging method. In Grammatical Inference and Applications (pp. 139-152). Springer.

de la Higuera, C. (2010). _Grammatical Inference: Learning Automata and Grammars_.  Cambridge University Press.

Dempster, A., N. Laird, and D. Rubin. (1977). Maximum likelihood estimation from incomplete data via the EM algorithm. _Journal of the Royal Statistical Society B_, 39:1-38.

Hulden, M. (2012). Treba: Efficient Numerically Stable EM for PFA. Journal of Machine Learning Research Workshop and Conference Proceedings, Vol. 21: ICGI 2012: 249-253.

MacKay, D. J. (1997). Ensemble learning for hidden Markov models. Technical report, Cavendish Laboratory, University of Cambridge

Rabiner, L. R. (1989). A tutorial on Hidden Markov Models and selected applications in speech recognition.  _Proc. of the IEEE_, 77(2):257-286.

Rao,  A.  and  K. Rose. (2001). Deterministically annealed design of Hidden Markov Model speech recognizers. _IEEE Transactions on Speech and Audio Processing_, 9(2):111-126.

Rose, K. (1998). Deterministic annealing for clustering, compression, classification, regression, and  related  optimization  problems.  _Proc.  of  the  IEEE_, 86(11):2210-2239.

Shibata, C., & Yoshinaka, R. (2012). Marginalizing Out Transition Probabilities for Several Subclasses of PFAs. Journal of Machine Learning Research Workshop and Conference Proceedings, Vol. 21: ICGI 2012: 259-263.

Thollard, F., P. Dupont, and C. de la Higuera. (2000). Probabilistic DFA inference using Kullback-Leibler divergence and minimality. In _Proceedings of the 17th International Conference on Machine Learning_:975–982. Morgan Kaufmann: San Francisco.

Ueda, N. and Nakano, R. (1998). Deterministic annealing EM algorithm. _Neural Networks_, 11(2):271-282.

Verwer, S., Eyraud, R., and de la Higuera, C. (2013). PAutomaC: a PFA/HMM learning competition. _Machine Learning_. Springer.

Smith, N. A. and Eisner, J. (2004). Annealing techniques for unsupervised statistical language learning. In _Proc. of the ACL_, pages 486-493.