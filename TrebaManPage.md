<pre>
treba(1)                                                              treba(1)<br>
<br>
<br>
<br>
NAME<br>
treba - probabilistic finite-state automaton training and decoding<br>
<br>
SYNOPSIS<br>
treba [options] [ OBSERVATIONS-FILE ]<br>
<br>
DESCRIPTION<br>
treba  is  a tool for training, decoding, and calculating with weighted<br>
(probabilistic) finite state automata and hidden Markov models  (HMMs).<br>
Training  algorithms  include  Baum-Welch (EM), Viterbi training, Gibbs<br>
sampling, merge-based algorithms (ALERGIA, MDI et al.), and  Baum-Welch<br>
augmented with deterministic annealing.  Training algorithms can be run<br>
multi-threaded (EM/hard EM/Variational Bayes) on hardware with multiple<br>
cores/CPUs,  or  with  NVIDIA GPUs (Gibbs sampler).  Forward, backward,<br>
and Viterbi decoding are supported.  Automata/HMMs for  training/decod-<br>
ing  are  read  from  a text file, or can be generated randomly or with<br>
uniform transition probabilities with different topologies (ergodic  or<br>
fully  connected,  Bakis or left-to-right, or deterministic).  Observa-<br>
tions used for training or decoding are  read  from  text  files.   The<br>
resulting  automata, path calculations, or probability calculations are<br>
printed to standard output.<br>
<br>
<br>
BASIC OPTIONS<br>
The program has four main modes of  operation:  training  (the  --train<br>
option),  decoding  (the  --decode option), likelihood calculation (the<br>
--likelihood option), and sequence generation (the --generate  option).<br>
All modes except generation depend on an observations file and possibly<br>
a finite-state automaton or HMM file (if not initialized  randomly  for<br>
training).  The observations file is a text file and is assumed to con-<br>
sist of whitespace-separated sequences  of  integers,  one  observation<br>
sequence  on each line.  Empty lines correspond to the empty string and<br>
are acceptable.  All output is directed to the  standard  output:  when<br>
run  in  training  mode, the output will consist of a FSA/HMM in a text<br>
format (see below); when run in decoding mode, the  output  will  be  a<br>
string  of whitespace-separated integers representing the most probable<br>
path through a given FSA or HMM, one path for  each  observation;  when<br>
run  in likelihood calculation mode, the output will be one probability<br>
for each line in the observations file.<br>
<br>
<br>
OBSERVATION FILES<br>
Observations files are assumed to contain one observation on each line,<br>
separated  by  whitespace.  Each observation is a sequence of integers,<br>
representing symbols.  It is assumed that the lowest numbered symbol is<br>
0, and that the "alphabet" contains no gaps.  The following illustrates<br>
an observations file with three observations:<br>
<br>
<br>
1 7 4 3 6 3 9 7 8 10 11 1 1<br>
3 4 3 2<br>
5 5 5 1 0 0 1 0<br>
<br>
Lines beginning with the #-symbols are assumed to be comments  and  are<br>
ignored in all inputs (observations or automata/HMM specifications).<br>
<br>
<br>
DETAILED OPTIONS<br>
--hmm  Use  HMMs  instead  of  the default input/output which is proba-<br>
bilistic finite automata<br>
<br>
--cuda Enable NVIDIA-CUDA  for  Gibbs  sampling  (if  compiled  in  and<br>
NVIDIA-card is present)<br>
<br>
<br>
--train=merge|mdi|bw|dabw|gs|vb|vit|vitbw<br>
Train  a  model  with  one  of the algorithms merge (ALERGIA and<br>
variants), mdi (MDI), bw  (Baum-Welch),  dabw  (Baum-Welch  with<br>
deterministic  annealing),  gs  (Gibbs  sampler), vb Variational<br>
Bayes (EM), vit (Viterbi training/hard EM),  or  vitbw  (Viterbi<br>
training followed by Baum-Welch).<br>
<br>
--decode=f|b|vit[,p]<br>
Decode  (find  the best path) through the automaton/HMM for each<br>
word in obervation-file using either the forward path (  f  )  ,<br>
the  backward  path  (  b  ),  or the Viterbi path ( vit ).  The<br>
optional ,p will  also  print  out  the  respective  probability<br>
together with the path.  Note that forward and backward decoding<br>
chooses the most probable state for each point in  time  and  so<br>
the  path may or may not correspond to an actually valid path in<br>
the automaton.  For example, --decode=vit,p will  calculate  the<br>
Viterbi path and print its probability.<br>
<br>
--likelihood=f|vit<br>
Calculate  the  likelihood (probability) for each observation in<br>
observation-file using forward probability, or the Viterbi prob-<br>
ability.<br>
<br>
--generate=NUM<br>
Generate  NUM  random  sequences  from  FSA/HMM.   Randomness is<br>
weighted by transition probabilities.  The sequences are  output<br>
in three TAB-separated fields: (1) the sequence probability; (2)<br>
the symbol sequence itself; (3) the state sequence.<br>
<br>
<br>
--input-format=FORMAT<br>
Set format of probabilities in input automata (real  numbers  or<br>
logs or negative logs in various bases), one of real, log10, ln,<br>
log2, nlog10, nln, nlog2.  Default is real.<br>
<br>
--output-format=FORMAT<br>
Set format  of  probabilities  related  to  output  automata  or<br>
results of decoding and likelihood calculations (real numbers or<br>
logs or negative logs in various bases), one of real, log10, ln,<br>
log2, nlog10, nln, nlog2.  Default is real.<br>
<br>
--file=FILENAME<br>
Specify  finite  state  automaton or HMM file.  Each line in the<br>
automaton file consists of one to four  numbers:  a  four-number<br>
line  S1 S2 A P indicates a transition from S1 to S2 with symbol<br>
A and probability P whereas a line of the format S  P  indicates<br>
the  final  probability at state S.  Three-number and one-number<br>
lines are identical to the above, with an implicit probabibility<br>
P of one.  The initial state is always state number 0.  The for-<br>
mat is identical to the text formats accepted by  the  AT&T  FSM<br>
toolkit  or  OpenFST,  with  the  exception that strings are not<br>
allowed to represent symbols or states: all symbols  and  states<br>
need to be integers.<br>
<br>
The  following  snippet  illustrates  a  typical FSA file of two<br>
states with an alphabet size of three using  real-valued  proba-<br>
bilities:<br>
<br>
<br>
0 0 0 0.25<br>
0 0 1 0.25<br>
0 1 0 0.2<br>
0 1 1 0.1<br>
0 1 2 0.1<br>
1 0 0 0.15<br>
1 0 1 0.15<br>
1 1 0 0.3<br>
1 1 1 0.1<br>
1 1 2 0.1<br>
0 0.1<br>
1 0.2<br>
<br>
<br>
The following illustrates a HMM:<br>
<br>
<br>
0 > 1 0.44<br>
0 > 2 0.03<br>
0 > 3 0.53<br>
1 > 1 0.11<br>
1 > 2 0.08<br>
1 > 3 0.81<br>
2 > 1 0.48<br>
2 > 2 0.25<br>
2 > 3 0.27<br>
1 0 0.21<br>
1 1 0.79<br>
2 0 0.92<br>
2 1 0.08<br>
<br>
<br>
That  is,  transitions  in  HMMs are specified with lines of the format<br>
SOURCE-STATE > TARGET-STATE TRANSITION-PROBABILITY, and emissions  with<br>
lines of the format STATE SYMBOL EMISSION-PROBABILITY.<br>
<br>
<br>
--initialize=TYPE-NUMSTATES[,NUMSYMBOLS]<br>
Generate  an  initial  automaton of type type for training.  The<br>
type is a combination of an option letter and the number of  the<br>
states and the alphabet of the format [b|d]#states(,#numsymbols)<br>
Here, b will generate a Bakis  (left-to-right)  automaton,  d  a<br>
deterministic  automaton,  and the default with no specification<br>
generates a fully connected (ergodic)  nondeterministic  automa-<br>
ton.  If numsymbols is omitted, the necessary alphabet size will<br>
be deduced from the observation file.  For  example,  --initial-<br>
ize=20 will generate an initial fully connected random automaton<br>
with 20 states, inferring the necessary alphabet size  from  the<br>
data  while  --initialize=20,10  will  generate an initial fully<br>
connected random automaton with 20 states and force the alphabet<br>
to  be of size 10 symbols whereas --initialize=b10 will generate<br>
a Bakis automaton with 10  states,  and  the  alphabet  size  is<br>
deduced from the largest symbol in the observations file.<br>
<br>
--uniform<br>
This  option  sets  all  initially  generated  automata with the<br>
--initialize option to have  uniform  probabilities  instead  of<br>
randomly assigned ones.<br>
<br>
--help print help and exit<br>
<br>
--version<br>
print version and exit<br>
<br>
<br>
TRAINING OPTIONS<br>
--max-delta=MAXDELTA<br>
Maximum  delta (change in log likelihood between training itera-<br>
tions) to use for convergence.  Default is 0.1.  Note that treba<br>
will output the result of training calculations so far if CTRL-C<br>
is pressed without the need to wait for convergence.<br>
<br>
--max-iter=NUM<br>
Maximum number of iterations in  training.   Default  is  20000.<br>
Note  that treba will output the result of training calculations<br>
so far if CTRL-C is pressed without the need to wait for conver-<br>
gence.<br>
<br>
--prior=PRIOR1[,PRIOR2]<br>
Priors  to use in various contexts: the Dirichlet prior in Gibbs<br>
sampling, a smoothing prior in merging  algorithms,  or  pseudo-<br>
counts  to use in Viterbi training to prevent paths of probabil-<br>
ity zero. For HMMs, two priors are usually specified, where  the<br>
first one is state-to-state transition prior, and the second one<br>
is the emission prior,  e.g.  --prior=0.1,0.2  --merge=ALGORITHM<br>
Choose   merging   algorithm,   one   of   alergia,chi2,lr,bino-<br>
mial,exactm,exact (ALERGIA default Hoeffding bound,  chi-squared<br>
test,  likelihood  ratio, binomial test, exact multinomial test,<br>
exact test).  --alpha=VALUE This is the main parameter for state<br>
merging  algorithms. Default is 0.05.  --t0=VALUE The t0-parame-<br>
ter for ALERGIA/MDI.  --recursive-merge Enables  recursive  com-<br>
patibility checks for state merging algorithms. Default: OFF.<br>
<br>
--restarts=OPT<br>
where OPT=numrestarts,iterations-per-restart. Sets Baum-Welch to<br>
restart itself restarts times running for iterations-per-restart<br>
iterations  each  restart.   After all random restarts have been<br>
run, the fsm with the best log likelihood is  chosen  and  Baum-<br>
Welch proceeds as normal.<br>
<br>
<br>
--annealopts=betamin,betamax,alpha<br>
Controls  the  parameters  for  deterministic annealing when run<br>
with -T dabw setting the initial beta value to betamin, the max-<br>
imum  beta  value  to  betamax and the multiplier alpha by which<br>
beta in increased each time Baum-Welch converges.   The  default<br>
values are 0.02, 1.0, and 1.01.<br>
<br>
--threads=NUM<br>
Number  of  threads to launch in Baum-Welch training.  The value<br>
num-threads can be optionally prefixed by c or  c/  to  use  c-n<br>
threads or c/n threads, where c is the number of logical CPUs on<br>
the system.  To use half  the  available  processors  one  would<br>
issue  the  flag  --threads=c/2 and likewise --threads=c1 to use<br>
all but one of the available CPUs/cores. Default value is 1.<br>
<br>
EXAMPLE USAGE<br>
treba --train=bw --initialize=10 sentences.txt<br>
Reads all the sentences from sentences.txt and trains a 10-state<br>
probabilistic  automaton  using  Baum-Welch  using  the  default<br>
training parameters.  The initial automaton has random probabil-<br>
ities on its transitions and fully connected.<br>
<br>
treba --hmm --train=bw --initialize=25 sentences.txt<br>
Reads all the sentences from sentences.txt and trains a 25-state<br>
HMM using Baum-Welch using the default training parameters.  The<br>
initial  automaton  has  random probabilities on its transitions<br>
and fully connected.<br>
<br>
treba --hmm --train=gs --initialize=n20 sentences.txt<br>
Reads all the sentences from sentences.txt and trains a 20-state<br>
HMM with Gibbs sampling using the default parameters of burn-in,<br>
lag, and maximum iterations.<br>
<br>
treba --cuda --train=gs --initialize=n40 --burnin=1000 --lag=100 --max-<br>
iter=10000 sentences.txt<br>
Reads all the sentences from sentences.txt and trains a 40-state<br>
PFSA  with Gibbs sampling using a burn-in of 1000, and then col-<br>
lecting samples each 100 iterations,  running  for  a  total  of<br>
10000  iterations.  The --cuda option forces the GPU implementa-<br>
tion to be used (only available if compiled  in  and  an  NVIDIA<br>
card is present).<br>
<br>
treba --train=bw --threads=c/2 --initialize=b100 sentences.txt<br>
Reads   all  the  sentences  from  sentences.txt  and  trains  a<br>
100-state probabilistic automaton using Baum-Welch.  The initial<br>
automaton is left-to-right.  During training, half the available<br>
CPUs will be used.<br>
<br>
treba --train=merge --alpha=0.01 sentences.txt<br>
Reads all the sentences from sentences.txt and infers an  deter-<br>
ministic  probabilistic finite automaton using the ALERGIA algo-<br>
rithm, with the alpha parameter set to 0.01.<br>
<br>
treba --train=dabw --threads=8 --file=initial.fsm sentences.txt<br>
Reads all the sentences from sentences.txt and trains  a  proba-<br>
bilistic  automaton  using Baum-Welch with deterministic anneal-<br>
ing.  The initial automaton is read from initial.fsm.   A  total<br>
of 8 threads will be launched in parallel for Baum-Welch.<br>
<br>
treba --train=vit --initialize=b25,5 --maxiter=10 sentences.txt<br>
Reads all the sentences from sentences.txt and trains a 25-state<br>
probabilistic automaton with an alphabet size of 5 using Viterbi<br>
training running a maximum of 10 iterations.  The initial autom-<br>
aton is random and left-to-right (Bakis).<br>
<br>
treba --likelihood=f --file=myfsm.fsm sentences.txt<br>
Reads sentences.txt and calculates for each observation line the<br>
forward probability in the automaton in myfsm.fsm.<br>
<br>
treba --decode=vit --file=myfsm.fsm sentences.txt<br>
Reads sentences.txt and calculates for each observation line the<br>
most  probable  path  (the  Viterbi  path)  in   the   automaton<br>
myfsm.fsm.<br>
<br>
treba --decode=vit,p --file=myfsm.fsm sentences.txt<br>
Reads sentences.txt and calculates for each observation line the<br>
most probable path through  the  automaton  in  myfsm.fsm.   The<br>
probability of the path is also printed.<br>
<br>
treba --generate=100 --output-format=log10 --file=myfsm.fsm<br>
Generate 100 random sequences (weighted by transition probabili-<br>
ties) from myfsm.fsm.  Output probability scores are  log10  and<br>
myfsm.fsm has real-valued transitions (default).<br>
<br>
treba --input-format=real --output-format=nln --file=myfsm.fsm<br>
No  real  action: myfsm.fsm is read (inputs are real-valued) and<br>
converted to negative logprobs (aka weights) with the base e and<br>
output  to  stdout.   This can be used to export an FSA for pro-<br>
cessing with e.g. the AT&T tools or OpenFST.  Not issuing any of<br>
the  flags --train, --decode or --likelihood, will simply output<br>
the input FSA, performing a possible conversion depending on the<br>
input and output specifiers.<br>
<br>
<br>
SEE ALSO<br>
http://treba.googlecode.com<br>
<br>
Baum, L. E., T. Petrie, G. Soules, and N. Weiss. (1970). A maximization<br>
technique occurring in the statistical analysis of probabilistic  func-<br>
tions of Markov chains. Annals of Mathematical Statistics, vol. 41, no.<br>
1, pages. 164-171.<br>
<br>
Carrasco, R. C., & Oncina, J. (1994). Learning stochastic regular gram-<br>
mars  by  means of a state merging method. In Grammatical Inference and<br>
Applications (pp. 139--152). Springer.<br>
<br>
de la Higuera, C. (2010). Grammatical Inference: Learning Automata  and<br>
Grammars.  Cambridge University Press.<br>
<br>
Dempster,  A., N. Laird, and D. Rubin. (1977). Maximum likelihood esti-<br>
mation from incomplete data via the EM algorithm. Journal of the  Royal<br>
Statistical Society B, 39:1-38.<br>
<br>
Hulden,  M.  (2012).  Treba:  Efficient  Numerically Stable EM for PFA.<br>
Journal of Machine Learning Research Workshop and  Conference  Proceed-<br>
ings, Vol. 21: ICGI 2012: 249--253.<br>
<br>
MacKay, D. J. (1997). Ensemble learning for hidden Markov models. Tech-<br>
nical report, Cavendish Laboratory, University of Cambridge.<br>
<br>
Rabiner, L. R. (1989). A tutorial on Hidden Markov Models and  selected<br>
applications  in speech recognition.  Proc. of the IEEE, 77(2):257-286.<br>
<br>
Rao, A. and K. Rose. (2001). Deterministically annealed design of  Hid-<br>
den  Markov  Model  speech recognizers. IEEE Transactions on Speech and<br>
Audio Processing, 9(2):111-126.<br>
<br>
Rose, K. (1998). Deterministic annealing for  clustering,  compression,<br>
classification, regression, and related optimization problems. Proc. of<br>
the IEEE, 86(11):2210-2239.<br>
<br>
Shibata, C., & Yoshinaka, R. (2012). Marginalizing Out Transition Prob-<br>
abilities  for  Several Subclasses of PFAs. Journal of Machine Learning<br>
Research Workshop and  Conference  Proceedings,  Vol.  21:  ICGI  2012:<br>
259--263.<br>
<br>
Thollard,  F.,  P.  Dupont, and C. de la Higuera. (2000). Probabilistic<br>
DFA inference using Kullback-Leibler divergence and minimality. In Pro-<br>
ceedings  of  the  17th  International  Conference  on  Machine  Learn-<br>
ing:975-a982. Morgan Kaufmann: San Francisco.<br>
<br>
Ueda, N. and Nakano, R. (1998). Deterministic annealing  EM  algorithm.<br>
Neural Networks, 11(2):271-282.<br>
<br>
Verwer,  S.,  Eyraud,  R.,  and  de  la Higuera, C. (2012). Pautomac: a<br>
PFA/HMM learning competition. In Proceedings of the 11th  International<br>
Conference on Grammatical Inference.<br>
<br>
Smith,  N.  A. and Eisner, J. (2004). Annealing techniques for unsuper-<br>
vised statistical  language  learning.  In  Proc.  of  the  ACL,  pages<br>
486-493.<br>
<br>
fsm(1)<br>
<br>
<br>
BUGS<br>
Many and hairy. This is academic code. Don't build Mars rovers with it.<br>
<br>
<br>
AUTHOR<br>
Mans Hulden <mans.hulden@gmail.com><br>
<br>
<br>
<br>
November 27, 2013                      treba(1)<br>
</pre>


&lt;hr&gt;

