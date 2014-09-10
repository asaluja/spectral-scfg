`spectral-scfg` is a package for latent synchronous context-free grammars (L-SCFGs) as applied to machine translation (MT). It is a direct implementation of the following paper:

> Latent-variable Synchronous CFGs for Hierarchical Translation.  Avneesh Saluja, Chris Dyer, and Shay B. Cohen. In *Proceedings of EMNLP*, 2014. 

The source code for the implementation of the minimal grammar extractor used in this paper (Zhang et al., COLING 2008) is
also in this repository, and was generously provided by Dan Gildea.  It is also available [here](http://www.cs.rochester.edu/u/gildea/mt/factorize-alignment.tgz).  

## System Requirements
- NumPy (at least v1.6) and SciPy (at least v0.13)
- a relatively recent vereion of `cdec` ([link](https://github.com/redpony/cdec) which includes `pycdec`.
- MATLAB for fast SVD computation

In addition, a valid C++ compiler is needed to compile the grammar extractor.  

## Prerequisites

- Compiled version of `cdec` along with compiled Python extensions (see [here](http://www.cdec-decoder.org/guide/compiling.html) for instructions).  Make sure the location of the `pycdec` files is in your `$PYTHONPATH` environment variable. 
- Compiled version of the minimal grammar extractor (see instructions within `ZGC-extractor` sub-directory for more information)

## End-to-end Pipeline Instructions

Input: parallel sentence corpus (tokenized, lower-cased). Assume all code is run from the repository root. 

0. Preprocess the sentence corpus to escape special characters (e.g., '[' and ']').

   ```
   python scripts/escape_special_chars.py < parallel_corpus_in > filtered_corpus_out
   ```

   Do the same with the development and test sets.  

1. Use `fast_align` in `cdec` along with `atools` to get symmetrized bi-directional word alignments. Also, compile the parallel data into a suffix array, which we can use to get lexical scores when featurizing our grammars. 

   ```
   $CDEC/word-aligner/fast_align -i filtered_corpus_out -d -v -o > fwd_align
   $CDEC/word-aligner/fast_align -i filtered_corpus_out -d -v -o -r > rev_align
   $CDEC/utils/atools -i fwd_align -j rev_align -c grow-diag-final-and > align.gdfa
   python -m cdec.sa.compile -b filtered_corpus_out -align.gdfa -c extract.ini -o training.sa
   ```

2. Preprocess word alignments to format that minimal rule extractor can take in.

   ```
   python minrule-extraction/preproc_minrule_input.py < align.gdfa > minrule.input
   ```

3. Run minimal rule extraction.

   ```
   ZGC-extractor/rc < minrule.input > minrule.log 2> minrule.output
   ```

4. Convert output of minimal rule extraction to per-sentence grammars.  These grammars can be in one of two formats, either a) a "Hiero" variant containing only lexicalized rules (with maximum arity 2) or b) a "full" variant, which encodes the entire minimal derivation tree.  Both are needed. 

   ```
   python minrule-extraction/tree_to_rule.py -z minrule.hiero.grammar/ filtered_corpus_out align.gdfa < minrule.output &> invalid.rules
   python minrule-extraction/tree_to_rule.py -d -z minrule.full.grammar filtered_corpus_out align.gdfa < minrule.output
   ```

5. From the "full" per-sentence grammars, create a counts file, that simply estimates count(e,f) (source RHS, target RHS) for minimal rules.

   ```
   python parameter-estimation/count_estimation.py minrule.full.grammar/ counts-file
   ```

   Note: by adding the `-n` flag, you can also estimate the counts from the Hiero-style rules; the argument to the `-n` flag needs to be the list of invalid rules from the previous conversion. 

6. As an additional step to make the grammar compatible with `cdec`, you need to featurize the minimal rules and write it all out in one large grammar. 

   ```
   python minrule-extraction/featurize_rules.py minrule.hiero.grammar/ minimal.featurized.grammar.gz counts-file training.sa/lex.bin numProcesses
   ```

   Flags:
   - `-a`: add one to counts when writing out features
   - `-s`: output a per-sentence grammar instead of a combined grammar
   - `-f X`: filter rules to top X sorted by P(e|f)
   At this stage, you can use the minimal grammar (and refer to it it appropriately in a `cdec` decoder configuration file) to decode. 

7. Parameter estimation can be done in one of three ways:

  1. SVD-based estimation:

   ```
   python parameter-estimation/svd_estimation.py minrule.entire.grammar/ feature_list rank parameters > training.rules.with.rowIdxs
   ```

   Flags (rule indicator features always used):
     - `-o`: estimate OOV probabilities (recommended)
     - `-f X`: filter rules to top X sorted by P(e|f) (recommended)
     - `-a`: arity feature
     - `-l`: lexical feature
     - `-c X`: class-based lexical feature (word IDs are replaced by class IDs); classes for both languages should be provided in the format `-c source-file:target-file`
     - `-L`: span length feature
     - `-r X`: real-valued features (not fully supported)
     - `-s X`: smoothing of parameters, with hyperparameter value X

  2. MLE estimation:

   ```
   python parameter-estimation/svd_estimation.py -m minrule.entire.grammar/ dummy 1 mle.parameters > training.rules.with.rowIdxs
   ```

   Flags: 
     - `-o`: estimate OOV probabilities (recommended)
     - `-f X`: filter rules to top X sorted by P(e|f) (recommended)

  3. EM estimation:

   ```
   parameter-estimation/em_estimation.py minrule.entire.grammar/ rank numIterations outDirForParameters scaling
   ```

   The scaling argument is used so that we don't underflow, and should be set to something reasonably high (e.g., 10^5 or 10^6). Flags:
     - `-o`: whether to do OOV estimation or not; like svd_estimation.py, OOV estimation here is based on singletons (recommended)
     - `-f X`: filter rules, need to provide as an argument a dictionary of parameters where the rules are already filtered
     - `-m X`: Matsuzaki-style initialization: need to provide a dictionary of MLE parameters
     - `-n X`: number of cores to use; default is 8

8. To decode an evaluation set (make sure it has also been preprocessed to escape special characters), run the parser with the estimated parameters.

   ```
   python parser/compute_hg.py -f parameters rank test_sentences_escaped numProcesses outDirForPerSentenceGrammars
   ```

   Flags:
	- `-m`: required for MLE estimation (i.e. if parameters are MLE parameters)
	- `-s`: write out rule marginal normalized by source RHS as well
	- `-t`: write out rule marginal normalized by target RHS as well
	- `-d`: for debugging (traverses the hypergraph and prints out nodes)
	- `-n`: write-out node-marginals
	- `-x`: only write marginals for non-lexical rules on source-side

9. After decoding, the per-sentence grammars need to be decorated with other common MT features for effective performance in an MT system. 

   ```
   python minrule-extraction/featurize_rules.py -m outDirForPerSentenceGrammars/ outDirForDecoratedGrammars/ counts-file training.sa/lex.bin numProcesses
   ```

   Flags:
	- `-a`: add one to counts

The final output is a directory of per-sentence grammars, each rule in the grammar decorated with the rule marginal probabilities. 
For downstream processing in MT, the following command/script may be useful to convert the evaluation set to a .sgm file where each sentence points to the location of the per-sentence grammar.

```
python scripts/corpus2sgm.py outDirForDecoratedGrammars < test_sentences_escaped > test_sent.sgm
```

## Things to add

- add SciPy SVD computation as an option to those who do not have access to MATLAB
- option to output inside and outside probability vectors as features in addition or instead of the rule marginal features
- finish adding support for real-valued features

## Citation

If you make use of this package, please cite:
Latent-variable Synchronous CFGs for Hierarchical Translation.  Avneesh Saluja, Chris Dyer, and Shay B. Cohen. In *Proceedings of EMNLP*, 2014. [[bibtex](http://www.cs.cmu.edu/~avneesh/EMNLP_2014.pdf)][[pdf](http://www.cs.cmu.edu/~avneesh/EMNLP_2014.pdf)]
