# Compatibility

This version of the estimator uses [Lattigo v6.1.0](https://github.com/tuneinsight/lattigo).

# How to run the experiments

Each experiment runs the same circuit with both estimator and Lattigo in on the same input and then prints statistic about the resulting precision of the output.

1. `$cd estimator/experiments/<experiment_name>/`
2. Optional: open `<file>.go`, set the parameters and number of iterations to the desired values
3. `$go run <file>.go` to run the experiment

# Citing

Please use the following BibTex entry for citing the estimator:

    @misc{cryptoeprint:2024/853,
      author = {Jean-Philippe Bossuat and Anamaria Costache and Christian Mouchet and Lea Nürnberger and Juan Ramón Troncoso-Pastoriza},
      title = {Practical q-{IND}-{CPA}-D-Secure Approximate Homomorphic Encryption},
      howpublished = {Cryptology {ePrint} Archive, Paper 2024/853},
      year = {2024},
      url = {https://eprint.iacr.org/2024/853}
    }
    
