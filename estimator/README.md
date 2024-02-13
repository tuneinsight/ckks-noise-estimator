# Dependencies

The estimator requires a modified version of [Lattigo v5.0.2](https://github.com/tuneinsight/lattigo). The small modifications were notably done to improve how the precision statistics are computed and enable access to private fields.

The modified version of Lattigo v5.0.2 is located in the folder `./lattigo`.

# How to run the experiments

Each experiment runs the same circuit with both estimator and Lattigo in on the same input and then prints statistic about the resulting precision of the output.

1. `$cd experiments`
2. Optional: open `<file>.go`, set the parameters and number of iterations to the desired values
3. `$go run <file>.go` to run the experiment