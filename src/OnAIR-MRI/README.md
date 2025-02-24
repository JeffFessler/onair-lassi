# OnAIR on MRI data

The `onlineDls_run.m` file in this folder runs OnAIR on MRI data.

This file was designed to be used in the context of parameter sweeps.
The included `*_par*.m` files show how to define range(s) of parameter(s) that
you want to try:

```
pincat_onlineDls_par3g.m
otazo_onlineDls_par3g.m
invivo_onlineDls_par3g.m
```

You can see the available range of `idx` values like so:

```
onlineDls_run(pincat_onlineDls_par3g);
onlineDls_run(otazo_onlineDls_par3g);
onlineDls_run(invivo_onlineDls_par3g);
```

You can run an experiment with a particular set of hyperparameters by passing
the appropriate `idx`:

```
onlineDls_run(pincat_onlineDls_par3g, idx);
onlineDls_run(otazo_onlineDls_par3g, idx);
onlineDls_run(invivo_onlineDls_par3g, idx);
```

The experiment outputs are written to `*.mat` files as configured in the
`*_par*.m` files.

## Dependencies

The input data can be downloaded from
https://drive.google.com/file/d/158cnwBAnKqFgX6yrcADme9wBHMNWtGfj/view?usp=sharing.

The source files for the `deps_dls/` directory are located in the `src/` folder
of this repository.
