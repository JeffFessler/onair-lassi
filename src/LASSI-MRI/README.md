# LASSI on MRI data

The `*_run.m` files in this folder run LASSI on MRI data.

These files were designed to be used in the context of parameter sweeps.
The included `*_par*.m` files show how to define range(s) of parameter(s) that
you want to try:

```
invivo_LASSI_full_par4.m
invivo_LASSI_R8_par4.m
LASSI_full_par4c.m
LASSI_R8_par10.m
```

You can see the available range of `idx` values like so:

```
invivo_LASSI_full_run(invivo_LASSI_full_par4);
invivo_LASSI_R8_run(invivo_LASSI_R8_par4);
LASSI_full_run(LASSI_full_par4c);
LASSI_R8_run(LASSI_R8_par10);
```

You can run an experiment with a particular set of hyperparameters by passing
the appropriate `idx`:

```
invivo_LASSI_full_run(invivo_LASSI_full_par4, idx);
invivo_LASSI_R8_run(invivo_LASSI_R8_par4, idx);
LASSI_full_run(LASSI_full_par4c, idx);
LASSI_R8_run(LASSI_R8_par10, idx);
```

The experiment outputs are written to `*.mat` files as configured in the
`*_par*.m` files.

## Dependencies

The input data can be downloaded from
https://drive.google.com/file/d/158cnwBAnKqFgX6yrcADme9wBHMNWtGfj/view?usp=sharing.

The source files for the `deps_lassi/` directory are located in the `src/`
folder of this repository.
