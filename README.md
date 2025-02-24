## OnAIR and LASSI code

This repo contains code for algorithms
OnAIR, LASSI, SOUP-DIL,
and some other variations on them (like OptShrink-powered variations).

### OnAIR

OnAIR for video denoising: `src/demo_onlineDls.m`

OnAIR for MRI: `src/OnAIR-MRI/`

### LASSI

LASSI for video denoising: `src/demo_drpca.m`

LASSI for MRI: `src/LASSI-MRI/`

### Other methods

There is a `README` file in the `src/` directory
that gives starting points for other files.

### Notes

To keep this repo small,
the data files (`.png` and `.mat` files)
are stored separately at
https://web.eecs.umich.edu/~fessler/irt/reproduce/20/moore-20-oai/

You will likely need to copy that data
into the `src/` directory
to run any of the scripts.

### Citation

The code here is related to the following papers:

* Brian E Moore, Saiprasad Ravishankar, Raj Rao Nadakuditi, J A Fessler,
"Online adaptive image reconstruction using (OnAIR) using dictionary models".
IEEE Trans. on Computational Imaging 6:153-66, 2020.
http://doi.org/10.1109/TCI.2019.2931092

* Saiprasad Ravishankar, Brian E Moore, Raj Rao Nadakuditi, J A Fessler,
"Low-rank and adaptive sparse signal (LASSI) models
for highly accelerated dynamic imaging".
IEEE Trans. on Medial Imaging 36(5):1116-28, May 2017.
http://doi.org/10.1109/TMI.2017.2650960
(See also https://gitlab.com/ravsa19/lassi)

Please cite the above paper(s) if you use this code.
