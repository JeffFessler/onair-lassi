## OnAIR and LASSI code

This repo contains code for algorithms
OnAIR, LASSI, SOUP-DIL,
and some other variations on them (like OptShrink-powered variations).

The code here is related to the following papers

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

In the `src/` directory,
the script that runs OnAIR is called
`demo_onlineDls.m`.
Note that this is only for video denoising.
We have not yet had a chance to put together encapsulated scripts
that use sensing matrices for MRI, for example.

There is a `README` file in the `src/` directory
that gives starting points for other files.

To keep this repo small,
the data files (`.png` and `.mat` files)
are stored separately at
https://web.eecs.umich.edu/~fessler/irt/reproduce/20/moore-20-oai/

You will likely need to copy that data
into the `src/` directory
to run any of the scripts.

We hope to provide more complete code and data in the future.

Please cite the above paper(s) if you use this code.


For further code and data:
* Here's a link to a 160MB `.zip` download that (should) contain a standalone set of files. There are `demo_XXX` files that run demos for various algorithms.
https://drive.google.com/file/d/1m-swPi6jLsfoaLeGesYdSQaCmw234luq

This folder contains OnAIR, LASSI, SOUP-DIL, and some other variations on them (like OptShrink-powered variations).
We haven't gone through the files and rename things or separate out into per-paper code.  That `.zip` file was the starting point for this repo.

The script that runs OnAIR is called `demo_onlineDls.m` and is only for video denoising.  We haven't had a chance to put together encapsulated scripts that use sensing matrices for MRI, for example.

* We also have a 13GB dump of everything done by the first author, including scripts and probably both input and output data for experiments used to generate figures in our various LASSI/OnAIR conference/journal papers.
https://drive.google.com/file/d/12a3pHIQbJxVVuWQEnv1ic42fYzGFvZml
