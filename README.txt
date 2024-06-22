Numerical continuation for finding Turing base state solutions to heterogeneous RD equations.

directory hierarchy:

src: All source code, including LaTeX code for thesis and Python code for plots.

    continuation: All the code for the numerics, excluding pde2path.
        common: this has any code which is relevant to multiple different codes.
                For example, all of the reaction terms are defined in here.
            gm: Gierer-Meinhardt reaction with saturating term
            gm_orig: Gierer-Meinhardt reaction without saturating term
            schnakenberg: Schnakenberg reaction term. Untested.
            schnakenberg_krause: The Schnakenberg reaction term used by the Krause2020 paper.

    lib: external libraries. pde2path, ilupack...
        pde2path: Numerical Continuation Library.
        ilupack: Not used, but an incomplete lu linear equation solver.


The continuation directory has multiple subdirectories for different numerical experiments.

**IMPORTANT**
To run any experimental results -
First change the matlab directory to the one you want to run.
Then run 'setpath'. This will initialise pde2path and the environment.
Then you can run the script that you want.

domainLength:
    This is for the experiments for the Schnakenberg critical domain length.
    scanCriticalEps.m runs the code to find the derivative of the critical domain length with respect to the size of the heterogeneity. These are the results presented in the thesis.
    scanCriticalEps2.m runs the code to find the derivative of the critical domain length for a range of theta and beta0 values. This was presented in the paper.

domainLength2_gm:
    This is for the experiments for the Gierer-Meinhardt critical domain length.
    scanCritGam.m runs the code to find the derivative of the critical domain length with respect to the size of the heterogeneity. These are the results presented in the thesis.
    scanCritGam2.m runs the code to find the derivative of the critical domain length for a range of theta and a0 values. This was presented in the paper.

folds:
    This is for the experiments for the Schnakenberg base state existence.
    foldSize.m runs the code to find the base state existence while changing the production rate beta0 and gamma.
    foldSize2.m runs the code to find the base state existence while changing the form of the heterogeneity and gamma.
    qualitative_results.m runs the code to evaluate the case studies listed in the thesis/paper.

folds_gm_orig:
    This is for the experiments for the Gierer-Meinhardt base state existence.
    foldSize.m runs the code to find the base state existence while changing the production rate beta0 and gamma.
    foldSize2.m runs the code to find the base state existence while changing the form of the heterogeneity and gamma.

krauseSteadyState:
    Reproducing the steady states from Krause2020.
    range_eps and range_krause are the ones which run the code.


