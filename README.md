# EigenKernel

## About
EigenKernel is a package of hybrid parallel solvers for (standard / generalized) real-symmetric eigenvalue problems.
Here 'hybrid solver' means combination of high-performance routines from the distributed parallel numerical linear algebra libraries, ScaLAPACK, ELPA and EigenExa. You can test various combination of routines from a unique interface and find the best one for your problem and machine environment.

EigenKernel also works as a mini-application (stand alone eigenvalue problem solver), named EigenKernel_App. This note shows usage of EigenKernel_App shortly.

EigenKernel_App requires input matrix data in [the MatrixMarket file format](http://math.nist.gov/MatrixMarket/). Sample data files are stored in the `matrix/` directory. These are part of [ELSES matrix library](http://www.elses.jp/matrix/).

## Quick start
You can build EigenKernel_App only with ScaLAPACK (without ELPA and EigenExa).
Test following commands to solve a generalized eigenvalue problem with the matrix size of M = 30. The sample Makefile.inc supposes there are mpif90 over gfortran as a fortran compiler and libgomp as an OpenMP library.

    tar zxvf eigenkernel-*.tar.gz
    cd eigenkernel-*
    cp Makefile.inc.gfortran.noext Makefile.inc
    make

The mini-application appears as bin/eigenkernel_app. You can run the application, for example, as

    mpirun -np 4 bin/eigenkernel_app -s general_scalapack matrix/ELSES_MATRIX_BNZ30_A.mtx matrix/ELSES_MATRIX_BNZ30_B.mtx

and obtain the output files named `eigenvalues.dat`, `ipratios.dat`, `log.json`. eigenvalues.dat and ipratios.dat contain computed eigenvalues and inversed participation ratios of computed eigenvectors respectively. log.json contains execution information such as given commandline options in the JSON format.

When you need eigenvectors,
you should specify the index range to be output by `-p` option. For example,

    mpirun -np 4 bin/eigenkernel_app -s general_scalapack -d vector/ -p 1-30 matrix/ELSES_MATRIX_BNZ30_A.mtx matrix/ELSES_MATRIX_BNZ30_B.mtx

We should note that the directory `vector/` must be created before execution.

In the execution command `-s <solver>` is a mandatory option to specify the solver routine. The general_scalapack solver is, of course, a pure ScaLAPACK solver. The last two arguments are paths to input matrix files in the Matrix Market format. Note that the input matrices must be symmetric and moreover the latter one must be positive definite (only real-valued matrices are supported now).

You can also solve standard eigenvalue problems.

    mpirun -np 4 bin/eigenkernel_app -s scalapack matrix/ELSES_MATRIX_VCNT400std_A.mtx


## Link with ELPA and EigenExa
You must build and link ELPA and EigenExa to utilize full functions of EigenKernel_App. Follow installation instructions in each library.

[ELPA official page](http://elpa.rzg.mpg.de/)

[EigenExa official page](http://www.aics.riken.jp/labs/lpnctrt/en/projects/eigenexa/)

After installing them, edit Makefile.inc (sample found in Makefile.inc.gfortran.withext) and set `LIBS` variable properly to indicate the paths where .a and .mod are installed.
Then you should rebuild EigenKernel_App with ELPA and EigenExa like below.

    emacs Makefile.inc  # Edit $LIBS properly
    make clean
    make WITH_EIGENEXA=1 WITH_ELPA=1

If the rebuild is succeeded, now you can select truly hybrid solvers in the `-s` option.

    mpirun -np 4 bin/eigenkernel_app -s general_elpa_eigensx matrix/ELSES_MATRIX_VCNT900_A.mtx matrix/ELSES_MATRIX_VCNT900_B.mtx

In log.json (this default filename can be changed with `-l <filename>` option) you can find whole computation time in an element of the array `events` whose `name` is `main`. More detailed time consumption for each subprocedure is also reported. You can compare the solvers with these information.


## Available solvers
Names of available solvers are listed below. `general_scalapack` and `general_elpa_eigensx` are already used above. 'standard' solvers are for standard eigenvalue problems and take one matrix file for the last argument. 'generalized' solvers are for generalized ones and take two matrix files for the last arguments. The 'selecting' solvers can compute a part of eigenvectors. The number of computed eigenvectors can be specified with the `-n <num>` option.

### pure ScaLAPACK solvers
- scalapack (standard)  -- PDSYEVD
- scalapack_select (standard, selecting) -- PDSYEVX
- general_scalapack (generalized) -- reduction with PDPOTRF & PDSYGST; The solver 'A' in our papers [1,2]
- general_scalapack_select (generalized, selecting)
- general_scalapacknew (generalized) -- reduction with PDPOTRF & PDSYNGST

### solvers need ELPA
- general_elpa_scalapack (generalized) -- reduction with ELPA, solve SEP with PDSYEVD; The solver 'C' in our papers [1,2]
- general_elpa1 (generalized); The solver 'E' in our papers [1,2]
- general_elpa2 (generalized); The solver 'D' in our papers [1,2]

### solvers need EigenExa
- eigensx (standard)
- general_scalapack_eigensx (generalized) -- reduction with PDPOTRF & PDSYGST, solve SEP with eigen_sx; The solver 'B' in our papers [1,2]
- general_scalapack_eigens (generalized) -- reduction with PDPOTRF & PDSYGST, solve SEP with eigen_s
- general_scalapacknew_eigensx (generalized) -- reduction with PDPOTRF & PDSYNGST, solve SEP with eigen_sx
- general_scalapacknew_eigens (generalized) -- reduction with PDPOTRF & PDSYNGST, solve SEP with eigen_s

### solvers need both of ELPA and EigenExa
- general_elpa_eigensx (generalized) -- reduction with ELPA, solve SEP with eigen_sx; The solver 'G' in our papers [1,2]
- general_elpa_eigens (generalized) -- reduction with ELPA, solve SEP with eigen_s; The solver 'F' in our papers [1,2]


## Useful commandline options (output eigenvectors, check accuracy, change default output filename)
You can see a help message for commandline options by `eigenkernel_app -h`.

- `-n <num>`  This option can be used only with 'selecting' solvers. Compute only &lt;num&gt; eigenpairs in ascending order of their eigenvalues.
- `-p <num1>[-<num2>][,<num3>[-<num4>]]...`  Specify index ranges of eigenvectors to be output. No eigenvectors are output in default. Eigenvectors are output from all the MPI processes evenly. The maximum number of index ranges is 100.
- `-c <num>`  Calculate residual norm of eigenvectors whose index is from 1 to &lt;num&gt; (included). If num == -1, EigenKernel_App calculates all the eigenvectors. No eigenvectors are checked in default.
- `-t <num1>,<num2>`  Calculate orthogonality of eigenvectors whose index is from &lt;num1&gt; to &lt;num2&gt; (included). No eigenvectors are checked in default.
- `-o <file>`  Set output file name for eigenvalues to &lt;file&gt;.
- `-i <file>`  Set output file name for ipratios to &lt;file&gt;.
- `-d <dir>`  Set output files directory for eigenvectors to &lt;dir&gt;.
- `-l <file>`  Set output file name for elapse time log to &lt;file&gt;
- `--block-size <n>`  Change block size in the block cyclic distribution. (default:64)
- `--dry-run`  Exit before starting eigensolver. Used for testing matrix read and broadcast.


## Supported ELPA and EigenExa versions
The current master branch (v．2018) supports
- ELPA: elpa-2018.05.001
- EigenExa: EigenExa 2.4b

EigenKernel v.2017 supports
- ELPA: elpa-2014.06.001
- EigenExa: EigenExa 2.3c



## Reference
[1] H. Imachi and T. Hoshi, '[Hybrid numerical solvers for massively parallel eigenvalue computation and their benchmark with electronic structure calculations](https://www.jstage.jst.go.jp/article/ipsjjip/24/1/24_164/_article)', J. Inf. Process. 24, pp. 164 -- 172 (2016).

[2] K.Tanaka, H. Imachi, T. Fukumoto, T. Fukaya, Y. Yamamoto, T. Hoshi, 'EigenKernel - A middleware for parallel generalized eigenvalue solvers to attain high scalability and usability', submitted; Preprint: http://arxiv.org/abs/1806.00741

If you use the present software in an academic publication, please cite the above reference.


## License
This software is released under the MIT License. See the file LICENSE.
