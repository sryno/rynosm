Installation:

1) Anaconda 3.6 or later must be installed and in PATH
2) gmake and gcc or similar must be installed
3) cd to source directory
4) run `make all`
5) Run program from any directory via `python /path/to/getConnectionsMolecules.py -h`
	e.g., python ~/files/getConnectionsMolecules.py test-system.gro test-out -conn -mol

Note:
	If the compiler gives an error that omp.h is missing please comment out the line #include <omp.h> in c_nm_conn_molID.c
