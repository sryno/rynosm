This is an omega tuning script based on the one originally written by Cheng Zhong in the Bredas group at KAUST. 

I have rewritten the program to provide additional information due to errors, allow compatibility with Python 
2.7 and 3.6, and allow any version of Gaussian to be used, if the version supports the necessary options. 

If you have any questions do not hesitate to ask.

-Sean


This program requires either the Anaconda Python distribution to be installed or the following packages:
	- numpy
	- scipy

All Gaussian input files must have a gjf or com extension.

`python aw_tuning.py -h` will return a list of possible options for running the script.

Example input: `python aw_tuning input.gjf > input.log`
This input will run aw_tuning for gap tuning with the default optional keywords. However, these have shown
to be a good combination for difficult to converge systems.

input.gjf:
%chk=ppv.chk
%nprocshared=16
%mem=8GB
#p lc-wpbe/may-cc-pvdz			(define the range separated functional)

aw_tuning example

0 1
C	34.084552	1.193596	-0.981740
C	32.704045	1.025247	-0.856374
...



###########################################################################################################

usage: aw_tuning.py [-h] [-g GAUSS] [-l LEVEL] [-m MEMORY] [-n NPROC]
                    [-p PRECISION] [-P ALPHAPRECISION] [-b BOUND]
                    [-B ALPHABOUND] [-k KEYWORDS] [-c CRITERION] [-t TUNING]
                    [-a ALPHA] [-w OMEGA]
                    [gaussianInput [gaussianInput ...]]

A script used with Gaussian09/16 to auto optimize omega and alpha for long-range corrected DFT methodologies. This program would not be possible without the original version by Cheng Zhong.

positional arguments:
  gaussianInput

optional arguments:
  -h, --help            show this help message and exit
  -g GAUSS, --gauss GAUSS
                        Set the executable name to be used. Default: g09
  -l LEVEL, --level LEVEL
                        Set the theoretical level. Default: read from gjf file. If not found, use LC-wPBE/6-31g(d)
  -m MEMORY, --memory MEMORY
                        Set the memory used by Gaussian (with unit).  Default: read from gjf file. If not found, use: 500MB
  -n NPROC, --nproc NPROC
                        Set the number of cores used by Gaussian. Default: read from gjf file. If not found, use 1
  -p PRECISION, --precision PRECISION
                        Set the precision of omega. Default: 4, means omega will converge to 0.0001
  -P ALPHAPRECISION, --alphaPrecision ALPHAPRECISION
                        Set the precision of alpha. Default: equal to precision of omega
  -b BOUND, --bound BOUND
                        Set the boundary for omega, e.g. -b 0.05,0.3. Default: 0.05,0.5
  -B ALPHABOUND, --alphaBound ALPHABOUND
                        Set the boundary for alpha, e.g. -B 0.05,0.3. Default: 0.0,0.5
  -k KEYWORDS, --keywords KEYWORDS
                        Set additional Gaussian keywords. Default: "scf=(xqc,fermi,noincfock,ndamp=35,conver=6,vshift=500,novaracc)"
  -c CRITERION, --criterion CRITERION
                        Set the optimization criterion (the value to minimize), available options are:
                        J2---((HOMO-IP)^2+(A_HOMO+EA)^2)
                        Jh---(HOMO-IP)
                        Jl(LUMO+EA)
                        Jn2---((HOMO-IP)^2+(LUMO+EA)^2)
                        O2---((A_HOMO-LUMO)^2+(C_LUMO-HOMO)^2)
                        or any other valid experssion, e.g. -c "(A_HOMO-N_LUMO)**2"
                        If multiple criterion is used, seperate them by comma. The order and number of criterions must be the same to that of structures. Default: J2
  -t TUNING, --tuning TUNING
                        Set the tuning scheme. Available options are: w(omega),a(alpha),aw(alpha and omega),s(scan alpha or omege, must be used with -a or -w option)
                        Put m in front of above options (mw,ma,maw,ms) means tuning with multiple structure,
                        e.g. "-t mw" means optimize the omega to minimize sum of the criterion of all the input structure.
                        If different criterion is used for different structure, use -c option to designate criterion explicitly for each structure,
                        e.g. "-t ma -c Jh,Jl" means optimize alpha to minimize sum of Jh of structure 1 and Jl of structure 2.  Default: w
  -a ALPHA, --alpha ALPHA
                        Set the value of alpha. Can be a single value or multiple value seperated by space, e.g. -a '0.1 0.3 0.5 0.7' Default: None
  -w OMEGA, --omega OMEGA
                        Set the value of omega. Can be a single value or multiple value seperated by space. Default: None