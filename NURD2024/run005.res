Run Start Date:  17/05/2023 
Run Start Time: 12:01:03.43 
Run Stop Date:   17/05/2023 
Run Stop Time:  12:01:22.96 
  
************************************************************************************************************************ 
********************                            CONTROL STREAM                                      ******************** 
************************************************************************************************************************ 
  
$PROB Warfarin PKPD; run005

;O'Reilly RA, Aggeler PM. Studies on coumarin anticoagulant drugs
;Initiation of warfarin therapy without a loading dose.
;Circulation 1968;38:169-177
;
;O'Reilly RA, Aggeler PM, Leong LS. Studies of the coumarin anticoagulant
;drugs: The pharmacodynamics of warfarin in man.
;Journal of Clinical Investigation 1963;42(10):1542-1551


$INPUT ID TIME AMT DV EVID WT AGE SEX SPARSE
$DATA warfarin_PKS.csv IGNORE=@

$SUBR ADVAN13 TOL=9
$MODEL
   COMP=GUT
   COMP=CENTRAL
$PK
	POP_ka  = THETA(1)
	POP_cl  = THETA(2)
	POP_v   = THETA(3)

    ka = exp(POP_ka+ETA(1))
    cl = exp(POP_cl+ETA(2))
    v = exp(POP_v+ETA(3))
    
$DES
    gut   = A(1)
    center= A(2)

    DADT(1) = -ka * gut
    DADT(2) =  ka * gut - cl / v * center

$ERROR
    CONC   = A(2)/v
	IPRED  = CONC
    RESCV  = THETA(4) 
	RESADD = THETA(5)
    W      = SQRT(IPRED*IPRED*RESCV*RESCV+RESADD*RESADD)	
    IRES   = DV-IPRED
    IWRES  = IRES/W
    Y      = IPRED + W*EPS(1)

$THETA 0.15 ; log ka
$THETA  -2  ; log cl
$THETA   2  ; log v
$THETA 0.2  ; EPS_prop
$THETA 0.5  ; EPS_pkadd

$OMEGA
1   ; ETA_ka
1   ; ETA_cl
1   ; ETA_v

$SIGMA   1 FIX


$EST METHOD=1 INTERACTION MAXEVAL=9999 SIGL=5 PRINT=10 NOHABORT
$COV UNCONDITIONAL 
$TABLE ID TIME AMT DV EVID WT AGE SEX SPARSE Y ka cl v IRES IWRES IPRED ETAS(1:LAST)
ONEHEADER NOPRINT FILE=run005.csv
   
   
************************************************************************************************************************ 
********************                            SYSTEM INFORMATION                                  ******************** 
************************************************************************************************************************ 
   
Operating system: 

Microsoft Windows [Version 10.0.22621.1702]
   
Compiler: Intel(R) Parallel Studio XE 2015 Update 3 Composer Edition (package 208)  
   
Compiler settings: /Gs /nologo /nbs /w /fp:strict   
   
************************************************************************************************************************ 
********************                            NMTRAN MESSAGES                                     ******************** 
************************************************************************************************************************ 
   
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.
   
************************************************************************************************************************ 
********************                            NONMEM EXECUTION                                    ******************** 
************************************************************************************************************************ 
   
License Registered to: Occams Cooperatie UA
Expiration Date:    14 JUN 2024
Current Date:       17 MAY 2023
Days until program expires : 392
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 Warfarin PKPD; run005
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      283
 NO. OF DATA ITEMS IN DATA SET:  10
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   4
 MDV DATA ITEM IS DATA ITEM NO.: 10
0INDICES PASSED TO SUBROUTINE PRED:
   5   2   3   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME AMT DV EVID WT AGE SEX SPARSE MDV
0(NONBLANK) LABELS FOR PRED-DEFINED ITEMS:
 ka cl v IPRED IRES IWRES Y
0FORMAT FOR DATA:
 (9E5.0,1F2.0)

 TOT. NO. OF OBS RECS:      251
 TOT. NO. OF INDIVIDUALS:       32
0LENGTH OF THETA:   5
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   3
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
   0.1500E+00 -0.2000E+01  0.2000E+01  0.2000E+00  0.5000E+00
0INITIAL ESTIMATE OF OMEGA:
 0.1000E+01
 0.0000E+00   0.1000E+01
 0.0000E+00   0.0000E+00   0.1000E+01
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
0TABLES STEP OMITTED:    NO
 NO. OF TABLES:           1
 SEED NUMBER (SEED):    11456
 RANMETHOD:             3U
 MC SAMPLES (ESAMPLE):    300
 WRES SQUARE ROOT TYPE (WRESCHOL): EIGENVALUE
0-- TABLE   1 --
0RECORDS ONLY:    ALL
04 COLUMNS APPENDED:    YES
 PRINTED:                NO
 HEADER:                YES
 FILE TO BE FORWARDED:   NO
 FORMAT:                S1PE11.4
 IDFORMAT:
 LFORMAT:
 RFORMAT:
 FIXED_EFFECT_ETAS:
0USER-CHOSEN ITEMS:
 ID TIME AMT DV EVID WT AGE SEX SPARSE Y ka cl v IRES IWRES IPRED ETA1 ETA2 ETA3
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 GENERAL NONLINEAR KINETICS MODEL WITH STIFF/NONSTIFF EQUATIONS (LSODA, ADVAN13)
0MODEL SUBROUTINE USER-SUPPLIED - ID NO. 9999
0MAXIMUM NO. OF BASIC PK PARAMETERS:   3
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         GUT          ON         YES        YES        YES        NO
    2         CENTRAL      ON         YES        YES        NO         YES
    3         OUTPUT       OFF        YES        NO         NO         NO
 INITIAL (BASE) TOLERANCE SETTINGS:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   9
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            *           *           *           *           *
    3            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      5
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   3

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
0ERROR SUBROUTINE INDICATES THAT DERIVATIVES OF COMPARTMENT AMOUNTS ARE USED.
0DES SUBROUTINE USES COMPACT STORAGE MODE.
1
 
 
 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction
 
 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            9999
 NO. OF SIG. FIGURES REQUIRED:            3
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      5
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     5
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): run005.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE

 TOLERANCES FOR ESTIMATION/EVALUATION STEP:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   9
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
 TOLERANCES FOR COVARIANCE STEP:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   9
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
 TOLERANCES FOR TABLE/SCATTER STEP:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   9
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
 
 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI
 
 MONITORING OF SEARCH:

 
0ITERATION NO.:    0    OBJECTIVE VALUE:   558.484996855942        NO. OF FUNC. EVALS.:   7
 CUMULATIVE NO. OF FUNC. EVALS.:        7
 NPARAMETR:  1.5000E-01 -2.0000E+00  2.0000E+00  2.0000E-01  5.0000E-01  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01 -1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   1.9152E+01  1.4744E+02 -3.9857E+01  8.9320E+02 -2.1284E+01 -6.7409E+00  5.7892E+01  5.7793E+01
 
0ITERATION NO.:   10    OBJECTIVE VALUE:   429.590743212616        NO. OF FUNC. EVALS.:  89
 CUMULATIVE NO. OF FUNC. EVALS.:       96
 NPARAMETR: -6.5735E-01 -2.0124E+00  2.1386E+00  1.2454E-01  7.2135E-01  3.7464E-01  4.9359E-02  5.6451E-02
 PARAMETER: -4.3824E-01 -1.0062E-01  1.0693E-01  6.2271E-02  1.4427E-01 -3.9089E-01 -1.4043E+00 -1.3372E+00
 GRADIENT:  -3.1658E+00  2.6029E+01  1.4748E+03  1.0327E+02  3.3994E+01 -3.0748E+00 -2.5798E+01  3.6696E-01
 
0ITERATION NO.:   20    OBJECTIVE VALUE:   423.327278202764        NO. OF FUNC. EVALS.: 153
 CUMULATIVE NO. OF FUNC. EVALS.:      249
 NPARAMETR: -5.6879E-01 -2.0072E+00  2.0585E+00  1.1704E-01  7.2249E-01  4.7490E-01  7.9029E-02  4.8275E-02
 PARAMETER: -3.7920E-01 -1.0036E-01  1.0293E-01  5.8522E-02  1.4450E-01 -2.7232E-01 -1.1690E+00 -1.4154E+00
 GRADIENT:  -2.0655E-03 -7.3082E-01  1.2216E-01  7.6988E-02  3.1440E-03  1.8932E-03 -9.1520E-03  2.6215E-03
 
 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      249
 NO. OF SIG. DIGITS IN FINAL EST.:  4.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.
 
 ETABAR:        -3.2106E-02 -2.1099E-05 -1.0301E-02
 SE:             6.4953E-02  4.7081E-02  3.3117E-02
 N:                      32          32          32
 
 P VAL.:         6.2110E-01  9.9964E-01  7.5575E-01
 
 ETASHRINKSD(%)  4.5829E+01  3.7444E+00  1.3373E+01
 ETASHRINKVR(%)  7.0655E+01  7.3487E+00  2.4957E+01
 EBVSHRINKSD(%)  4.4421E+01  5.0670E+00  1.3986E+01
 EBVSHRINKVR(%)  6.9110E+01  9.8773E+00  2.6016E+01
 RELATIVEINF(%)  3.0657E+01  8.9418E+01  7.3964E+01
 EPSSHRINKSD(%)  1.3028E+01
 EPSSHRINKVR(%)  2.4360E+01
 
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          251
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    461.307143668746     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:    423.327278202764     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       884.634421871510     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                            96
  
 #TERE:
 Elapsed estimation  time in seconds:     6.48
 Elapsed covariance  time in seconds:     3.08
 Elapsed postprocess time in seconds:     0.05
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************      423.327       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5     
 
        -5.69E-01 -2.01E+00  2.06E+00  1.17E-01  7.22E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        4.75E-01
 
 ETA2
+        0.00E+00  7.90E-02
 
 ETA3
+        0.00E+00  0.00E+00  4.83E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        6.89E-01
 
 ETA2
+        0.00E+00  2.81E-01
 
 ETA3
+        0.00E+00  0.00E+00  2.20E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5     
 
         2.46E-01  5.38E-02  4.60E-02  2.65E-02  1.54E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        2.24E-01
 
 ETA2
+       .........  2.40E-02
 
 ETA3
+       ......... .........  1.44E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        1.63E-01
 
 ETA2
+       .........  4.26E-02
 
 ETA3
+       ......... .........  3.29E-02
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      OM11      OM12      OM13      OM22      OM23      OM33      SG11  

 
 TH 1
+        6.05E-02
 
 TH 2
+       -2.94E-03  2.90E-03
 
 TH 3
+        2.74E-04  7.83E-04  2.12E-03
 
 TH 4
+       -2.33E-03  4.56E-04  1.85E-04  7.04E-04
 
 TH 5
+       -1.58E-02  1.71E-04 -3.87E-03 -6.25E-04  2.37E-02
 
 OM11
+        2.41E-02  1.20E-03  5.13E-04  9.63E-04 -1.07E-02  5.03E-02
 
 OM12
+       ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+        2.59E-04  3.74E-04 -2.32E-04 -1.36E-04 -4.62E-04  9.92E-04 ......... .........  5.74E-04
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+        2.24E-05 -4.24E-05  9.67E-05 -6.28E-05  2.19E-04 -3.15E-04 ......... .........  8.99E-05 .........  2.08E-04
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      OM11      OM12      OM13      OM22      OM23      OM33      SG11  

 
 TH 1
+        2.46E-01
 
 TH 2
+       -2.22E-01  5.38E-02
 
 TH 3
+        2.42E-02  3.16E-01  4.60E-02
 
 TH 4
+       -3.57E-01  3.19E-01  1.51E-01  2.65E-02
 
 TH 5
+       -4.16E-01  2.06E-02 -5.47E-01 -1.53E-01  1.54E-01
 
 OM11
+        4.37E-01  9.90E-02  4.98E-02  1.62E-01 -3.10E-01  2.24E-01
 
 OM12
+       ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+        4.40E-02  2.90E-01 -2.11E-01 -2.14E-01 -1.25E-01  1.85E-01 ......... .........  2.40E-02
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+        6.30E-03 -5.46E-02  1.46E-01 -1.64E-01  9.86E-02 -9.73E-02 ......... .........  2.60E-01 .........  1.44E-02
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      OM11      OM12      OM13      OM22      OM23      OM33      SG11  

 
 TH 1
+        4.66E+01
 
 TH 2
+       -7.36E+01  1.04E+03
 
 TH 3
+        1.58E+02 -1.15E+03  2.28E+03
 
 TH 4
+        2.79E+02 -1.17E+03  1.59E+03  3.83E+03
 
 TH 5
+        6.27E+01 -3.15E+02  5.88E+02  5.68E+02  2.14E+02
 
 OM11
+       -2.12E+01  2.57E+01 -5.35E+01 -1.51E+02 -1.61E+01  3.32E+01
 
 OM12
+       ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+        2.80E+02 -1.94E+03  2.98E+03  3.23E+03  8.76E+02 -1.54E+02 ......... .........  6.58E+03
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       -2.28E+02  1.61E+03 -2.82E+03 -2.07E+03 -8.01E+02  1.20E+02 ......... ......... -4.83E+03 .........  8.94E+03
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.04
 #CPUT: Total CPU Time in Seconds,        2.734
