10/04/2024 
09:54
$PROB  BRV Adult popPK model ; run025 analysis simulated data

$DATA DAT_A.csv IGNORE=@ 
$INPUT ID RND OCC TIME MDV AMT SS II DV WT AED

$SUB ADVAN2
$PK

;Allometrically scaled to 70kg WT patient
TVCL   = THETA(2)
TVVC   = THETA(3)
TVKA   = THETA(4)
ALLOCL = THETA(5)
ALLOVC = THETA(6)
AEDCL  = THETA(7)

CL    = EXP(LOG(TVCL) + ALLOCL*(LOG(WT/70)) + AEDCL*AED + ETA(1)) ;Clearance (L/hr)
VC    = EXP(LOG(TVVC) + ALLOVC*(LOG(WT/70)) + ETA(2)) ;Volume (L)
KA    = EXP(LOG(TVKA) + ETA(3)) ;Ka (1/hr) 

S2    = VC
K     = CL/VC

$ERROR
IPRED = F      
RESCV = THETA(1)
IRES  = DV-IPRED
W     = IPRED*RESCV
IWRES = (DV-IPRED)/W 
Y     = IPRED+W*EPS(1)

$THETA  (0.001,0.22);Res.error (SD/mean) 
$THETA  (0.001,3.7) ;CL (L/hr)
$THETA  (0.01,50)   ;V (L)
$THETA  (0.01,1.5)  ;Ka (1/hr)
$THETA  0.75 FIX    ;WT on Cl
$THETA  1    FIX    ;WT on V
$THETA  0.2         ;AED effect on CL

$OMEGA .1 .3 1.3

$SIGMA 1 FIXED 

$ESTIMATION PRINT=10 MAXEVAL=9999 NOABORT METHOD=1 
            POSTHOC INTERACTION NOOBT 
$COVARIANCE PRINT=E UNCONDITIONAL
$TABLE 
 ID RND WT AED CL VC KA
 ONEHEADER NOPRINT NOAPPEND FIRSTONLY FILE=run025.csv
  
NM-TRAN MESSAGES 
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.
  
License Registered to: Occams Cooperatie UA
Expiration Date:    14 JUN 2025
Current Date:       10 APR 2024
Days until program expires : 429
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 BRV Adult popPK model ; run025 analysis simulated data
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:    10712
 NO. OF DATA ITEMS IN DATA SET:  12
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   9
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
  12   4   6   0   7   8   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID RND OCC TIME MDV AMT SS II DV WT AED EVID
0(NONBLANK) LABELS FOR PRED-DEFINED ITEMS:
 CL VC KA
0FORMAT FOR DATA:
 (E5.0,2E4.0,E8.0,E2.0,E9.0,E2.0,E3.0,E8.0,E5.0,E2.0,1F2.0)

 TOT. NO. OF OBS RECS:     5820
 TOT. NO. OF INDIVIDUALS:     1377
0LENGTH OF THETA:   7
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   3
0DEFAULT OMEGA BOUNDARY TEST OMITTED:   YES
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.1000E-02     0.2200E+00     0.1000E+07
  0.1000E-02     0.3700E+01     0.1000E+07
  0.1000E-01     0.5000E+02     0.1000E+07
  0.1000E-01     0.1500E+01     0.1000E+07
  0.7500E+00     0.7500E+00     0.7500E+00
  0.1000E+01     0.1000E+01     0.1000E+01
 -0.1000E+07     0.2000E+00     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.1000E+00
 0.0000E+00   0.3000E+00
 0.0000E+00   0.0000E+00   0.1300E+01
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:             YES
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
0RECORDS ONLY:    FIRSTONLY
04 COLUMNS APPENDED:    NO
 PRINTED:                NO
 HEADERS:               ONE
 FILE TO BE FORWARDED:   NO
 FORMAT:                S1PE11.4
 IDFORMAT:
 LFORMAT:
 RFORMAT:
 FIXED_EFFECT_ETAS:
0USER-CHOSEN ITEMS:
 ID RND WT AED CL VC KA
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 ONE COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN2)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   3
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   ELIMINATION RATE (K) IS BASIC PK PARAMETER NO.:  1
   ABSORPTION RATE (KA) IS BASIC PK PARAMETER NO.:  3

0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            4           *           *           *           *
    3            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:     12
   TIME DATA ITEM IS DATA ITEM NO.:          4
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   6
   STEADY STATE DATA ITEM IS DATA ITEM NO.:  7
   INTERVAL DATA ITEM IS DATA ITEM NO.:      8

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
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
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): run025.ext
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

 
 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI
 
 MONITORING OF SEARCH:

 
0ITERATION NO.:    0    OBJECTIVE VALUE:  -13679.6482216533        NO. OF FUNC. EVALS.:   7
 CUMULATIVE NO. OF FUNC. EVALS.:        7
 NPARAMETR:  2.2000E-01  3.7000E+00  5.0000E+01  1.5000E+00  2.0000E-01  1.0000E-01  3.0000E-01  1.3000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   1.0714E+03 -7.3621E+02 -4.1032E+02 -9.6185E+00 -1.9285E+03  5.2119E+02  3.6210E+02  1.7527E+02
 
0ITERATION NO.:   10    OBJECTIVE VALUE:  -13939.8917862056        NO. OF FUNC. EVALS.:  87
 CUMULATIVE NO. OF FUNC. EVALS.:       94
 NPARAMETR:  2.0718E-01  3.6183E+00  5.1284E+01  1.4027E+00  3.1857E-01  6.9422E-02  1.1313E-01  8.2288E-01
 PARAMETER:  3.9695E-02  7.7673E-02  1.2536E-01  3.2463E-02  1.5928E-01 -8.2481E-02 -3.8763E-01 -1.2865E-01
 GRADIENT:  -7.3536E+01  5.0754E+01 -2.2816E+00  2.9292E+00  4.0398E+00  1.3384E+01  1.3571E+01  9.7388E+00
 
0ITERATION NO.:   17    OBJECTIVE VALUE:  -13940.5864188288        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      189
 NPARAMETR:  2.0836E-01  3.6053E+00  5.1148E+01  1.3909E+00  3.2463E-01  6.9036E-02  1.0563E-01  7.9337E-01
 PARAMETER:  4.5405E-02  7.4073E-02  1.2270E-01  2.3984E-02  1.6232E-01 -8.5268E-02 -4.2193E-01 -1.4691E-01
 GRADIENT:   6.1571E-03 -1.3778E-02 -6.6406E-03  6.2196E-03 -5.3090E-03  2.7865E-03 -3.7489E-03  6.1538E-03
 
 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      189
 NO. OF SIG. DIGITS IN FINAL EST.:  3.8

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.
 
 ETABAR:         1.7686E-03 -2.8121E-02 -4.6123E-02
 SE:             6.5225E-03  4.2441E-03  1.1551E-02
 N:                    1248        1248        1248
 
 P VAL.:         7.8627E-01  3.4686E-11  6.5288E-05
 
 ETASHRINKSD(%)  1.2268E+01  5.3849E+01  5.4168E+01
 ETASHRINKVR(%)  2.3031E+01  7.8701E+01  7.8994E+01
 EBVSHRINKSD(%)  1.2442E+01  5.5206E+01  5.3046E+01
 EBVSHRINKVR(%)  2.3336E+01  7.9935E+01  7.7953E+01
 RELATIVEINF(%)  7.5700E+01  1.7116E+01  1.9005E+01
 EPSSHRINKSD(%)  1.2421E+01
 EPSSHRINKVR(%)  2.3299E+01
 
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):         5820
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    10696.4445265024     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -13940.5864188288     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -3244.14189232643     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                          3744
  
 #TERE:
 Elapsed estimation  time in seconds:    10.19
 Elapsed covariance  time in seconds:     8.12
 Elapsed postprocess time in seconds:     0.11
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************   -13940.586       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         2.08E-01  3.61E+00  5.11E+01  1.39E+00  7.50E-01  1.00E+00  3.25E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        6.90E-02
 
 ETA2
+        0.00E+00  1.06E-01
 
 ETA3
+        0.00E+00  0.00E+00  7.93E-01
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        2.63E-01
 
 ETA2
+        0.00E+00  3.25E-01
 
 ETA3
+        0.00E+00  0.00E+00  8.91E-01
 


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


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         2.40E-03  4.21E-02  1.42E+00  7.61E-02 ......... .........  1.68E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        3.75E-03
 
 ETA2
+       .........  1.95E-02
 
 ETA3
+       ......... .........  1.34E-01
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        7.13E-03
 
 ETA2
+       .........  3.00E-02
 
 ETA3
+       ......... .........  7.53E-02
 


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
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        5.75E-06
 
 TH 2
+        5.16E-06  1.77E-03
 
 TH 3
+        5.48E-04 -9.80E-03  2.01E+00
 
 TH 4
+        1.09E-05 -3.66E-04  4.24E-02  5.80E-03
 
 TH 5
+       ......... ......... ......... ......... .........
 
 TH 6
+       ......... ......... ......... ......... ......... .........
 
 TH 7
+       -1.26E-06 -4.83E-04  2.17E-03  9.30E-05 ......... .........  2.83E-04
 
 OM11
+       -2.00E-08  1.16E-05 -1.81E-04 -1.86E-05 ......... ......... -5.41E-06  1.41E-05
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       -9.62E-07 -8.14E-05  1.54E-02  3.97E-04 ......... .........  5.85E-06 -1.29E-05 ......... .........  3.79E-04
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       -5.23E-05 -7.17E-05 -5.71E-02 -1.56E-04 ......... ......... -5.05E-05 -3.06E-06 ......... ......... -1.01E-03 .........
          1.80E-02
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        2.40E-03
 
 TH 2
+        5.11E-02  4.21E-02
 
 TH 3
+        1.61E-01 -1.64E-01  1.42E+00
 
 TH 4
+        5.99E-02 -1.14E-01  3.93E-01  7.61E-02
 
 TH 5
+       ......... ......... ......... ......... .........
 
 TH 6
+       ......... ......... ......... ......... ......... .........
 
 TH 7
+       -3.13E-02 -6.81E-01  9.09E-02  7.26E-02 ......... .........  1.68E-02
 
 OM11
+       -2.22E-03  7.32E-02 -3.41E-02 -6.50E-02 ......... ......... -8.57E-02  3.75E-03
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       -2.06E-02 -9.93E-02  5.58E-01  2.68E-01 ......... .........  1.78E-02 -1.76E-01 ......... .........  1.95E-02
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       -1.63E-01 -1.27E-02 -3.00E-01 -1.53E-02 ......... ......... -2.23E-02 -6.09E-03 ......... ......... -3.87E-01 .........
          1.34E-01
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        1.88E+05
 
 TH 2
+       -6.23E+02  1.09E+03
 
 TH 3
+       -7.38E+01  2.93E+00  8.61E-01
 
 TH 4
+       -1.96E+02  1.03E+01 -4.06E+00  2.10E+02
 
 TH 5
+       ......... ......... ......... ......... .........
 
 TH 6
+       ......... ......... ......... ......... ......... .........
 
 TH 7
+        4.97E+02  1.83E+03 -2.90E-01 -2.00E+01 ......... .........  6.69E+03
 
 OM11
+        4.77E+03 -7.45E+00 -2.30E+01  1.18E+02 ......... .........  1.53E+03  7.51E+04
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+        5.32E+03  1.41E+02 -2.87E+01 -9.43E+01 ......... .........  5.26E+02  3.81E+03 ......... .........  4.55E+03
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+        6.10E+02  2.50E+01  8.78E-01 -1.69E+01 ......... .........  5.62E+01  1.73E+02 ......... .........  1.82E+02 .........
          7.05E+01
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                      EIGENVALUES OF COR MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

             1         2         3         4         5         6         7         8
 
         3.05E-01  3.69E-01  5.50E-01  9.55E-01  9.77E-01  1.06E+00  1.63E+00  2.16E+00
 
 Elapsed finaloutput time in seconds:     0.06
 #CPUT: Total CPU Time in Seconds,        8.109
Stop Time: 
10/04/2024 
09:54
