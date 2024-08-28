10/04/2024 
09:57
$PROB  BRV Paediatric popPK model ; run625 analysis simulated data

$DATA  DAT_P.csv IGNORE=@ 
$INPUT ID TIME MDV AMT RATE EVID SS II DV WT AGE AED 

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
 ID WT AGE AED CL VC KA
 ONEHEADER NOPRINT NOAPPEND FIRSTONLY FILE=run625.csv
  
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
 BRV Paediatric popPK model ; run625 analysis simulated data
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:     1971
 NO. OF DATA ITEMS IN DATA SET:  12
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   9
 MDV DATA ITEM IS DATA ITEM NO.:  3
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   5   7   8   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME MDV AMT RATE EVID SS II DV WT AGE AED
0(NONBLANK) LABELS FOR PRED-DEFINED ITEMS:
 CL VC KA
0FORMAT FOR DATA:
 (E4.0,E7.0,E2.0,E8.0,3E2.0,E3.0,E8.0,2E6.0,E2.0)

 TOT. NO. OF OBS RECS:      855
 TOT. NO. OF INDIVIDUALS:      235
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
 ID WT AGE AED CL VC KA
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
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4
   DOSE RATE DATA ITEM IS DATA ITEM NO.:     5
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
 RAW OUTPUT FILE (FILE): run625.ext
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

 
0ITERATION NO.:    0    OBJECTIVE VALUE:  -430.815172229268        NO. OF FUNC. EVALS.:   7
 CUMULATIVE NO. OF FUNC. EVALS.:        7
 NPARAMETR:  2.2000E-01  3.7000E+00  5.0000E+01  1.5000E+00  2.0000E-01  1.0000E-01  3.0000E-01  1.3000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:  -8.2420E+02 -3.7458E+02 -2.1001E+02  4.6565E+00 -3.0675E+02  5.2884E+01  9.3778E+01  4.7388E+01
 
0ITERATION NO.:   10    OBJECTIVE VALUE:  -724.350329706708        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:       90
 NPARAMETR:  3.1152E-01  4.2078E+00  6.5032E+01  2.1666E+00  2.3482E-01  5.5410E-02  1.8103E-02  1.1728E-01
 PARAMETER:  4.4917E-01  2.2865E-01  3.6291E-01  4.6976E-01  1.1741E-01 -1.9520E-01 -1.3038E+00 -1.1028E+00
 GRADIENT:   1.2339E+01  1.8700E+01 -1.1165E+01  5.4447E+00  1.2710E+01  7.6194E+00  6.0724E-01  3.7533E+00
 
0ITERATION NO.:   20    OBJECTIVE VALUE:  -725.996453186689        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      170
 NPARAMETR:  3.1054E-01  4.1829E+00  6.4850E+01  2.0918E+00  2.3585E-01  5.3076E-02  1.8723E-02  2.1819E-04
 PARAMETER:  4.4603E-01  2.2271E-01  3.6009E-01  4.3448E-01  1.1792E-01 -2.1672E-01 -1.2870E+00 -4.2462E+00
 GRADIENT:  -2.6457E-01  3.8835E-02 -1.3416E-01 -5.6574E-02  3.9423E-01  3.4344E-02 -3.4799E-02  2.8283E-03
 
0ITERATION NO.:   30    OBJECTIVE VALUE:  -726.000104013998        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      309
 NPARAMETR:  3.1073E-01  4.1864E+00  6.4916E+01  2.0966E+00  2.3571E-01  5.3085E-02  1.8845E-02  8.2904E-08
 PARAMETER:  4.4664E-01  2.2354E-01  3.6112E-01  4.3677E-01  1.1786E-01 -2.1664E-01 -1.2838E+00 -8.1840E+00
 GRADIENT:  -4.3289E-03  1.4322E-02 -4.8941E-03 -5.0519E-04  1.0110E-02 -1.0968E-03 -5.9629E-04  1.0839E-06
 
0ITERATION NO.:   34    OBJECTIVE VALUE:  -726.000104626168        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      382
 NPARAMETR:  3.1073E-01  4.1864E+00  6.4917E+01  2.0966E+00  2.3571E-01  5.3085E-02  1.8847E-02  2.9293E-09
 PARAMETER:  4.4664E-01  2.2353E-01  3.6113E-01  4.3677E-01  1.1786E-01 -2.1664E-01 -1.2837E+00 -9.8554E+00
 GRADIENT:   7.2318E-05 -1.2994E-03  1.3995E-04 -7.1730E-05 -2.8882E-04  1.7500E-04  1.2095E-04  1.1626E-06
 
 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      382
 NO. OF SIG. DIGITS IN FINAL EST.:  3.9

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.
 
 ETABAR:         6.5740E-04 -1.8695E-03  2.3404E-11
 SE:             1.1168E-02  2.8461E-03  1.0295E-10
 N:                     232         232         232
 
 P VAL.:         9.5306E-01  5.1127E-01  8.2017E-01
 
 ETASHRINKSD(%)  2.6010E+01  6.8354E+01  9.9997E+01
 ETASHRINKVR(%)  4.5255E+01  8.9986E+01  1.0000E+02
 EBVSHRINKSD(%)  2.6406E+01  6.8863E+01  9.9997E+01
 EBVSHRINKVR(%)  4.5839E+01  9.0305E+01  1.0000E+02
 RELATIVEINF(%)  5.3913E+01  7.6237E+00  7.1776E-08
 EPSSHRINKSD(%)  7.4167E+00
 EPSSHRINKVR(%)  1.4283E+01
 
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          855
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1571.38489177999     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -726.000104626168     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       845.384787153822     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           696
  
 #TERE:
 Elapsed estimation  time in seconds:     3.23
 Elapsed covariance  time in seconds:     1.19
 Elapsed postprocess time in seconds:     0.02
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -726.000       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         3.11E-01  4.19E+00  6.49E+01  2.10E+00  7.50E-01  1.00E+00  2.36E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        5.31E-02
 
 ETA2
+        0.00E+00  1.88E-02
 
 ETA3
+        0.00E+00  0.00E+00  2.93E-09
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        2.30E-01
 
 ETA2
+        0.00E+00  1.37E-01
 
 ETA3
+        0.00E+00  0.00E+00  5.41E-05
 


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
 
         9.14E-03  1.10E-01  2.53E+00  2.92E-01 ......... .........  3.87E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        8.11E-03
 
 ETA2
+       .........  1.21E-02
 
 ETA3
+       ......... .........  4.40E-09
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        1.76E-02
 
 ETA2
+       .........  4.39E-02
 
 ETA3
+       ......... .........  4.07E-05
 


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
+        8.35E-05
 
 TH 2
+        1.63E-04  1.21E-02
 
 TH 3
+        2.15E-03  1.57E-02  6.38E+00
 
 TH 4
+        4.84E-05 -2.12E-03  2.48E-01  8.55E-02
 
 TH 5
+       ......... ......... ......... ......... .........
 
 TH 6
+       ......... ......... ......... ......... ......... .........
 
 TH 7
+       -6.99E-06 -2.58E-03  8.77E-04  8.66E-04 ......... .........  1.50E-03
 
 OM11
+       -1.21E-05 -9.52E-05 -2.68E-04 -3.08E-04 ......... .........  2.67E-05  6.58E-05
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       -2.68E-05 -6.50E-05  1.53E-02 -2.56E-04 ......... .........  1.64E-05 -6.84E-06 ......... .........  1.45E-04
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+        2.90E-12 -4.09E-11  1.07E-09  1.56E-10 ......... .........  7.94E-13  1.09E-12 ......... ......... -5.12E-12 .........
          1.94E-17
 
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
+        9.14E-03
 
 TH 2
+        1.62E-01  1.10E-01
 
 TH 3
+        9.31E-02  5.64E-02  2.53E+00
 
 TH 4
+        1.81E-02 -6.58E-02  3.36E-01  2.92E-01
 
 TH 5
+       ......... ......... ......... ......... .........
 
 TH 6
+       ......... ......... ......... ......... ......... .........
 
 TH 7
+       -1.97E-02 -6.05E-01  8.96E-03  7.65E-02 ......... .........  3.87E-02
 
 OM11
+       -1.64E-01 -1.07E-01 -1.31E-02 -1.30E-01 ......... .........  8.49E-02  8.11E-03
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       -2.43E-01 -4.90E-02  5.03E-01 -7.28E-02 ......... .........  3.52E-02 -7.00E-02 ......... .........  1.21E-02
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+        7.21E-02 -8.46E-02  9.60E-02  1.21E-01 ......... .........  4.66E-03  3.06E-02 ......... ......... -9.65E-02 .........
          4.40E-09
 
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
+        1.49E+04
 
 TH 2
+       -2.04E+02  1.39E+02
 
 TH 3
+       -1.89E+01 -7.13E-01  2.96E-01
 
 TH 4
+        7.50E+01  3.22E+00 -9.89E-01  1.56E+01
 
 TH 5
+       ......... ......... ......... ......... .........
 
 TH 6
+       ......... ......... ......... ......... ......... .........
 
 TH 7
+       -4.28E+02  2.34E+02 -2.87E-01 -6.06E+00 ......... .........  1.07E+03
 
 OM11
+        3.43E+03  8.53E+01 -1.15E+01  1.06E+02 ......... ......... -2.22E+02  1.68E+04
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+        4.96E+03  9.55E+01 -3.79E+01  1.51E+02 ......... ......... -6.64E+01  2.86E+03 ......... .........  1.24E+04
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       -1.09E+09  3.49E+08 -1.64E+07 -4.12E+07 ......... .........  5.74E+08 -7.36E+08 ......... .........  3.44E+09 .........
          5.47E+16
 
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
 
         2.68E-01  3.69E-01  8.29E-01  8.83E-01  1.03E+00  1.33E+00  1.58E+00  1.70E+00
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,        1.906
Stop Time: 
10/04/2024 
09:57
