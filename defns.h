#define SUCCESS 0
#define ERROR   1
#define NO_MEM  2
#define SEGY_NOFILE 3
#define FILENOTFOUND	4
#define SHORT_FILE 57
#define WRONG_DOMAIN (-2)
#define VALUE_NOT_SET_FLOAT -9999999.
/*
*** Data types
*/
#define SEISMIC_DATA   1
#define FFT_DATA       2
#define LFFT_DATA      3
#define CEPT_DATA      4
#define RMS_DATA       5
#define GABOR_DATA     6
#define SWFT_DATA      7
#define WMA_DATA       8
#define WAVELET_DATA   9
#define WAVEFRQ_DATA  10
#define CEPSTRUM_DATA 11
#define WGRAPH	      31
#define XYGRAPH	      32
/*
*** Data domain
*/
#define REAL	      1
#define COMPLEX_POLAR 2
#define COMPLEX_CART  3
/*
*** Data units
*/
#define DEAD         0
#define TIME         1
#define FOURIER      2
#define QUEFREQUENCY 3
/*
*** Wavelet Precision
*/
#define PRECISION_MIN   1
#define PRECISION_AVG   2
#define PRECISION_MAX   3
#define PRECISION_8192  4
#define PRECISION_16384 5
/*
*** Decay types and Taper options
*/
#define DECAY_SPHERICAL   1
#define DECAY_EXPONENTIAL 2
#define ORIGIN_MIDDLE     1
#define ORIGIN_RIGHT      2
#define ORIGIN_LEFT       3
#define COSINE_QUARTER    1
#define COSINE_HALF       2
/*
*** Local Attributes
*/
#define LOCAL_AMPLITUDE     1
#define LOCAL_FREQUENCY     2
#define LOCAL_DOM_FREQUENCY 3
#define LOCAL_CORRELATION   4
#define LOCAL_PERIOD        5
#define LOCAL_PHASE         6
#define LOCAL_PHASE_SHIFT   7
/*
*** Picking parameters
*/
#define SNAP_NONE              1
#define SNAP_PICK              2
#define SNAP_TROUGH            3
#define SNAP_PLUS_MINUS        4
#define SNAP_MINUS_PLUS        5
#define SNAP_ZERO              6
#define SNAP_MINIMUM           7
#define SNAP_MAXIMUM           8
#define SNAP_MINIMUM_ABSOLUTE  9
#define SNAP_MAXIMUM_ABSOLUTE 10
/*
*** Eaa modes
*/
#define Q_NONE	    0
#define INTEGRAL   44
#define DERIVATIVE 45
#define RELATIVE   46
#define ABSOLUTE   47
#define Q_LINEAR   48
/*
*** Eaa f1-f2 options
*/
#define F1FLOAT    1
#define F1FIX      2
#define F1F2FLOAT  3
#define F2FLOAT    4
#define F2MAXFLOAT 5
#define F2ABS_LOW  6
#define F2ABS_HI   7
#define F2WINDOW   8
/*
*** LHA types
*/
#define LHA_NONE     0
#define PLASTIC      1
#define PLAS_ELASTIC 2
#define SPEC_SMOOTH  3
/*
*** Picking directions
*/
#define FIRST        1
#define NEXT         2
#define LAST         3
#define FIRST_LAST   4
/*
*** Display Clip Types
*/
#define PERC  1
#define VCLIP 2
#define BCLIP 3
/*
*** Trend analysis types
*/
#define TREND_RMS 1
#define TREND_ARITH 2
#define TREND_GEOM 3
#define TREND_HARM 4
#define TREND_MEDIAN 5
#define TREND_SPIKE 6
/*
*** DFQ output options
*/
#define DOMFRQ   0
#define EAAVAL   1
#define A1VAL    2
#define AMPMAX   3
#define WAVMAX   4
#define WIDTH    5
#define ENERGY   6
#define RCVAL    7
#define SPECTRA  8
#define FREQ2    9
#define DELTAF  10
#define MAXOUT  11
/*
*** General
*/
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif
#define YES 1
#define NO  0
#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif
#ifndef MAX
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#define ISIZE sizeof(int)
#define FSIZE sizeof(float)
#define DSIZE sizeof(double)

#define	STREQ(s,t) (strcmp(s,t) == 0)
#define	STRLT(s,t) (strcmp(s,t) < 0)
#define	STRGT(s,t) (strcmp(s,t) > 0)

#define OK_STAT 0
#define OK     1
#define CANCEL 2
#define VALUE_NOT_SET_FLOAT  -9999999.
/*
*** Work flags
*/
#define STOP 999
#define RUN  1
#define RUNNING 2
#define IDLE 5
