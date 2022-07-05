*> \brief \b SGETRS_EC
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SGETRS_EC + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgetrs_ec.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgetrs_ec.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgetrs_ec.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGETRS_EC( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, 
*      $                     INFO, FLAG_REPORT, INFO_ARRAY, CONTEXT )
*
*       .. Scalar Arguments ..
*       CHARACTER          TRANS
*       INTEGER            INFO, LDA, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       INTEGER            FLAG_REPORT( 2 ), INFO_ARRAY( 7 )
*       REAL               A( LDA, * ), B( LDB, * )
*       ..
*       .. Pointer Arguments ..
*       POINTER            CONTEXT ... advice requested
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGETRS_EC solves a system of linear equations
*>    A * X = B  or  A**T * X = B
*> with a general N-by-N matrix A using the LU factorization computed
*> by SGETRF_EC.
*>
*> SGETRS_EC also provides new exception handling and 
*> reporting capabilities.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          Specifies the form of the system of equations:
*>          = 'N':  A * X = B  (No transpose)
*>          = 'T':  A**T* X = B  (Transpose)
*>          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of columns
*>          of the matrix B.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          The factors L and U from the factorization A = P*L*U
*>          as computed by SGETRF.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          The pivot indices from SGETRF; for 1<=i<=N, row i of the
*>          matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is REAL array, dimension (LDB,NRHS)
*>          On entry, the right hand side matrix B.
*>          On exit, the solution matrix X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          INFO is defined below depending on FLAG_REPORT
*> \endverbatim
*>
*> \param[in] FLAG_REPORT
*> \verbatim
*>          FLAG_REPORT is INTEGER array, dimension(2)
*>          FLAG_REPORT(1) defines what kinds of exceptions to report 
*>          using INFO and possibly also INFO_ARRAY for more details.
*>          FLAG_REPORT(1)
*>           <= -1 turns off all error checking, so INFO=0 is 
*>                 returned.
*>            =  0 does standard argument checking:
*>                 INFO = 0  means successful exit
*>                 INFO = -i means the i-th (non-floating point) 
*>                           argument had an illegal value 
*>                           (first error found is reported)
*>                 Using INFO to report the above errors has priority 
*>                 over reporting any of the errors described below. 
*>                 More generally, an error that would be found with 
*>                 a lower value of FLAG_REPORT(1) has priority to 
*>                 report using INFO than an error that would only
*>                 be found with a higher value of FLAG_REPORT(1).  
*>            =  1 also checks for Infs and NaNs in inputs and 
*>                 outputs, if INFO is not already nonzero:
*>                 INFO = -4 means A contained an Inf or NaN on 
*>                           input; execution continues.
*>                 INFO = -7 means B contained an Inf or NaN on 
*>                           input but A did not; 
*>                           execution continues.
*>                 INFO = 1  means B contained an Inf or NaN on 
*>                           output but neither A nor B did on input. 
*>                 Since A is an input variable, it is not checked
*>                 on output.
*>                 Each input and output will also be checked if 
*>                 FLAG_REPORT(2) = 1, 2, or 3 and reported in 
*>                 INFO_ARRAY as described below.
*>           >=  2 has the same behavior as 1, since there are no 
*>                 internal calls to LAPACK routines with INFO 
*>                 parameters to be checked.
*>
*>          FLAG_REPORT(2) defines how to report the exceptions 
*>          requested by FLAG_REPORT(1).
*>            If FLAG_REPORT(1) <= -1, FLAG_REPORT(2) is ignored and 
*>            INFO=0 is returned. 
*>            Otherwise, FLAG_REPORT(2)
*>             <=  0 only returns the value of INFO described above.
*>              =  1 also returns INFO_ARRAY, as described below.
*>              =  2 means that SGETRS_EC will also call 
*>                   REPORT_EXCEPTIONS to report INFO_ARRAY, if INFO 
*>                   is nonzero.
*>              =  3 has the same behavior as 2. (If SGETRS called
*>                   any LAPACK routines with INFO parameters 
*>                   internally then they would call
*>                   REPORT_EXCEPTIONS too, but there are none.)
*>             >=  4 means that SGETRS_EC will call 
*>                   GET_FLAGS_TO_REPORT(CONTEXT,FLAG_REPORT) 
*>                   to get values of FLAG_REPORT to use, overriding 
*>                   input values. The user needs to have called 
*>                   SET_FLAGS_TO_REPORT(CONTEXT,FLAG_REPORT) 
*>                   before calling SGETRS_EC in order to set 
*>                   FLAG_REPORT, otherwise the default is 
*>                   FLAG_REPORT = [0, 0]. The input array
*>                   FLAG_REPORT will not be overwritten.
*> \endverbatim
*>
*> \param[in,out] INFO_ARRAY
*> \verbatim
*>          INFO_ARRAY is INTEGER array, dimension( 8 )
*>          If FLAG_REPORT(1) <= -1 or FLAG_REPORT(2) <= 0, 
*>          INFO_ARRAY is not accessed. Otherwise:
*>          INFO_ARRAY(1)
*>              = value of INFO from standard argument checking 
*>                (as defined by FLAG_REPORT(1) = 0)
*>          INFO_ARRAY(2)
*>              = value of FLAG_REPORT(1) used to determine the rest 
*>                of INFO_ARRAY
*>          INFO_ARRAY(3)
*>              = value of FLAG_REPORT(2) used to determine the rest 
*>                of INFO_ARRAY
*>          INFO_ARRAY(4)
*>              = value of INFO as specified by FLAG_REPORT(1) above
*>          INFO_ARRAY(5)
*>              = 2 = number of input/output arguments reported on
*>          INFO_ARRAY(6)
*>              = 0 = number of internal LAPACK calls reported on
*>          INFO_ARRAY(7) reports exceptions in A, as specified by 
*>            FLAG_REPORT
*>              = -1 if not checked (default)
*>              =  0 if checked and contains no Infs or NaNs
*>              =  1 if checked and contains an Inf or NaN on input
*>            If INFO_ARRAY(7) = 0 or 1 on input, then A will not
*>            be checked again on input. Input values < -1 or > 1 
*>            will be treated the same as -1, i.e. not checked.
*>          INFO_ARRAY(8) reports exceptions in B, as specified by
*>            FLAG_REPORT
*>              = -1 if not checked (default)
*>              =  0 if checked and contains no Infs or NaNs
*>              =  1 if checked and contains an Inf or NaN on input
*>                   but not output
*>              =  2 if checked and contains an Inf or NaN on output 
*>                   but not input
*>              =  3 if checked and contains an Inf or NaN on input 
*>                   and output
*>            As above, if INFO_ARRAY(8) = 0 or 1 on input, then B
*>            will not be checked again on input. Input values < -1
*>            or > 1 will be treated the same as -1, 
*>            i.e. not checked.
*> \endverbatim
*>
*> \param[in] CONTEXT
*> \verbatim
*>           CONTEXT is POINTER to an "opaque object"
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup realGEcomputational
*
*  =====================================================================
      SUBROUTINE SGETRS_EC( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, 
     $                      INFO, FLAG_REPORT, INFO_ARRAY, CONTEXT )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), B( LDB, * )
      INTEGER            FLAG_REPORT( 2 ), INFO_ARRAY( * )
*     ..
*     .. Pointer Arguments
      POINTER            CONTEXT
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      CHARACTER, DIMENSION(6), PARAMETER :: 
     $   ROUTINE_NAME = (/ 'S','G','E','T','R','S' /)
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN, CALL_REPORT_EXCEPTIONS
      INTEGER            WHAT, HOW
*     ..
*     .. Local Arrays ..
      INTEGER            FLAG_REPORT_INTERNAL(2), FLAG_REPORT_CALL(2)    
      INTEGER            INFO_INTERNAL( 2 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLASWP, STRSM, XERBLA
      EXTERNAL           CHECKINIT1, CHECKINIT2
      EXTERNAL           SGECHECKARG
      EXTERNAL           UPDATE_INFO
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*     Initialize error checking flags
      CALL CHECKINIT1(FLAG_REPORT, FLAG_REPORT_INTERNAL,  
     $                FLAG_REPORT_CALL, CALL_REPORT_EXCEPTIONS, 
     $                CONTEXT ) 
      WHAT = FLAG_REPORT_INTERNAL( 1 )
      HOW = FLAG_REPORT_INTERNAL( 2 )     
      IF (WHAT .EQ. -1 ) GOTO 100
*
*     Check for standard input errors
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
*
*     Initialize error flags
      CALL CHECKINIT2( FLAG_REPORT_INTERNAL, INFO, INFO_INTERNAL, 
     $                 INFO_ARRAY, 2, 0)  
*
      IF( INFO.NE.0 ) THEN
         IF (CALL_REPORT_EXCEPTIONS) 
     $      CALL REPORT_EXCEPTIONS(CONTEXT,6,ROUTINE_NAME,INFO_ARRAY)
         CALL XERBLA( 'SGETRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
*     Check for exceptional inputs in A 
      CALL SGECHECKARG(FLAG_REPORT_INTERNAL, N, N, A, LDA, 
     $                 INFO, INFO_INTERNAL, INFO_ARRAY, 4, 0, 0, 7)
*     Check for exceptional inputs in B 
      CALL SGECHECKARG(FLAG_REPORT_INTERNAL, N, NRHS, B, LDB, 
     $                 INFO, INFO_INTERNAL, INFO_ARRAY, 7, 2, 0, 8)
*
100   CONTINUE
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( NOTRAN ) THEN
*
*        Solve A * X = B.
*
*        Apply row interchanges to the right hand sides.
*
         CALL SLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
*
*        Solve L*X = B, overwriting B with X.
*
         CALL STRSM('Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve U*X = B, overwriting B with X.
*
         CALL STRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
*
*        Solve A**T * X = B.
*
*        Solve U**T *X = B, overwriting B with X.
*
         CALL STRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
*
*        Solve L**T *X = B, overwriting B with X.
*
         CALL STRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Apply row interchanges to the solution vectors.
*
         CALL SLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
*
*     Check for errors before returning
*
      IF (WHAT .EQ. -1) RETURN
*     Check for exceptional outputs in B 
      CALL SGECHECKARG(FLAG_REPORT_INTERNAL, N, NRHS, B, LDB, 
     $                 INFO, INFO_INTERNAL, INFO_ARRAY, 7, 3, 1, 8)  
*
*     Update INFO and INFO_ARRAY
      CALL UPDATE_INFO( INFO, INFO_ARRAY, INFO_INTERNAL )
      IF (CALL_REPORT_EXCEPTIONS .AND. INFO .NE. 0) 
     $   CALL REPORT_EXCEPTIONS(CONTEXT, 6, ROUTINE_NAME, INFO_ARRAY)
      RETURN
*
*     End of SGETRS_EC
*
      END
