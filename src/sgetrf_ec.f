*> \brief \b SGETRF_EC
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SGETRF_EC + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgetrf_ec.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgetrf_ec.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgetrf_ec.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGETRF_EC( M, N, A, LDA, IPIV, 
*      $                   INFO, FLAG_REPORT, INFO_ARRAY, CONTEXT )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       INTEGER            FLAG_REPORT( 2 ), INFO_ARRAY( * )
*       REAL               A( LDA, * )
*       ..
*       .. Pointer Arguments ..
*       POINTER            CONTEXT ... advice requested
*
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGETRF_EC computes an LU factorization of a general M-by-N matrix A
*> using partial pivoting with row interchanges.
*>
*> The factorization has the form
*>    A = P * L * U
*> where P is a permutation matrix, L is lower triangular with unit
*> diagonal elements (lower trapezoidal if m > n), and U is upper
*> triangular (upper trapezoidal if m < n).
*>
*> This is the right-looking Level 3 BLAS version of the algorithm.
*>
*> SGETRF_EC also provides new exception handling and 
*> reporting capabilities.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          On entry, the M-by-N matrix to be factored.
*>          On exit, the factors L and U from the factorization
*>          A = P*L*U; the unit diagonal elements of L are not stored.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (min(M,N))
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the
*>          matrix was interchanged with row IPIV(i).
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
*>          FLAG_REPORT(1) defines what kinds of exceptions to report,
*>          using INFO and possibly also INFO_ARRAY for more details.
*>          FLAG_REPORT(1) 
*>           <= -1 turns off all error checking, so INFO=0 is 
*>                 returned.
*>            =  0 does standard argument checking:
*>                 INFO = 0  means successful exit
*>                 INFO = -i means the i-th (non-floating point) 
*>                           argument had an illegal value 
*>                           (first error found is reported)
*>                 INFO = i means U(i,i) = 0. The factorization
*>                          has been completed, but the factor U is 
*>                          exactly singular, and division by zero 
*>                          will occur if it is used to solve a system 
*>                          of equations.
*>                 Using INFO to report the above errors has priority 
*>                 over reporting any of the errors described below. 
*>                 More generally, an error that would be found with a 
*>                 lower value of FLAG_REPORT(1) has priority to 
*>                 report using INFO than an error that would only
*>                 be found with a higher value of FLAG_REPORT(1).  
*>            =  1 also checks for Infs and NaNs in inputs and  
*>                 outputs, if INFO not already nonzero:
*>                 INFO = -3 means A contained an Inf or NaN on 
*>                           input; execution continues.
*>                 INFO = N+1 means that A contained an Inf or NaN on 
*>                        output but not input.
*>                 A will also be checked on input and output if 
*>                 FLAG_REPORT(2) = 1, 2 or 3 and reported in 
*>                 INFO_ARRAY as described below.
*>            >= 2 also checks for Infs and NaNs as inputs or outputs 
*>                 in all internal LAPACK calls in the call chain, if 
*>                 INFO is not already nonzero, and 
*>                 FLAG_REPORT(2) = 1, 2 or 3. In this case: 
*>                 INFO = N+2 means that either the first call
*>                        to SGETRF2_EC had an Inf or NaN as an 
*>                        input or output as above, or a subroutine
*>                        in its call tree did.
*>                 INFO = N+3 means that either the second call
*>                        to SGETRF2_EC had an Inf or NaN as an 
*>                        input or output as above, or a subroutine
*>                        in its call tree did.
*>                 Each input, output and internal LAPACK call will
*>                 also be checked if FLAG_REPORT(2) = 1, 2 or 3 and 
*>                 reported in INFO_ARRAY as described below.
*>
*>          FLAG_REPORT(2) defines how to report the exceptions 
*>          requested by FLAG_REPORT(1).
*>            If FLAG_REPORT(1) <= -1, FLAG_REPORT(2) is ignored and 
*>            INFO=0 is returned. 
*>            Otherwise, FLAG_REPORT(2)
*>             <=  0 only returns the value of INFO described above.
*>              =  1 also returns INFO_ARRAY, as described below.
*>              =  2 means that SGETRF_EC will also call 
*>                   REPORT_EXCEPTIONS to report INFO_ARRAY, if INFO 
*>                   is nonzero.
*>              =  3 means that all calls in the call tree will also
*>                   call REPORT_EXCEPTIONS, if the value of INFO they
*>                   return is nonzero.
*>             >=  4 means that SGETRF_EC will call 
*>                   GET_FLAGS_TO_REPORT(CONTEXT,FLAG_REPORT) 
*>                   to get values of FLAG_REPORT to use, overriding 
*>                   input values. The user needs to have called 
*>                   SET_FLAGS_TO_REPORT(CONTEXT,FLAG_REPORT) 
*>                   before calling SGETRF_EC in order to set 
*>                   FLAG_REPORT, otherwise the default is 
*>                   FLAG_REPORT = [0, 0]. The input array FLAG_REPORT 
*>                   will not be overwritten.
*> \endverbatim
*>
*> \param[out] INFO_ARRAY
*> \verbatim
*>          INFO_ARRAY is INTEGER array, dimension( 9 )
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
*>          INFO_SGETRF(4)
*>              = value of INFO as specified by FLAG_REPORT(1) above
*>          INFO_ARRAY(5)
*>              = 1 = number of input/output arguments reported on
*>          INFO_ARRAY(6)
*>              = 2 = number of internal LAPACK calls reported on
*>          INFO_ARRAY(7) reports exceptions in A, as specified by 
*>            FLAG_REPORT
*>              = -1 if not checked (default)
*>              =  0 if checked and contains no Infs or NaNs
*>              =  1 if checked and contains an Inf or NaN on input
*>                   but not output
*>              =  2 if checked and contains an Inf or NaN on output 
*>                   but not input
*>              =  3 if checked and contains an Inf or NaN on input 
*>                   and output
*>            If INFO_ARRAY(7) = 0 or 1 on input, then A will not
*>            be checked again. Input values < -1 or > 1 will be 
*>            treated the same as -1, i.e. not checked.
*>          INFO_ARRAY(8) reports exceptions in the first
*>            call to SGETRF2_EC, as specified by FLAG_REPORT
*>              = -1 if not checked (default)
*>              =  0 if checked and no Infs or NaNs reported
*>              =  1 if checked and no input or output contains an Inf 
*>                   or NaN, but some LAPACK call deeper in the call 
*>                   chain signaled an Inf or NaN
*>              =  2 if checked and an input argument contains an Inf 
*>                   or NaN, but not an output
*>              =  3 if checked and an output argument contains an Inf 
*>                   or NaN, but not an input
*>              =  4 if checked and an argument contains an Inf or NaN 
*>                   on input and output
*>          INFO_ARRAY(9) reports exceptions in the second
*>            call to SGETRF2_EC, analogously to INFO_ARRAY(8),
*>            reporting the maximum over all calls. 
*> \endverbatim
*>
*> \param[in] CONTEXT
*> \verbatim
*>           CONTEXT is POINTER to an "opaque object"
*> \endverbatim
*>
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
      SUBROUTINE SGETRF_EC( M, N, A, LDA, IPIV,
     $                      INFO, FLAG_REPORT, INFO_ARRAY, CONTEXT )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * )
      INTEGER            FLAG_REPORT( 2 ), INFO_ARRAY( * )
*     ..
*     .. Pointer Arguments
      POINTER            CONTEXT
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      CHARACTER, DIMENSION(6), PARAMETER :: 
     $   ROUTINE_NAME = (/ 'S','G','E','T','R','F' /)
*     ..
*     .. Local Scalars ..
      LOGICAL            CALL_REPORT_EXCEPTIONS
      INTEGER            I, IINFO, J, JB, NB
      INTEGER            WHAT, HOW
*     ..
*     .. Local Arrays ..
      INTEGER            FLAG_REPORT_INTERNAL(2), FLAG_REPORT_CALL(2)
      INTEGER            INFO_SGETRF2_TMP1(9), INFO_SGETRF2_TMP2(9)
      INTEGER            INFO_INTERNAL(2)
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM, SGETRF2_EC, SLASWP, STRSM, XERBLA
      EXTERNAL           CHECKINIT1, CHECKINIT2
      EXTERNAL           SGECHECKARG, CHECKCALL
      EXTERNAL           UPDATE_INFO
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
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
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
*
*     Initialize error flags
      CALL CHECKINIT2( FLAG_REPORT_INTERNAL, INFO, INFO_INTERNAL, 
     $                 INFO_ARRAY, 1, 2)  
*
      IF( INFO.NE.0 ) THEN
         IF (CALL_REPORT_EXCEPTIONS) 
     $     CALL REPORT_EXCEPTIONS(CONTEXT,6,ROUTINE_NAME,INFO_ARRAY)
         CALL XERBLA( 'SGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Check for exceptional inputs in A
      CALL SGECHECKARG( FLAG_REPORT_INTERNAL, M, N, A, LDA, 
     $                  INFO, INFO_INTERNAL, INFO_ARRAY, 3, 2, 0, 7 )
100   CONTINUE
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'SGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
*        Indicate if input already checked for Infs and NaNs
         IF (WHAT .GE. 1 .AND. HOW .GE. 1) THEN
            INFO_SGETRF2_TMP1(7) = -1
            IF (INFO_ARRAY(7) .NE. 0)
     $          INFO_SGETRF2_TMP1(7) = INFO_ARRAY(7)
         ENDIF
*
         CALL SGETRF2_EC( M, N, A, LDA, IPIV, 
     $             INFO, FLAG_REPORT_CALL, INFO_SGETRF2_TMP1, CONTEXT)
*   
*        Check inputs, outputs and internal calls of
*        first call to SGETRF2_EC
         CALL CHECKCALL( FLAG_REPORT_INTERNAL, INFO_INTERNAL,
     $                   INFO_SGETRF2_TMP1, INFO_ARRAY, N+2, 8 )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
*           Indicate that input not already checked
            IF (WHAT .GE. 1 .AND. HOW .GE. 1) 
     $         INFO_SGETRF2_TMP2(7) = -1
*
            CALL SGETRF2_EC( M-J+1, JB, A(J,J), LDA, IPIV(J), 
     $            IINFO, FLAG_REPORT_CALL, INFO_SGETRF2_TMP2, CONTEXT)
*   
*           Check inputs, outputs and internal calls of
*           second call to SGETRF2_EC
            CALL CHECKCALL( FLAG_REPORT_INTERNAL, INFO_INTERNAL,
     $                      INFO_SGETRF2_TMP2, INFO_ARRAY, N+3, 9 )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL SLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL SLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit',JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A(J, J+JB),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL SGEMM( 'No transpose','No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
*
*     Check for errors before returning
*
      IF (WHAT.EQ.-1) RETURN
*     Check for exceptional outputs in A
      CALL SGECHECKARG(FLAG_REPORT_INTERNAL, M, N, A, LDA, 
     $                 INFO, INFO_INTERNAL, INFO_ARRAY, 3, 3, N+1, 7)
*
*     Update INFO and INFO_ARRAY
      CALL UPDATE_INFO( INFO, INFO_ARRAY, INFO_INTERNAL )
      IF (CALL_REPORT_EXCEPTIONS .AND. INFO .NE. 0) 
     $   CALL REPORT_EXCEPTIONS(CONTEXT, 6, ROUTINE_NAME, INFO_ARRAY)
      RETURN
*
*     End of SGETRF_EC
*
      END
