*> \brief <b> SGESV_EC computes the solution to system of linear equations A * X = B for GE matrices</b> (simple driver)
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SGESV_EC + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgesv_ec.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgesv_ec.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgesv_ec.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGESV_EC( N, NRHS, A, LDA, IPIV, B, LDB, 
*      $                     INFO, FLAG_REPORT, INFO_ARRAY, CONTEXT )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       INTEGER            FLAG_REPORT( 2 ), INFO_ARRAY( * )
*       REAL               A( LDA, * ), B( LDB, * )
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
*> SGESV_EC computes the solution to a real system of linear equations
*>    A * X = B,
*> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
*>
*> The LU decomposition with partial pivoting and row interchanges is
*> used to factor A as
*>    A = P * L * U,
*> where P is a permutation matrix, L is unit lower triangular, and U 
*> is upper triangular.  The factored form of A is then used to solve
*> the system of equations A * X = B.
*>
*> SGESV_EC also provides new exception handling and 
*> reporting capabilities.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of linear equations, i.e., the order of the
*>          matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of columns
*>          of the matrix B.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          On entry, the N-by-N coefficient matrix A.
*>          On exit, the factors L and U from the factorization
*>          A = P*L*U; the unit diagonal elements of L are not stored.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          The pivot indices that define the permutation matrix P;
*>          row i of the matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is REAL array, dimension (LDB,NRHS)
*>          On entry, the N-by-NRHS matrix of right hand side matrix B.
*>          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
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
*>                 INFO = -6 means B contained an Inf or NaN on 
*>                           input, but not A; execution continues.
*>                 INFO = N+1 means that A contained an Inf or NaN on 
*>                        output but neither A nor B did on input.
*>                 INFO = N+2 means that B contained an Inf or NaN on 
*>                        output but neither A nor B did on input nor
*>                        A on output.
*>                 A and B will also be checked on input and output if 
*>                 FLAG_REPORT(2) = 1, 2 or 3 and reported in 
*>                 INFO_ARRAY as described below.
*>            >= 2 also checks for Infs and NaNs as inputs or outputs 
*>                 in all internal LAPACK calls in the call chain, if 
*>                 INFO is not already nonzero, and 
*>                 FLAG_REPORT(2) = 1, 2 or 3. In this case: 
*>                 INFO = N+3 means that either the call to
*>                        SGETRF_EC had an Inf or NaN as an 
*>                        input or output as above, or a subroutine
*>                        in its call tree did.
*>                 INFO = N+4 means that the call to SGETRS_EC
*>                        had an Inf or NaN as an 
*>                        input or output as above.
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
*>              =  2 means that SGESV_EC will also call 
*>                   REPORT_EXCEPTIONS to report INFO_ARRAY, if INFO 
*>                   is nonzero.
*>              =  3 means that all calls in the call tree will also
*>                   call REPORT_EXCEPTIONS, if the value of INFO they
*>                   return is nonzero.
*>             >=  4 means that SGESV_EC will call 
*>                   GET_FLAGS_TO_REPORT(CONTEXT,FLAG_REPORT) 
*>                   to get values of FLAG_REPORT to use, overriding 
*>                   input values. The user needs to have called 
*>                   SET_FLAGS_TO_REPORT(CONTEXT,FLAG_REPORT) 
*>                   before calling SGESV_EC in order to set 
*>                   FLAG_REPORT, otherwise the default is 
*>                   FLAG_REPORT = [0, 0]. The input array FLAG_REPORT 
*>                   will not be overwritten.
*> \endverbatim
*>
*> \param[out] INFO_ARRAY
*> \verbatim
*>          INFO_ARRAY is INTEGER array, dimension( 10 )
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
*>              = 2 = number of internal LAPACK calls reported on
*>          INFO_ARRAY(7), reports exceptions in A, as specified by 
*>            FLAG_REPORT
*>              = -1 if not checked (default)
*>              =  0 if checked and contains no Infs or NaNs
*>              =  1 if checked and contains an Inf or NaN on input
*>                   but not output
*>              =  2 if checked and contains an Inf or NaN on output 
*>                   but not input
*>              =  3 if checked and contains an Inf or NaN on input 
*>                   and output
*>          INFO_ARRAY(8), reports exceptions in B, analogously to 
*>            INFO_ARRAY(7)
*>          INFO_ARRAY(9), reports exceptions in SGETRF_EC, as 
*>            specified by FLAG_REPORT
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
*>          INFO_ARRAY(10), reports exceptions in SGETRS_EC, 
*>            analogously to INFO_ARRAY(9)
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
*> \ingroup realGEsolve
*
*  =====================================================================
      SUBROUTINE SGESV_EC( N, NRHS, A, LDA, IPIV, B, LDB, 
     $                     INFO, FLAG_REPORT, INFO_ARRAY, CONTEXT ) 
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), B( LDB, * )
      INTEGER            FLAG_REPORT( 2 ), INFO_ARRAY( * )
*     ..
*     .. Pointer Arguments
      POINTER            CONTEXT
*
*  =====================================================================
*
*     .. Parameters ..
      CHARACTER, DIMENSION(5), PARAMETER :: 
     $   ROUTINE_NAME = (/ 'S','G','E','S','V' /)
*     ..
*     .. Local Scalars ..
      LOGICAL            CALL_REPORT_EXCEPTIONS
      INTEGER            WHAT, HOW
*     ..
*     .. Local Arrays ..
      INTEGER            FLAG_REPORT_INTERNAL(2), FLAG_REPORT_CALL(2)
      INTEGER            INFO_SGETRF_TMP(9), INFO_SGETRS_TMP(9)
      INTEGER            INFO_INTERNAL(2)
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGETRF_EC, SGETRS_EC, XERBLA
      EXTERNAL           CHECKINIT1, CHECKINIT2
      EXTERNAL           SGECHECKARG, CHECKCALL
      EXTERNAL           UPDATE_INFO
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*     Initialize error checking flags
      CALL CHECKINIT1(FLAG_REPORT, FLAG_REPORT_INTERNAL, 
     $                FLAG_REPORT_NEXT, CALL_REPORT_EXCEPTIONS, 
     $                CONTEXT ) 
      WHAT = FLAG_REPORT_INTERNAL( 1 )
      HOW = FLAG_REPORT_INTERNAL( 2 )
      IF (WHAT .EQ. -1 ) GOTO 100
*
*     Check for standard input errors
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
*
*     Initialize error flags
      CALL CHECKINIT2( FLAG_REPORT_INTERNAL, INFO, INFO_INTERNAL, 
     $                 INFO_ARRAY, 2, 2) 
*
      IF( INFO.NE.0 ) THEN
         IF (CALL_REPORT_EXCEPTIONS) 
     $       CALL REPORT_EXCEPTIONS(CONTEXT,5,ROUTINE_NAME,INFO_ARRAY)
         CALL XERBLA( 'SGESV ', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*     Check for exceptional inputs in A 
      CALL SGECHECKARG( FLAG_REPORT_INTERNAL, N, N, A, LDA, 
     $                  INFO, INFO_INTERNAL, INFO_ARRAY, 3, 2, 0, 7 )
*     Check for exceptional inputs in B 
      CALL SGECHECKARG( FLAG_REPORT_INTERNAL, N, NRHS, B, LDB, 
     $                  INFO, INFO_INTERNAL, INFO_ARRAY, 6, 2, 0, 8 )
*
100   CONTINUE
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*     Compute the LU factorization of A.
*
*     Indicate if input A already checked for Infs and NaNs
      IF (WHAT .GE. 1 .AND. HOW .GE. 1) 
     $    INFO_SGETRF_TMP(7) = INFO_ARRAY(7)
      CALL SGETRF_EC( N, N, A, LDA, IPIV, 
     $               INFO, FLAG_REPORT_CALL, INFO_SGETRF_TMP, CONTEXT)
*
*     Check inputs, outputs and internal calls of SGETRF_EC
      CALL CHECKCALL( FLAG_REPORT_INTERNAL, INFO_INTERNAL,
     $                INFO_SGETRF_TMP, INFO_ARRAY, N+3, 9 )
      IF( INFO.EQ.0 ) THEN
*
*        Solve the system A*X = B, overwriting B with X.
*
*        Indicate if input B already checked for Infs and NaNs
         IF (WHAT .GE. 1 .AND. HOW .GE. 1) 
     $      INFO_SGETRS_TMP(8) = INFO_ARRAY(8)
         CALL SGETRS_EC( 'No transpose', N, NRHS, A, LDA, IPIV, B,LDB,
     $                  INFO,FLAG_REPORT_CALL,INFO_SGETRS_TMP,CONTEXT)
*
*        Check inputs, outputs and internal calls of call to SGETRS_EC
         CALL CHECKCALL( FLAG_REPORT_INTERNAL, INFO_INTERNAL,
     $                   INFO_SGETRS_TMP, INFO_ARRAY, N+4, 9 )
      ENDIF
*
*     Check for errors before returning
*
      IF (WHAT.EQ.-1) RETURN
*     Check for exceptional outputs in A
      CALL SGECHECKARG(FLAG_REPORT_INTERNAL, N, N, A, LDA, 
     $                 INFO, INFO_INTERNAL, INFO_ARRAY, 3, 3, N+1, 7)
*     Check for exceptional outputs in B 
      CALL SGECHECKARG(FLAG_REPORT_INTERNAL, N, NRHS, B, LDB, 
     $                 INFO, INFO_INTERNAL, INFO_ARRAY, 7, 3, N+2, 8)
*
*     Update INFO and INFO_ARRAY
      CALL UPDATE_INFO( INFO, INFO_ARRAY, INFO_INTERNAL )
      IF (CALL_REPORT_EXCEPTIONS .AND. INFO .NE. 0) 
     $   CALL REPORT_EXCEPTIONS(CONTEXT, 5, ROUTINE_NAME, INFO_ARRAY)
      RETURN
*
*     End of SGESV_EC
*
      END
