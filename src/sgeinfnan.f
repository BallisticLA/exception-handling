*> \brief \b SGEINFNAN
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SGEINFNAN + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeinfnan.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeinfnan.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeinfnan.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGEINFNAN( M, N, A, LDA, INFNAN, I, J )
*
*       .. Scalar Arguments ..
*       INTEGER            LDA, M, N, INFNAN, I, J
*       ..
*       .. Array Arguments ..
*       REAL               A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGEINFNAN  returns INFNAN = 1 if the matrix A contains an Inf or NaN, and 0 otherwise,
*>            The location of the Inf or NaN is at A(I,J)
*> \endverbatim
*>
*> \return SGEINFNAN
*
*  Arguments:
*  ==========
*
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.  When M = 0,
*>          INFNAN is set to 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.  When N = 0,
*>          INFNAN is set to 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          The M by N matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(M,1).
*> \endverbatim
*>
*> \param[out] INFNAN
*> \verbatim
*>          INFNAN is INTEGER
*>          INFNAN = 1 if A contains an Inf or NaN, and 0 otherwise.
*> \endverbatim
*>
*> \param[out] I
*> \verbatim
*>          I is INTEGER
*>          If A contains an Inf or NaN, one is located at A(I,J),
*>          otherwise I = 0.
*> \endverbatim
*>
*> \param[out] J
*> \verbatim
*>          J is INTEGER
*>          If A contains an Inf or NaN, one is located at A(I,J),
*>          otherwise J = 0. 
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup realGEauxiliary
*
*  =====================================================================
      SUBROUTINE SGEINFNAN( M, N, A, LDA, INFNAN, I, J )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            LDA, M, N, INFNAN, I, J
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            MAXEXP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          EXPONENT, HUGE
*     ..
*     .. Executable Statements ..
*
      INFNAN = 0
      I = 0
      J = 0
      IF( MIN( M, N ).EQ.0 ) RETURN
*     MAXEXP is the exponent of an Inf or NaN
      MAXEXP = EXPONENT(HUGE(ONE)) + 1
      DO 20 J = 1, N
         DO 10 I = 1, M
*           It could be faster to use the function IEEE_IS_FINITE 
*           provided by the intrinsic module IEEE_ARITHMETIC, but 
*           whether this module is provided is processor dependent, 
*           according to the Fortran 2008 standard.
            IF( EXPONENT(A(I,J)) .EQ. MAXEXP ) THEN
*              A(I,J) is an Inf or NaN
               INFNAN = 1
               RETURN
            ENDIF
   10    CONTINUE
   20 CONTINUE
      I = 0
      J = 0
      RETURN
*
*     End of SGEINFNAN
*
      END
