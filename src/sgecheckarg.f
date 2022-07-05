*> \brief \b SGECHECKARG
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SGECHECKARG + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgecheckarg.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgecheckarg.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgecheckarg.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGECHECKARG( FLAG_REPORT_INTERNAL, M, N, A, LDA, 
*      $                        INFO, INFO_INTERNAL, INFO_ARRAY,
*      $                        ARGNUM, INOUT, ERRFLAG, LOC )
*
*       .. Scalar Arguments ..
*       INTEGER            M, N, LDA, ARGNUM, INOUT, INFO, ERRFLAG 
*       INTEGER            LOC
*       ..
*       .. Array Arguments ..
*       REAL               A( LDA, * )
*       INTEGER            FLAG_REPORT_INTERNAL( 2 ) 
*       INTEGER            INFO_INTERNAL( 2 ), INFO_ARRAY( * )
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGECHECKARG checks whether argument A contains Infs or NaNs,
*> depending on the input variable FLAG_REPORT_INTERNAL, and 
*> accordingly updates reporting variables INFO_INTERNAL and 
*> INFO_ARRAY.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] FLAG_REPORT_INTERNAL
*> \verbatim
*>          FLAG_REPORT_INTERNAL is INTEGER array, dimension(2)
*>          FLAG_REPORT_INTERNAL(1) defines what kinds of exceptions 
*>          to report, using INFO and possibly also INFO_ARRAY for 
*>          more details.
*>          FLAG_REPORT_INTERNAL(1)
*>           <= -1 turns off all error checking
*>            =  0 standard error checks only (eg LDA < 0)
*>            =  1 also check inputs and outputs for Infs and NaN
*>           >=  2 also check input and output arguments of internal 
*>               LAPACK routines (not performed unless 
*>               FLAG_REPORT_INTERNAL(2) >= 1)
*>           FLAG_REPORT_INTERNAL(2) determines how errors are 
*>           reported:
*>           FLAG_REPORT_INTERNAL(2)
*>           <= 0 only returns INFO
*>            = 1 also return INFO_ARRAY with more details
*>            = 2 also calls REPORT_EXCEPTIONS if there are errors
*>           >= 3 also calls REPORT_EXCEPTIONS in internal LAPACK
*>                routines
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows in A. M >= 0
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns in A. M >= 0
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          Matrix to be checked for containing Infs and NaNs
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in] INFO
*> \verbatim
*>          INFO is INTEGER
*>          INFO contains the standard value that LAPACK would return
*> \endverbatim
*>
*> \param[in,out] INFO_INTERNAL
*> \verbatim
*>          INFO_INTERNAL is INTEGER array, dimension(2)
*>          INFO_INTERNAL( 1 ) is used to track potential changes
*>          to INFO from errors detected by checking input and output
*>          arguments.
*>          INFO_INTERNAL( 2 ) is not accessed.
*> \endverbatim
*>
*> \param[in,out] INFO_ARRAY
*> \verbatim
*>          INFO_ARRAY is INTEGER ARRAY
*>          INFO_ARRAY can report error checks for each
*>          floating point input and output, and each 
*>          internal subroutine call.
*> \endverbatim
*>
*> \param[in] ARGNUM
*> \verbatim
*>          ARGNUM is INTEGER
*>          A is ARGNUM-th argument of routine calling SGECHECKARG
*> \endverbatim
*>
*> \param[in] INOUT
*> \verbatim
*>          INOUT is INTEGER
*>          INOUT = 0 if A is an input-only argument
*>                    in the routine calling SGECHECKARG.
*>          INOUT = 1 if A is an output-only argument.
*>          INOUT = 2 if A is both an input and output argument,
*>                    and is being checked on input.
*>          INOUT = 3 if A is both an input and output argument,
*>                    and is being checked on output.
*> \endverbatim
*>
*> \param[in] ERRFLAG
*> \verbatim
*>          ERRFLAG is INTEGER
*>          INFO_INTERNAL(1) = ERRFLAG is used to indicate that
*>          A contains an Inf or NaN on output, but not input.
*>          Not accessed if A is checked on input (INOUT = 0 or 2).
*> \endverbatim
*>
*> \param[in] LOC
*> \verbatim
*>          LOC is INTEGER
*>          LOC points to the entry of INFO_ARRAY used to report on A
*> \endverbatim
*>
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
      SUBROUTINE SGECHECKARG( FLAG_REPORT_INTERNAL, M, N, A, LDA, 
     $                        INFO, INFO_INTERNAL, INFO_ARRAY,
     $                        ARGNUM, INOUT, ERRFLAG, LOC )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            M, N, LDA, ARGNUM, INOUT, INFO, ERRFLAG 
      INTEGER            LOC
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * )
      INTEGER            FLAG_REPORT_INTERNAL( 2 ) 
      INTEGER            INFO_INTERNAL( 2 ), INFO_ARRAY( * )
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            WHAT, HOW, INFO_A, INFNAN, II, JJ      
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEINFNAN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Check for exceptional entries in A 
      WHAT = FLAG_REPORT_INTERNAL( 1 )
      IF (WHAT .GE. 1) THEN
*        Check input or output for Infs and NaNs
         HOW = FLAG_REPORT_INTERNAL( 2 )
         IF (HOW .GE. 1) THEN
*           INFO_ARRAY used for reporting, as well as INFO
            INFO_A = INFO_ARRAY(LOC)
            IF (INOUT .EQ. 0 .OR. INOUT .EQ. 2) THEN
*              Checking A on input
               IF (INFO_A .NE. 0 .AND. INFO_A .NE. 1) THEN
*                 A not checked yet, so check
                  CALL SGEINFNAN( M, N, A, LDA, INFNAN, II, JJ )
                  INFO_ARRAY(LOC) = INFNAN
               ENDIF
*              Update INFO_INTERNAL to point to A if not already set
               IF (INFO_INTERNAL(1).EQ.0 .AND. INFO_ARRAY(LOC).EQ.1)
     $            INFO_INTERNAL(1) = -ARGNUM
            ELSEIF (INOUT .EQ. 1) THEN
*              Checking output-only variable A on output
               CALL SGEINFNAN( M, N, A, LDA, INFNAN, II, JJ )
               INFO_ARRAY(LOC) = 2*INFNAN
*              Update INFO_INTERNAL to point to A if not already set
               IF (INFO_INTERNAL(1).EQ.0 .AND. INFO_ARRAY(LOC).EQ.2)
     $            INFO_INTERNAL(1) = ERRFLAG
            ELSE
*              Checking input-output variable A on output
               CALL SGEINFNAN( M, N, A, LDA, INFNAN, II, JJ )
               IF (INFNAN .EQ. 1) INFO_ARRAY(LOC) = INFO_ARRAY(LOC)+2
*              Update INFO_INTERNAL to point to A if not already set
               IF (INFO_INTERNAL(1).EQ.0 .AND. INFO_ARRAY(LOC).EQ.2)
     $            INFO_INTERNAL(1) = ERRFLAG
            ENDIF
         ELSEIF (INFO_INTERNAL(1) .EQ. 0 .AND. INFO .EQ. 0) THEN
*           Only INFO used for reporting, not set yet
            CALL SGEINFNAN( M, N, A, LDA, INFNAN, II, JJ )
            IF (INOUT .EQ. 0 .OR. INOUT .EQ. 2) THEN
*              Checking A on input
               IF (INFNAN .EQ. 1) INFO_INTERNAL(1) = -ARGNUM
            ELSE
*              Checking A on output
               IF (INFNAN .EQ. 1) INFO_INTERNAL(1) = ERRFLAG
            ENDIF
         ENDIF
      ENDIF
*     
      RETURN
*
*     End of SGECHECKARG
*
      END
