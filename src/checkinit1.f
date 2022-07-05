*> \brief \b CHECKINIT1
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download CHECKINIT1 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/checkinit1.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/checkinit1.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/checkinit1.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE CHECKINIT1( FLAG_REPORT, FLAG_REPORT_INTERNAL, 
*                              FLAG_REPORT_CALL, 
*                              CALL_REPORT_EXCEPTIONS, CONTEXT )
*
*       .. Scalar Arguments ..
*       LOGICAL            CALL_REPORT_EXCEPTIONS 
*       ..
*       .. Array Arguments ..
*       INTEGER            FLAG_REPORT( 2 ), FLAG_REPORT_INTERNAL( 2 )
*       INTEGER            FLAG_REPORT_CALL( 2 )
*       ..
*       .. Pointer Arguments ..
*       POINTER            CONTEXT ... advice requested
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CHECKINIT1 initializes variables used to determine how
*> to check and report errors.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] FLAG_REPORT
*> \verbatim
*>          FLAG_REPORT is INTEGER array, dimension(2)
*>          FLAG_REPORT(1) defines what kinds of exceptions to report,
*>          using INFO and possibly also INFO_ARRAY for more details.
*>          FLAG_REPORT(1)
*>           <= -1 turns off all error checking
*>            =  0 standard error checks only (eg LDA < 0)
*>            =  1 also check inputs and outputs for Infs and NaNs
*>           >=  2 also check input and output arguments of internal 
*>               LAPACK routines (not performed unless 
*>               FLAG_REPORT(2) >= 1)
*>           FLAG_REPORT(2) determines how errors are reported:
*>           FLAG_REPORT(2)
*>           <= 0 only returns INFO
*>            = 1 also returns INFO_ARRAY with more details
*>            = 2 also calls REPORT_EXCEPTIONS if there are errors
*>            = 3 also calls REPORT_EXCEPTIONS in internal LAPACK
*>                routines
*>           >= 4 will call GET_FLAGS_TO_REPORT(CONTEXT,FLAG_REPORT)
*>                to get values of FLAG_REPORT to use, overriding
*>                input values. The user needs to have called 
*>                SET_FLAGS_TO_REPORT(CONTEXT,FLAG_REPORT) 
*>                before calling SCHECKINIT in order to set 
*>                FLAG_REPORT, otherwise the default is 
*>                FLAG_REPORT = [0, 0]. The input array FLAG_REPORT 
*>                will not be overwritten.
*> \endverbatim
*>
*> \param[out] FLAG_REPORT_INTERNAL
*> \verbatim
*>          FLAG_REPORT_INTERNAL is INTEGER array, dimension(2)
*>          FLAG_REPORT_INTERNAL contains updated values of
*>          FLAG_REPORT on return.
*> \endverbatim
*>
*> \param[in] FLAG_REPORT_CALL
*> \verbatim
*>          FLAG_REPORT_CALL is INTEGER array, dimension(2)
*>          FLAG_REPORT_CALL contains the values of FLAG_REPORT
*>          to be used in internal LAPACK calls.
*>
*> \param[out] CALL_REPORT_EXCEPTIONS
*> \verbatim
*>          CALL_REPORT_EXCEPTIONS is LOGICAL
*>          CALL_REPORT_EXCEPTIONS is true if REPORT_EXCEPTIONS should 
*>          be called, and false if not.
*> \endverbatim
*>
*> \param[in] CONTEXT
*> \verbatim
*>          CONTEXT is POINTER to an "opaque object"
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
      SUBROUTINE CHECKINIT1( FLAG_REPORT, FLAG_REPORT_INTERNAL, 
     $                       FLAG_REPORT_CALL,
     $                       CALL_REPORT_EXCEPTIONS, CONTEXT )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            CALL_REPORT_EXCEPTIONS 
*     ..
*     .. Array Arguments ..
      INTEGER            FLAG_REPORT( 2 ), FLAG_REPORT_INTERNAL( 2 )
      INTEGER            FLAG_REPORT_CALL( 2 )
*     ..
*     .. Pointer Arguments ..
      POINTER            CONTEXT 
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER, DIMENSION(4), PARAMETER ::
     $   WHAT_NEXT = (/ -1, 0, 0, 2 /)
      INTEGER, DIMENSION(4), PARAMETER ::
     $   HOW_NEXT = (/ 0, 1, 1, 3 /)
*
*     .. Local Scalars ..
      INTEGER            WHAT, HOW
*
*     .. Local Arrays ..
      INTEGER            FLAG_REPORT_TMP( 2 )       
*     ..
*     .. External Subroutines ..
      EXTERNAL           GET_FLAGS_TO_REPORT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     WHAT records what errors are to be reported
      WHAT = MAX( -1, MIN( FLAG_REPORT(1), 2 ) )
*     HOW records how those errors are to be reported
      HOW = 0
      CALL_REPORT_EXCEPTIONS = .FALSE.
      FLAG_REPORT_INTERNAL( 1 ) = WHAT
      FLAG_REPORT_INTERNAL( 2 ) = HOW
      FLAG_REPORT_CALL( 1 ) = WHAT_NEXT(WHAT+2)
      FLAG_REPORT_CALL( 2 ) = HOW_NEXT(HOW+1)
*     Check if error reporting turned off
      IF (WHAT .EQ. -1) RETURN
      HOW = MAX( 0, MIN(FLAG_REPORT(2), 4 ) )
      FLAG_REPORT_INTERNAL( 2 ) = HOW
      IF (HOW .EQ. 4) THEN
*        Get updated FLAG_REPORT
         CALL GET_FLAGS_TO_REPORT(CONTEXT, FLAG_REPORT_TMP)
         WHAT = MAX( -1, MIN( FLAG_REPORT_TMP(1), 2 ) )
         FLAG_REPORT_INTERNAL( 1 ) = WHAT
         FLAG_REPORT_CALL( 1 ) = WHAT_NEXT(WHAT+2)
         IF (WHAT .EQ. -1) THEN
*          Error reporting turned off
           HOW = 0
           FLAG_REPORT_INTERNAL( 2 ) = HOW
           FLAG_REPORT_CALL( 2 ) = HOW_NEXT(HOW+1)
           RETURN
         END IF
         HOW = MAX( 0, MIN(FLAG_REPORT_TMP(2), 3 ) )
         FLAG_REPORT_INTERNAL( 2 ) = HOW
         FLAG_REPORT_CALL( 2 ) = HOW_NEXT(HOW+1)
      END IF
*     Error reporting not turned off
*     Decide whether to call REPORT_EXCEPTIONS
      IF (HOW .GE. 2) CALL_REPORT_EXCEPTIONS = .TRUE.
*     
      RETURN
*
*     End of CHECKINIT1
*
      END