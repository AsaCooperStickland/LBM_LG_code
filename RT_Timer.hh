#ifndef RT_TIMER_HH
#define RT_TIMER_HH

//=======================================================
//
// = LIBRARY
//     Timer
//
// = FILENAME
//     RT_Timer.hh
//
// = AUTHOR(S)
//     Frederic Guidec (EPFL)
//     Alexandre Dupuis
//
// = VERSION
//     $Revision: 1.1 $ 
//
// = DATE RELEASED
//     $Date: 2000/05/29 14:33:41 $
//
// = COPYRIGHT
//     University of Geneva, Switzerland
//
//=======================================================


#include "Timer.hh"

class RT_Timer: public Timer
//=======================================================
//
// = DESCRIPTION
//
// A <{RT_Timer}> makes it possible to measure the real time elapsed
// between events during an execution. This time is returned in
// <{seconds}>.  The class <{CPU_Timer}> inherits most of its features
// from the class <{Timer}>.  See the manual page of this class for a
// detailed report on timers' interface.
//
// = SEE ALSO
//   <{Timer}>, <{CPU_Timer}>
//
// = INHERITS
//     <{Timer}>
//
// = NOTES 
//
// A <{RT_Timer}> measures <[real time]>, as returned by the local
// clock of a processing unit. It does <[not]> measure CPU time.
//=======================================================
{
private:

    virtual double currentTime();
    // Returns the current <[real]> time. (Implements the <{pure
    // virtual}> method inherited from class <{Timer}>.)

public:

  //=======================================================
  // = CONSTRUCTOR
  //=======================================================
   
  RT_Timer();
  // Constructor.
  
  RT_Timer(const RT_Timer &other);
  // Copy constructor.

  
  //=======================================================
  // = DESTRUCTOR
  //=======================================================

  virtual ~RT_Timer();
  // Destructor.

  //=======================================================
  // = OPERATORS
  //=======================================================
  
  RT_Timer& operator=(const RT_Timer& other);
  // Assignment operator.

};

#endif

