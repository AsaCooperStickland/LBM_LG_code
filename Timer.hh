#ifndef TIMER_HH
#define TIMER_HH

//=======================================================
//
// = LIBRARY
//     Timer
//
// = FILENAME
//     Timer.hh
//
// = AUTHOR(S)
//     Frederic Guidec (EPFL)
//     Alexandre Dupuis
//
// = VERSION
//     $Revision: 1.1 $
//
// = DATE RELEASED
//     $Date: 2000/05/29 14:33:52 $
//
// = COPYRIGHT
//     University of Geneva, Switzerland
//
//=======================================================

#include "Bool.hh"

class Timer
//=======================================================
//
// = DESCRIPTION
//
//   A <{Timer}> makes it possible to measure the time elapsed between
//   two events during an execution. This time is returned in
//   <{seconds}>.
//
//
// = SEE ALSO
//   <{CPU_Timer}>, <{RT_Timer}>
//
//=======================================================
{
private:

  bool isRunning_;
  double startTime_;
  double stopTime_;
  double elapsed_;

  virtual double currentTime() = 0;
  // Returns the current time. This is a <{pure virtual}> method. It
  // must be defined in a descendant class.

public:

  //=======================================================
  // = CONSTRUCTOR
  //=======================================================
   
  Timer();
  // Constructor.
  
  Timer(const Timer &other);
  // Copy constructor.

  //=======================================================
  // = DESTRUCTOR
  //=======================================================

  virtual ~Timer();
  // Destructor.

  //=======================================================
  // = OPERATORS
  //=======================================================
   
  Timer& operator=(const Timer& other);
  // Assignment operator.

  //=======================================================
  // = ACCESSORS
  //=======================================================
   
  double elapsed();
  // Returns the time elapsed while this timer was running, since the
  // last call to <{reset()}> (or if no such call was made, since the
  // creation of this timer). Note that a call to this routine does
  // <{not}> stop the timer.

  bool isRunning();
  // Is this timer currently running?

  //=======================================================
  // = MUTATORS
  //=======================================================

  void reset();
  // Resets this timer.

  void start();
  // Starts this timer.

  void stop();
  // Stops this timer.

};

#endif

