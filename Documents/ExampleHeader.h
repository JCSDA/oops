/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef DOCUMENTS_EXAMPLEHEADER_H_
#define DOCUMENTS_EXAMPLEHEADER_H_

// C++ header files
    // --- Add #includes here (use forward declarations if possible) --- //

#include <boost/noncopyable.hpp>

// Forward declarations
    // --- Add declarations (class SomeClass;) for non-oops classes here --- //
namespace oops {
    // --- Add declarations (class SomeOopsClass;) for oops classes here --- //
}

namespace oops {

  /// Doxygen brief comment. One line only. The comment must start with ///
  /*!
   *  Replace this comment block with longer documentation that describes
   *  the function of the class. The comment will be processed by doxygen
   *  and included in the documentation. Member functions should be documented
   *  where they are declared (e.g. the makeTea function below).
   */

  class Example : private boost::noncopyable {
   public:

  // -- Constructors
    Example();

  // -- Destructor
    ~Example()

  // -- Methods

  /// makeTea takes an argument indicating the number of teabags to be used.
    void makeTea(const int);

   private:
  // -- Members

  // -- Overridden methods
  };

}  // namespace oops
#endif  // DOCUMENTS_EXAMPLEHEADER_H_
