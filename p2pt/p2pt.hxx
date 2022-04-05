#ifndef p2pt_hxx_
#define p2pt_hxx_

#include <iostream>
#include "p2pt.h"
#include "rf_find_bounded_root_intervals.hxx"
#include "rf_sample_pose_poly.hxx"

namespace P2Pt {
  
template <typename F> 
void p2pt<F>::
hello()
{
  std::cout << "hello\n";
}
  
} // namespace p2pt


#endif // p2pt_hxx_

