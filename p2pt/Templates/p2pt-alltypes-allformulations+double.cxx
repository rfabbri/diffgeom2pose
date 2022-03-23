// Helps instantiate and reuse code,
// when including the .hxx directly causes inefficiencies and slowdown
#include <p2pt/p2pt.hxx>

namespace p2pt {
template class p2pt<double>;
}
