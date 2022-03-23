// 
// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date February 2019
//
// Tests more comprehensive runs of minus using the public interface
// 
#include <cstring>
#include <iostream>
#include <p2pt/p2pt.h>

using namespace p2pt;

static void
test_hello()
{
  p2pt<double>::hello();
  TEST("Rodou hello? ", true, true);
}

void
test_p2pt()
{
  test_hello();
}

TESTMAIN(test_p2pt);
