
#pragma once

#include "const.h"
#include "vars.h"

YAKL_INLINE real gammafff(real x){
  return exp(lgamma(x));
}

