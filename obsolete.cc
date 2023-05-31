#include <iostream>

#include "include/lcde/obsolete.h"
#include "test/rdata.h"

int main() {
  lcde::OLCDE o(rdata);
  
  o.build();
  o.printLCD();
  o.plot("all", false);
}