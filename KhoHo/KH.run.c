/*-*- compile-command: "cc -c -o ./KH.run.o -g -O3 -Wall -fomit-frame-pointer -fno-strict-aliasing -fPIC -I"/usr/include/x86_64-linux-gnu" ./KH.run.c && cc -o ./KH.run.so -shared -g -O3 -Wall -fomit-frame-pointer -fno-strict-aliasing -fPIC -Wl,-shared -Wl,-z,relro ./KH.run.o -lc -lm -L/usr/lib/x86_64-linux-gnu -lpari"; -*-*/
#include <pari/pari.h>
/*
GP;install("init_KH","v","init_KH","././KH.so");
*/
extern void init_KH();
/*End of prototype*/

void
init_KH(void)	  /* void */
{
  return;
}

