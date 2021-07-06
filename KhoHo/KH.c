/*-*- compile-command: "cc -c -o ./KH.o -g -O3 -Wall -fomit-frame-pointer -fno-strict-aliasing -fPIC -I"/usr/include/x86_64-linux-gnu" ./KH.c && cc -o ./KH.so -shared -g -O3 -Wall -fomit-frame-pointer -fno-strict-aliasing -fPIC -Wl,-shared -Wl,-z,relro ./KH.o -lc -lm -L/usr/lib/x86_64-linux-gnu -lpari"; -*-*/
#include <pari/pari.h>
/*
GP;install("init_KH","v","init_KH","././KH.so");
*/
void init_KH(void);
/*End of prototype*/

static GEN KhoHo;
static GEN KTable_Rolfsen;
static GEN KTable_11;
static GEN LTable_11;
/*End of global vars*/

void
init_KH(void)	  /* void */
{
  KhoHo = pol_x(fetch_user_var("KhoHo"));
  KTable_Rolfsen = pol_x(fetch_user_var("KTable_Rolfsen"));
  KTable_11 = pol_x(fetch_user_var("KTable_11"));
  LTable_11 = pol_x(fetch_user_var("LTable_11"));
  /*
  * Read the main program as well as lists of some knots and links.
  */
  
  gp_read_file(GENtostr_unquoted(KhoHo));
  gp_read_file(GENtostr_unquoted(KTable_Rolfsen));
  gp_read_file(GENtostr_unquoted(KTable_11));
  gp_read_file(GENtostr_unquoted(LTable_11));
  gp_read_file("KTable_12a.gz");
  gp_read_file("KTable_12n.gz");
  gp_read_file("KTable_13a.gz");
  gp_read_file("KTable_13n.gz");
  gp_read_file("KTable_14a.gz");
  gp_read_file("KTable_14n.gz");
  return;
}

