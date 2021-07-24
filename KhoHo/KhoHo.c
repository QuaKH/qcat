/*-*- compile-command: "cc -c -o ./KhoHo.o -g -O3 -Wall -fomit-frame-pointer -fno-strict-aliasing -fPIC -I"/usr/include/x86_64-linux-gnu" ./KhoHo.c && cc -o ./KhoHo.so -shared -g -O3 -Wall -fomit-frame-pointer -fno-strict-aliasing -fPIC -Wl,-shared -Wl,-z,relro ./KhoHo.o -lc -lm -L/usr/lib/x86_64-linux-gnu -lpari"; -*-*/
#include "Pari64-2-13-2/include/pari/pari.h"
/*
GP;install("init_KhoHo","v","init_KhoHo","././KhoHo.so");
GP;install("T_ranks_assign","vD0,G,D0,G,","T_ranks_assign","././KhoHo.so");
GP;install("D_ranks_tors","D0,G,D0,G,D0,G,D0,G,D0,G,","D_ranks_tors","././KhoHo.so");
GP;install("D_inv_factors","D0,G,D0,G,D0,G,","D_inv_factors","././KhoHo.so");
GP;install("Betti","D0,G,","Betti","././KhoHo.so");
GP;install("Torsion","D0,G,","Torsion","././KhoHo.so");
GP;install("matrix2pol","D0,G,D0,G,DGp","matrix2pol","././KhoHo.so");
GP;install("KhPol_Q","D0,G,DGp","KhPol_Q","././KhoHo.so");
GP;install("KhPol_T","D0,G,DGp","KhPol_T","././KhoHo.so");
GP;install("KhPol","D0,G,DGDGp","KhPol","././KhoHo.so");
GP;install("conj1_factor","D0,G,p","conj1_factor","././KhoHo.so");
GP;install("check_conj1","D0,G,DGp","check_conj1","././KhoHo.so");
GP;install("check_H_thick","D0,G,DGp","check_H_thick","././KhoHo.so");
GP;install("check_T_thick","D0,G,D0,G,DGp","check_T_thick","././KhoHo.so");
GP;install("pol_diags","D0,G,DGp","pol_diags","././KhoHo.so");
*/
void init_KhoHo(void);
void T_ranks_assign(GEN D_ID, GEN t_orders);
GEN D_ranks_tors(GEN datapos, GEN i, GEN j, GEN do_rank, GEN do_torsion);
GEN D_inv_factors(GEN D_ID, GEN do_rank, GEN do_torsion);
GEN Betti(GEN D_ID);
GEN Torsion(GEN D_ID);
GEN matrix2pol(GEN D_ID, GEN matr_name, GEN ret_vector, long prec);
GEN KhPol_Q(GEN D_ID, GEN ret_vector, long prec);
GEN KhPol_T(GEN D_ID, GEN ret_vector, long prec);
GEN KhPol(GEN D_ID, GEN split, GEN ret_vector, long prec);
GEN conj1_factor(GEN lmatr, long prec);
GEN check_conj1(GEN khpolQ, GEN cfactor, long prec);
GEN check_H_thick(GEN khpolQ, GEN cfactor, long prec);
GEN check_T_thick(GEN khpolQ, GEN khpolT, GEN cfactor, long prec);
GEN pol_diags(GEN khpolQ, GEN min_max, long prec);
/*End of prototype*/

static GEN H_TYPE;
static GEN MAX_DIAGRAM_NUM;
static GEN VERBOSE_LEVEL;
static GEN V_SILENT;
static GEN V_WHAT;
static GEN V_PROGRESS;
static GEN V_DEBUG;
static GEN CHECK_D2;
static GEN H_torsion_list;
static GEN KhoHo_data;
static GEN KhoHo_gvars;
static GEN KhoHo_diagr;
static GEN KhoHo_chain;
static GEN KhoHo_odd;
static GEN KhoHo_reduce;
static GEN KhoHo_sign;
static GEN KhoHo_print;
static GEN dummy;
static GEN t;
static GEN q;
static GEN Q;
/*End of global vars*/

void
init_KhoHo(void)	  /* void */
{
  KhoHo_data = pol_x(fetch_user_var("KhoHo_data"));
  KhoHo_gvars = pol_x(fetch_user_var("KhoHo_gvars"));
  KhoHo_diagr = pol_x(fetch_user_var("KhoHo_diagr"));
  KhoHo_chain = pol_x(fetch_user_var("KhoHo_chain"));
  KhoHo_odd = pol_x(fetch_user_var("KhoHo_odd"));
  KhoHo_reduce = pol_x(fetch_user_var("KhoHo_reduce"));
  KhoHo_sign = pol_x(fetch_user_var("KhoHo_sign"));
  KhoHo_print = pol_x(fetch_user_var("KhoHo_print"));
  dummy = pol_x(fetch_user_var("dummy"));
  t = pol_x(fetch_user_var("t"));
  q = pol_x(fetch_user_var("q"));
  Q = pol_x(fetch_user_var("Q"));
  /*
  *    KhoHo --- program for computing and studying Khovanov homology:
  *              main routines for computing homology and Khovanov polynomials
  *              from the chain complex and for performing tests on them.
  *
  * Copyright (C) 2002--2018 Alexander Shumakovitch <Shurik@gwu.edu>
  *
  * This program is free software; you can redistribute it and/or modify
  * it under the terms of the GNU General Public License as published by
  * the Free Software Foundation; either version 2, or (at your option)
  * any later version.
  *
  * This program  is distributed in the hope that it will be useful,
  * but WITHOUT ANY WARRANTY; without even the implied warranty of
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  * GNU General Public License for more details.
  *
  * You should have received a copy of the GNU General Public License
  * along with this program; see COPYING.gz. If not, write to the Free
  * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
  *
  *
  *    Please refer to README for more details.
  */
  
  /*
  * Type of the homology to compute:
  *    0 --> standard;  1 --> reduced;  2 --> reduced odd.
  *
  * We don't compute non-reduced odd homology, since they are uniquely
  * determined by the reduced one.
  */
  H_TYPE = pol_x(fetch_user_var("H_TYPE"));
  H_TYPE = gen_0;
  /*
  * Maximal number of knot diagrams that can have their data and the
  * corresponding computation results stored and processed simultaneously.
  * Setting it to a big value doesn't affect memory usage too much, as the
  * actual space allocation for all the data is performed on demand.
  */
  MAX_DIAGRAM_NUM = pol_x(fetch_user_var("MAX_DIAGRAM_NUM"));
  MAX_DIAGRAM_NUM = stoi(100);
  /*
  * Choice of the verbosity level.
  */
  VERBOSE_LEVEL = pol_x(fetch_user_var("VERBOSE_LEVEL"));
  V_SILENT = pol_x(fetch_user_var("V_SILENT"));
  V_WHAT = pol_x(fetch_user_var("V_WHAT"));
  V_PROGRESS = pol_x(fetch_user_var("V_PROGRESS"));
  V_DEBUG = pol_x(fetch_user_var("V_DEBUG"));
  V_SILENT = gen_0;
  /* don't show any messages */
  V_WHAT = gen_1;
  /* report what computations are needed to be performed */
  V_PROGRESS = gen_2;
  /* show progress in computations (if any) */
  V_DEBUG = stoi(5);
  /* print debugging messages (and do some debugging) */
  
  /* print everything but the debugging information by default */
  VERBOSE_LEVEL = gsubgs(V_DEBUG, 1);
  /* perform check that the d^2 is zero (used for debugging) */
  CHECK_D2 = pol_x(fetch_user_var("CHECK_D2"));
  CHECK_D2 = gen_0;
  /*
  * Load other pieces of KhoHo
  */
  gp_read_file(GENtostr_unquoted(KhoHo_data));
  gp_read_file(GENtostr_unquoted(KhoHo_gvars));
  gp_read_file(GENtostr_unquoted(KhoHo_diagr));
  gp_read_file(GENtostr_unquoted(KhoHo_chain));
  gp_read_file(GENtostr_unquoted(KhoHo_odd));
  gp_read_file(GENtostr_unquoted(KhoHo_reduce));
  gp_read_file(GENtostr_unquoted(KhoHo_sign));
  gp_read_file(GENtostr_unquoted(KhoHo_print));
  /*
  * A stupid trick to make t, q, and Q appear before others in the list of
  * variables, and in this order.
  */
  dummy = gaddgs(gsqr(t), 1);
  dummy = gaddgs(gsqr(q), 1);
  dummy = gaddgs(gsqr(Q), 1);
  set_H_type(H_TYPE);
  /*
  * List of torsion factors as returned by rank_SNF (local for each computation)
  */
  H_torsion_list = pol_x(fetch_user_var("H_torsion_list"));
  return;
}

/* ************************************************************************ */

/*
* Given the list of possible torsion orders, compute and assign torsion ranks.
*/
void
T_ranks_assign(GEN D_ID, GEN t_orders)	  /* void */
{
  GEN datapos = gen_0, num_t_orders = gen_0, tors_num = gen_0, tors_mask = gen_0, cur_tors = gen_0, new_entry = gen_0, H_torsion_factors = pol_x(fetch_user_var("H_torsion_factors")), H_torsion_ranks = pol_x(fetch_user_var("H_torsion_ranks")), H_torsion_rank_pols = pol_x(fetch_user_var("H_torsion_rank_pols")), H_torsion_vars = pol_x(fetch_user_var("H_torsion_vars"));
  GEN p1;	  /* vec */
  GEN T2 = pol_x(fetch_user_var("T2")), p2;
  GEN p3;	  /* vec */
  GEN p4, DStore = pol_x(fetch_user_var("DStore"));
  datapos = check_ID(D_ID);
  num_t_orders = stoi(glength(t_orders));
  gel(H_torsion_factors, gtos(datapos)) = emptyCmatrix(D_ID, cgetg(1, t_VEC));
  gel(H_torsion_ranks, gtos(datapos)) = emptyCmatrix(D_ID, cgetg(1, t_COL));
  gel(H_torsion_rank_pols, gtos(datapos)) = emptyCmatrix(D_ID);
  {
    long i;
    p1 = cgetg(gtos(num_t_orders)+1, t_VEC);
    for (i = 1; gcmpsg(i, num_t_orders) <= 0; ++i)
      gel(p1, i) = geval(gconcat(strtoGENstr("T"), gel(t_orders, i)));
  }
  gel(/* ... and the list of variables in the torsion Khovanov polynomial */
  H_torsion_vars, gtos(datapos)) = p1;
  /* nothing to do if there is no torsion */
  if (gequal0(num_t_orders))
    return;
  /* replace T2 with 1 (it must be the first entry) */
  if (gequal(gel(gel(H_torsion_vars, gtos(datapos)), 1), T2))
    gel(gel(H_torsion_vars, gtos(datapos)), 1) = gen_1;
  p2 = gcopy(gel(t_orders, gtos(num_t_orders)));
  {
    long l5;
    p3 = cgetg(gtos(p2)+1, t_VEC);
    for (l5 = 1; gcmpsg(l5, p2) <= 0; ++l5)
      gel(p3, l5) = gen_0;
  }
  /* entry numbers that the torsion of different orders has in t_orders;
  * the last entry corresponds to the maximum torsion order met */
  tors_num = p3;
  {
    GEN i;
    for (i = gen_1; gcmp(i, num_t_orders) <= 0; i = gaddgs(i, 1))
      gel(tors_num, gtos(gel(t_orders, gtos(i)))) = gcopy(i);
  }
  p4 = _.jSize(gel(DStore, gtos(D_ID)));
  {
    GEN j;
    for (j = gen_1; gcmp(j, p4) <= 0; j = gaddgs(j, 1))
    {
      GEN p6;
      p6 = _.iSize(gel(DStore, gtos(D_ID)));
      {
        GEN i;
        for (i = gen_1; gcmp(i, p6) <= 0; i = gaddgs(i, 1))
        {
          GEN p7;	  /* vec */
          long l8;
          GEN p9;	  /* vec */
          {
            long l10;
            p7 = cgetg(gtos(num_t_orders)+1, t_COL);
            for (l10 = 1; gcmpsg(l10, num_t_orders) <= 0; ++l10)
              gel(p7, l10) = gen_0;
          }
          new_entry = p7;
          tors_mask = gtoset(cgetg(1, t_VEC));
          l8 = glength(gcoeff(H_torsion_list, gtos(j), gtos(i)));
          {
            long k;
            for (k = 1; k <= l8; ++k)
            {
              GEN p11;	  /* vec */
              cur_tors = gcopy(gel(tors_num, gtos(gel(gcoeff(H_torsion_list, gtos(j), gtos(i)), k))));
              gel(new_entry, gtos(cur_tors)) = gaddgs(gel(new_entry, gtos(cur_tors)), 1);
              p11 = cgetg(2, t_VEC);
              gel(p11, 1) = gcopy(cur_tors);
              tors_mask = setunion(tors_mask, p11);
            }
          }
          tors_mask = vecsort0(geval(tors_mask), NULL, 0);
          p9 = cgetg(3, t_VEC);
          gel(p9, 1) = gtrans(t_orders);
          gel(p9, 2) = gcopy(new_entry);
          gcoeff(gel(H_torsion_factors, gtos(datapos)), gtos(j), gtos(i)) = extract0(gtrans(gtocol(gtomat(p9))), tors_mask, NULL);
          gcoeff(gel(H_torsion_ranks, gtos(datapos)), gtos(j), gtos(i)) = gcopy(new_entry);
          gcoeff(gel(H_torsion_rank_pols, gtos(datapos)), gtos(j), gtos(i)) = gmul(gel(H_torsion_vars, gtos(datapos)), new_entry);
        }
      }
    }
  }
  return;
}

/*
* Assign ranks and/or torsions of the differential (i,j) after the reduction.
* Return the list of orders of torsion basis elements.
*/
GEN
D_ranks_tors(GEN datapos, GEN i, GEN j, GEN do_rank, GEN do_torsion)
{
  GEN res = gen_0, last0 = gen_0, first1 = gen_0, d_rank = gen_0, d_snf = gen_0, reduced_matr = pol_x(fetch_user_var("reduced_matr")), reduced_D_ranks = pol_x(fetch_user_var("reduced_D_ranks"));
  if (!gequalgs(gcoeff(gel(reduced_matr, gtos(datapos)), gtos(j), gtos(i)), 0))
  {
    long l1;
    res = matsnf0(gcoeff(gel(reduced_matr, gtos(datapos)), gtos(j), gtos(i)), 0);
    last0 = gen_0;
    first1 = stoi(glength(res) + 1);
    l1 = glength(res);
    {
      long i;
      for (i = 1; i <= l1; ++i)
      {
        if (gequal0(gel(res, i)))
        {
          last0 = stoi(i);
          continue;
        }
        if (gequal1(gel(res, i)))
        {
          first1 = stoi(i);
          break;
        }
      }
    }
    d_rank = gsubsg(glength(res), last0);
    if (gcmpgs(gsub(first1, last0), 1) > 0)
      d_snf = extract0(res, Str(mkvec3(gaddgs(last0, 1), strtoGENstr(".."), gsubgs(first1, 1))), NULL);
    else
      d_snf = cgetg(1, t_VEC);
  }
  else
  {
    /* reduced_matr[j, i] is empty */
    d_rank = gen_0;
    d_snf = cgetg(1, t_VEC);
  }
  if (!gequal0(do_torsion))
    gcoeff(H_torsion_list, gtos(j), gtos(gaddgs(i, 1))) = gcopy(d_snf);
  if (!gequal0(do_rank))
    gcoeff(gel(reduced_D_ranks, gtos(datapos)), gtos(j), gtos(i)) = gcopy(d_rank);
  return d_snf;
}

/*
* Compute ranks and/or torsions of differentials for all grades.
*/
GEN
D_inv_factors(GEN D_ID, GEN do_rank, GEN do_torsion)
{
  GEN datapos = gen_0, i_size = gen_0, j_size = gen_0, all_tors = gen_0, I_HRANKS = pol_x(fetch_user_var("I_HRANKS")), I_TORSION = pol_x(fetch_user_var("I_TORSION")), DStore = pol_x(fetch_user_var("DStore")), I_REDUCED = pol_x(fetch_user_var("I_REDUCED")), reduced_D_ranks = pol_x(fetch_user_var("reduced_D_ranks")), chain_D_ranks = pol_x(fetch_user_var("chain_D_ranks")), H_ranks = pol_x(fetch_user_var("H_ranks")), reduced_ranks = pol_x(fetch_user_var("reduced_ranks")), chain_ranks = pol_x(fetch_user_var("chain_ranks")), p1 = gen_0;
  datapos = check_ID(D_ID);
  if (((gequal0(do_rank)) || (!gequal0(do_rank) && gequal(get_info(D_ID, I_HRANKS), strtoGENstr("computed")))) && ((gequal0(do_torsion)) || (!gequal0(do_torsion) && gequal(get_info(D_ID, I_TORSION), strtoGENstr("computed")))))
  {
    pari_printf("  already computed\n");
    return gen_0;
  }
  i_size = _.iSize(gel(DStore, gtos(D_ID)));
  j_size = _.jSize(gel(DStore, gtos(D_ID)));
  if (!gequal(get_info(D_ID, I_REDUCED), strtoGENstr("computed")))
  {
    message(V_WHAT, "Reducing the chain complex first ... ");
    reduce(D_ID);
    message(V_WHAT, "    done with the reduction.");
  }
  if (!gequal0(do_rank))
    gel(reduced_D_ranks, gtos(datapos)) = emptyCmatrix(D_ID);
  if (!gequal0(do_torsion))
    H_torsion_list = emptyCmatrix(D_ID, cgetg(1, t_VEC));
  /* collect orders of all the torsion factors */
  all_tors = cgetg(1, t_VEC);
  {
    GEN j;
    for (j = gen_1; gcmp(j, j_size) <= 0; j = gaddgs(j, 1))
    {
      GEN p2;
      p2 = gsubgs(i_size, 1);
      {
        GEN i;
        for (i = gen_1; gcmp(i, p2) <= 0; i = gaddgs(i, 1))
          all_tors = gconcat(all_tors, D_ranks_tors(datapos, i, j, do_rank, do_torsion));
      }
    }
  }
  if (!gequal0(do_rank))
  {
    gel(chain_D_ranks, gtos(datapos)) = gcopy(gel(reduced_D_ranks, gtos(datapos)));
    gel(H_ranks, gtos(datapos)) = gcopy(gel(reduced_ranks, gtos(datapos)));
    {
      GEN j;
      for (j = gen_1; gcmp(j, j_size) <= 0; j = gaddgs(j, 1))
      {
        gcoeff(gel(/* the first column: zero differential has rank zero */
        H_ranks, gtos(datapos)), gtos(j), 1) = gsub(gcoeff(gel(H_ranks, gtos(datapos)), gtos(j), 1), gcoeff(gel(reduced_D_ranks, gtos(datapos)), gtos(j), 1));
        gcoeff(gel(chain_D_ranks, gtos(datapos)), gtos(j), 1) = gadd(gcoeff(gel(chain_D_ranks, gtos(datapos)), gtos(j), 1), gsub(gcoeff(gel(chain_ranks, gtos(datapos)), gtos(j), 1), gcoeff(gel(reduced_ranks, gtos(datapos)), gtos(j), 1)));
        /* general case */
        {
          GEN i;
          for (i = gen_2; gcmp(i, i_size) <= 0; i = gaddgs(i, 1))
          {
            gcoeff(gel(H_ranks, gtos(datapos)), gtos(j), gtos(i)) = gsub(gcoeff(gel(H_ranks, gtos(datapos)), gtos(j), gtos(i)), gadd(gcoeff(gel(reduced_D_ranks, gtos(datapos)), gtos(j), gtos(i)), gcoeff(gel(reduced_D_ranks, gtos(datapos)), gtos(j), gtos(gsubgs(i, 1)))));
            gcoeff(gel(/* deduce ranks of the original chain complex
            * from the reduced ones */
            chain_D_ranks, gtos(datapos)), gtos(j), gtos(i)) = gadd(gcoeff(gel(chain_D_ranks, gtos(datapos)), gtos(j), gtos(i)), gadd(gsub(gsub(gcoeff(gel(chain_ranks, gtos(datapos)), gtos(j), gtos(i)), gcoeff(gel(reduced_ranks, gtos(datapos)), gtos(j), gtos(i))), gcoeff(gel(chain_D_ranks, gtos(datapos)), gtos(j), gtos(gsubgs(i, 1)))), gcoeff(gel(reduced_D_ranks, gtos(datapos)), gtos(j), gtos(gsubgs(i, 1)))));
          }
        }
        /* rank of the last (zero) differential must be zero */
        if (!gequalgs(gcoeff(gel(chain_D_ranks, gtos(datapos)), gtos(j), gtos(i_size)), 0))
          pari_err(e_MISC, "D_inv_factors: wrong complex ranks");
      }
    }
    set_info(D_ID, I_HRANKS, "computed");
  }
  if (!gequal0(do_torsion))
  {
    /* make the list of unique torsion orders
    * and assign torsion ranks */
    T_ranks_assign(D_ID, vecsort0(geval(gtoset(all_tors)), NULL, 0));
    p1 = set_info(D_ID, I_TORSION, "computed");
  }
  return p1;
}

/*
* Compute Betti numbers for all grades.
*/
GEN
Betti(GEN D_ID)
{
  return D_inv_factors(D_ID, gen_1, gen_0);
}

/*
* Compute torsion for all grades.
*/
GEN
Torsion(GEN D_ID)
{
  return D_inv_factors(D_ID, gen_0, gen_1);
}

/* ************************************************************************ */

GEN
matrix2pol(GEN D_ID, GEN matr_name, GEN ret_vector, long prec)
{
  GEN datapos = gen_0, q_vec = gen_0, t_vec = gen_0, res = gen_0, p1, DStore = pol_x(fetch_user_var("DStore"));
  GEN p2;	  /* vec */
  GEN p3;
  GEN p4;	  /* vec */
  if (!ret_vector)
    ret_vector = gen_0;
  datapos = check_ID(D_ID);
  p1 = _.jSize(gel(DStore, gtos(D_ID)));
  {
    long j;
    p2 = cgetg(gtos(p1)+1, t_VEC);
    for (j = 1; gcmpsg(j, p1) <= 0; ++j)
      gel(p2, j) = gpow(q, m2j(D_ID, j), prec);
  }
  q_vec = p2;
  p3 = _.iSize(gel(DStore, gtos(D_ID)));
  {
    long i;
    p4 = cgetg(gtos(p3)+1, t_COL);
    for (i = 1; gcmpsg(i, p3) <= 0; ++i)
      gel(p4, i) = gpow(t, m2i(D_ID, i), prec);
  }
  t_vec = p4;
  if (!gequal0(ret_vector))
  {
    t_vec = diagonal(t_vec);
    q_vec = diagonal(q_vec);
  }
  res = gmul(gmul(q_vec, extract0(geval(gconcat(matr_name, strtoGENstr("[datapos]"))), strtoGENstr("^-1"), strtoGENstr("^-1"))), t_vec);
  if (!gequal0(ret_vector))
    res = gtrans(gconcat1(gtovec(res)));
  return res;
}

/*
* Compute the rational Khovanov polynomial (over \Bbb Q) in t and q variables.
* If ret_vector is not 0, return the vector of all monomials instead.
*/
GEN
KhPol_Q(GEN D_ID, GEN ret_vector, long prec)
{
  GEN I_HRANKS = pol_x(fetch_user_var("I_HRANKS"));
  if (!ret_vector)
    ret_vector = gen_0;
  if (!gequal(get_info(D_ID, I_HRANKS), strtoGENstr("computed")))
  {
    message(V_WHAT, "Computing Betti numbers ... ");
    Betti(D_ID);
    message(V_WHAT, "   ... done with computing Betti numbers.");
  }
  return matrix2pol(D_ID, strtoGENstr("H_ranks"), ret_vector, prec);
}

/*
* Compute the torsion Khovanov polynomial in t, Q, and T{i} variables,
* where T{i} corresponds to i-th torsion for i > 2. 
* If ret_vector is not 0, return the vector of all monomials in t and Q.
*/
GEN
KhPol_T(GEN D_ID, GEN ret_vector, long prec)
{
  GEN I_TORSION = pol_x(fetch_user_var("I_TORSION"));
  if (!ret_vector)
    ret_vector = gen_0;
  if (!gequal(get_info(D_ID, I_TORSION), strtoGENstr("computed")))
  {
    message(V_WHAT, "Computing homology torsion ... ");
    Torsion(D_ID);
    message(V_WHAT, "   ... done with computing the torsion.");
  }
  return gsubst(matrix2pol(D_ID, strtoGENstr("H_torsion_rank_pols"), ret_vector, prec), gvar(q), Q);
}

/*
* Compute the ``expanded'' Khovanov polynomial in t, q, and Q{i} variables.
* If split is not 0, return the rational and torsion polynomials separately.
* If ret_vector is not 0, use vectors of monomials instead of polynomials.
*/
GEN
KhPol(GEN D_ID, GEN split, GEN ret_vector, long prec)
{
  GEN I_TORSION = pol_x(fetch_user_var("I_TORSION")), I_HRANKS = pol_x(fetch_user_var("I_HRANKS")), p1 = gen_0;
  if (!split)
    split = gen_1;
  if (!ret_vector)
    ret_vector = gen_0;
  check_ID(D_ID);
  if (!gequal(get_info(D_ID, I_TORSION), strtoGENstr("computed")) && !gequal(get_info(D_ID, I_HRANKS), strtoGENstr("computed")))
  {
    message(V_WHAT, "Computing differential ranks and torsions ... ");
    D_inv_factors(D_ID, gen_1, gen_1);
    message(V_WHAT, "   ... done with computing ranks and torsions.");
  }
  if (!gequal0(split))
  {
    GEN p2;	  /* vec */
    p2 = cgetg(3, t_VEC);
    gel(p2, 1) = KhPol_Q(D_ID, ret_vector, prec);
    gel(p2, 2) = KhPol_T(D_ID, ret_vector, prec);
    p1 = p2;
  }
  else
    p1 = gadd(KhPol_Q(D_ID, gen_0, prec), KhPol_T(D_ID, gen_0, prec));
  return p1;
}

/* ************************************************************************ */

/*
* Given a linking matrix, compute the main factor in the extended
* Bar-Natan's Conjecture 1 for links. The conjecture is proved by
* Eun Soo Lee for nonsplit alternating links. The factor is
* \sum_{E\subset\{2,\cdots,n\}} (tq^2)^{\sum_{j\in E,k\not\in E} 2lk_{jk}}
*/
GEN
conj1_factor(GEN lmatr, long prec)
{
  GEN msize = gen_0, res = gen_0, p1, x = pol_x(fetch_user_var("x"));
  /* lmatr should be square */
  msize = gcopy(gel(matsize(lmatr), 1));
  res = gen_0;
  /* E goes over all subsets of {2, ..., n} */
  p1 = gsubgs(gpow(gen_2, msize, prec), 1);
  {
    GEN E;
    for (E = gen_0; gcmp(E, p1) <= 0; E = gaddgs(E, 2))
      /* leave in lmatr only columns from E and rows not from E
      * and sum all the entries there up */
      res = gadd(res, gpow(x, sum_entries(extract0(lmatr, E, gsub(gsubgs(gpow(gen_2, msize, prec), 1), E))), prec));
  }
  return gsubst(res, gvar(x), gsqr(gmul(t, gsqr(q))));
}

/*
* Check whether the Khovanov polynomial of a link satisfies the extended
* Conjecture 1. cfactor should be computed as above (is always 1 for knots).
* Return a vector consisting of the exponent s of q and Kh',
* or an error message if the conjecture fails.
*/
GEN
check_conj1(GEN khpolQ, GEN cfactor, long prec)
{
  GEN s = gen_0, Kh_p = gen_0, qs_1 = gen_0, guy1 = gen_0, guy2 = gen_0, khden = gen_0, den_qdeg = gen_0, den_tdeg = gen_0, extra = pol_x(fetch_user_var("extra"));
  GEN p1;	  /* vec */
  if (!cfactor)
    cfactor = gen_1;
  guy1 = gaddsg(1, gmul(t, gpowgs(q, 4)));
  guy2 = gaddsg(1, gsqr(q));
  /* take everything modulo 1+tq^4 to figure out the term in front */
  qs_1 = simplify(gdiv(gmod(khpolQ, guy1), gmul(guy2, gmod(cfactor, guy1))));
  /* check that qs_1 is indeed q^{s-1} */
  s = gaddgs(gppoldegree(qs_1, gvar(q)), 1);
  if (!gequalgs(gsub(qs_1, gpow(q, gsubgs(s, 1), prec)), 0))
  {
    pari_printf("\n\nWrong conjecture 1: no s\n\n");
    return strtoGENstr("Conjecture 1 failed: no s");
  }
  /* check that Kh_p is a polynomial (in PARI presentation) */
  /* dirty trick to avoid a bug in Pari */
  extra = gmul(gpowgs(q, 100), gpowgs(t, 100));
  Kh_p = gdiv(gdiv(gmul(gsub(gmul(khpolQ, gpow(q, gsubsg(1, s), prec)), gmul(guy2, cfactor)), extra), guy1), extra);
  khden = denom(Kh_p);
  den_qdeg = gppoldegree(khden, gvar(q));
  den_tdeg = gppoldegree(khden, gvar(t));
  if (typ(simplify(gdiv(khden, gmul(gpow(q, den_qdeg, prec), gpow(t, den_tdeg, prec))))) != t_INT)
  {
    pari_printf("\n\nWrong conjecture 1: Kh_p is not a polynomial\n\n");
    return strtoGENstr("Conjecture 1 failed: Kh_p is not a polynomial");
  }
  p1 = cgetg(3, t_VEC);
  gel(p1, 1) = gcopy(s);
  gel(p1, 2) = gcopy(Kh_p);
  return p1;
}

/*
* Check whether the Khovanov polynomial is H-thick (cfactor is optional).
* Return Kh' if yes, 0 if no, or an error message if the conjecture 1 fails.
*/
GEN
check_H_thick(GEN khpolQ, GEN cfactor, long prec)
{
  GEN conj_answer = gen_0, Kh_p = gen_0;
  if (!cfactor)
    cfactor = gen_1;
  conj_answer = check_conj1(khpolQ, cfactor, prec);
  if (typ(conj_answer) != t_VEC)
    return conj_answer;
  Kh_p = gcopy(gel(conj_answer, 2));
  if (!gequalgs(gsub(Kh_p, gsubst(gsubst(Kh_p, gvar(q), gen_1), gvar(t), gmul(t, gsqr(q)))), 0))
    return Kh_p;
  else
    return gen_0;
  return gen_0;
}

/*
* Check whether the Khovanov polynomial is T-thick (cfactor is optional).
* Return the defect of torsion, or an error message if the conjecture 1 fails.
*/
GEN
check_T_thick(GEN khpolQ, GEN khpolT, GEN cfactor, long prec)
{
  GEN conj_answer = gen_0, s = gen_0, Kh_p = gen_0;
  if (!cfactor)
    cfactor = gen_1;
  conj_answer = check_conj1(khpolQ, cfactor, prec);
  if (typ(conj_answer) != t_VEC)
    return conj_answer;
  s = gcopy(gel(conj_answer, 1));
  Kh_p = gcopy(gel(conj_answer, 2));
  return gsub(gmul(gmul(Kh_p, t), gpow(q, gaddgs(s, 1), prec)), gsubst(khpolT, gvar(Q), q));
}

/*
* Given a Khovanov polynomial, figure out how many diagonals it occupies.
* If min_max is not zero, return degrees of the lowest and highest diagonals.
*/
GEN
pol_diags(GEN khpolQ, GEN min_max, long prec)
{
  GEN khflat = gen_0, maxdeg = gen_0, mindeg = gen_0, p1 = gen_0;
  if (!min_max)
    min_max = gen_0;
  /* make diagonals into rows */
  khflat = gsubst(khpolQ, gvar(t), gdiv(t, gsqr(q)));
  /* workaround for a PARI's bug in poldegree */
  khflat = gsubst(khflat, gvar(t), mppi(prec));
  /* maximal and minimal degree in q */
  maxdeg = gppoldegree(khflat, gvar(q));
  mindeg = gneg(gppoldegree(gsubst(khflat, gvar(q), ginv(q)), gvar(q)));
  if (!gequal0(min_max))
  {
    GEN p2;	  /* vec */
    p2 = cgetg(3, t_VEC);
    gel(p2, 1) = gcopy(mindeg);
    gel(p2, 2) = gcopy(maxdeg);
    p1 = p2;
  }
  else
    p1 = gaddgs(gdivgs(gsub(maxdeg, mindeg), 2), 1);
  return p1;
}

