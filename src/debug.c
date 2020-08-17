/*
    Copyright (C) 2016-2019 Tomas Flouri and Ziheng Yang

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <t.flouris@ucl.ac.uk>,
    Department of Genetics, Evolution and Environment,
    University College London, Gower Street, London WC1E 6BT, England
*/

#include "bpp.h"

#define SHRINK 1
#define EXPAND 2

void debug_snl_stage1(stree_t * stree,
                      gtree_t ** gtree_list,
                      snode_t * y,
                      snode_t * c,
                      snode_t * a,
                      long movetype,
                      int downwards,
                      double tauy,
                      double newtauy)
{
  long i;
  double taufactor;
  char * gtree_newick = NULL;
  char * stree_newick = NULL;
  char * (*cb_stree_branch)(const snode_t *) = NULL;
  char * (*cb_gtree_branch)(const gnode_t *) = NULL;
  char * stree_desc = NULL;
  char * gtree_desc = NULL;

  if (opt_debug_counter < opt_debug_start || opt_debug_counter > opt_debug_end)
    return;

  if (opt_debug_snl > 1)
  {
    fprintf(stdout,
            "\n[DBG-%ld] ------ Proposal: GSPR (SNL) Number: %ld ------\n",
            opt_debug_snl,
            opt_debug_counter);

    switch (opt_debug_snl)
    {
      case 2:
        cb_stree_branch = cb_stree_print_none;
        cb_gtree_branch = cb_gtree_print_none;
        stree_desc = (char *)xmalloc(sizeof(char)); stree_desc[0] = '\0';
        gtree_desc = (char *)xmalloc(sizeof(char)); gtree_desc[0] = '\0';
        break;
      case 3:
        cb_stree_branch = cb_stree_print_tau;
        cb_gtree_branch = cb_gtree_print_age;
        stree_desc = xstrdup(" (with taus)");
        gtree_desc = xstrdup(" (with ages)");
        break;
      case 4:
      default:
        cb_stree_branch = cb_stree_print_node_tau;
        cb_gtree_branch = cb_gtree_print_node_age;
        stree_desc = xstrdup(" (with node indices and taus)");
        gtree_desc = xstrdup(" (with node indices and ages)");
        break;
    }

    stree_newick = stree_export_newick(stree->root,cb_stree_branch);
    fprintf(stdout,
            "[DBG-%ld] Species tree before move%s: %s\n",
            opt_debug_snl, stree_desc, stree_newick);
    free(stree_newick);
    free(stree_desc);

    taufactor = newtauy / tauy;
    fprintf(stdout, "[DBG-%ld] Move: %s  |  nodes Y: %s C: %s A: %s  |  "
            "tau_Y: %f tau*_Y: %f tau_Y/tau*_Y: %f\n",
            opt_debug_counter,
            movetype == EXPAND && !downwards ? 
              "EXPAND" : movetype == SHRINK ? "SHRINK" : "EXPAND+SHRINK",
            y->label,
            c->label,
            a->label,
            tauy,
            newtauy,
            taufactor);

    fprintf(stdout,
            "[DBG-%ld] Gene trees before proposal%s:\n",
            opt_debug_snl, gtree_desc);
    free(gtree_desc);
    for (i = 0; i < opt_locus_count; ++i)
    {
      gtree_newick = gtree_export_newick(gtree_list[i]->root,
                                               cb_gtree_branch);
      printf("[DBG-%ld]  %ld  ->  %s\n", opt_debug_snl, i+1, gtree_newick);
      free(gtree_newick);
    }
  }
}

void debug_snl_stage2(stree_t * stree,
                      gtree_t ** gtree,
                      double logpr_notheta,
                      double r,
                      double lnacceptance)
{
  long i;
  char * gtree_newick = NULL;
  char * stree_newick = NULL;
  char * (*cb_stree_branch)(const snode_t *) = NULL;
  char * (*cb_gtree_branch)(const gnode_t *) = NULL;
  char * stree_desc = NULL;
  char * gtree_desc = NULL;

  if (opt_debug_counter < opt_debug_start || opt_debug_counter > opt_debug_end)
    return;

  if (opt_debug_snl == 1)
  {
    long rc = (lnacceptance >= -1e-10 || r < exp(lnacceptance));
    fprintf(stdout,
            "[DBG-%ld] Proposal: GSPR (SNL) Number: %ld  lnacceptance: %f  |  %s\n",
            opt_debug_snl, opt_debug_counter, lnacceptance,
            rc ? "ACCEPTED" : "REJECTED");
    return;
  }

  if (opt_debug_snl > 1)
  {
    switch (opt_debug_snl)
    {
      case 2:
        cb_stree_branch = cb_stree_print_none;
        cb_gtree_branch = cb_gtree_print_none;
        stree_desc = (char *)xmalloc(sizeof(char)); stree_desc[0] = '\0';
        gtree_desc = (char *)xmalloc(sizeof(char)); gtree_desc[0] = '\0';
        break;
      case 3:
        cb_stree_branch = cb_stree_print_tau;
        cb_gtree_branch = cb_gtree_print_age;
        stree_desc = xstrdup(" (with taus)");
        gtree_desc = xstrdup(" (with ages)");
        break;
      case 4:
      default:
        cb_stree_branch = cb_stree_print_node_tau;
        cb_gtree_branch = cb_gtree_print_node_age;
        stree_desc = xstrdup(" (with node indices and taus)");
        gtree_desc = xstrdup(" (with node indices and ages)");
        break;
    }

    stree_newick = stree_export_newick(stree->root,cb_stree_branch);
    fprintf(stdout,
            "[DBG-%ld] Proposed species tree%s: %s\n",
            opt_debug_snl, stree_desc, stree_newick);
    free(stree_newick);
    free(stree_desc);

    fprintf(stdout,
            "[DBG-%ld] Proposed gene trees%s:\n",
            opt_debug_snl, gtree_desc);
    free(gtree_desc);
    for (i = 0; i < opt_locus_count; ++i)
    {
      gtree_newick = gtree_export_newick(gtree[i]->root,cb_gtree_branch);
      fprintf(stdout,
              "[DBG-%ld]  %ld  ->  %s\n",
              opt_debug_snl, i+1, gtree_newick);
      free(gtree_newick);
    }

    fprintf(stdout,
            "[DBG-%ld] log-L and gene tree probabilities:\n",
            opt_debug_snl);
    for (i = 0; i < opt_locus_count; ++i)
    {
      if (opt_est_theta)
        fprintf(stdout,
                "[DBG-%ld]   %ld  ->  log-L: %f  old log-L: %f  log-L diff: %f   "
                "|   log-gP: %f  old log-gP: %f  log-gP diff: %f\n",
                opt_debug_snl, i+1,
                gtree[i]->logl, gtree[i]->old_logl,
                gtree[i]->logl - gtree[i]->old_logl,
                gtree[i]->logpr, gtree[i]->old_logpr,
                gtree[i]->logpr - gtree[i]->old_logpr);
      else
        fprintf(stdout,
                "[DBG-%ld]   %ld  ->  log-L: %f  old log-L: %f  log-L diff: %f\n",
                opt_debug_snl, i+1,
                gtree[i]->logl, gtree[i]->old_logl,
                gtree[i]->logl - gtree[i]->old_logl);
    }
    if (!opt_est_theta)
    {
      fprintf(stdout,
              "[DBG-%ld] log-gP: %f  old log-gP: %f  log-gP diff: %f\n",
              opt_debug_snl, logpr_notheta, stree->notheta_logpr,
              logpr_notheta - stree->notheta_logpr);
    }

    fprintf(stdout,
            "[DBG-%ld] GSPR (SNL) log-acceptance ratio (lnacceptance): %f\n",
            opt_debug_snl, lnacceptance);
    long rc = (lnacceptance >= -1e-10 || r < exp(lnacceptance));
    if (rc)
      fprintf(stdout, "[DBG-%ld] Proposal accepted", opt_debug_snl);
    else
      fprintf(stdout, "[DBG-%ld] Proposal rejected", opt_debug_snl);
    fprintf(stdout,
           "\n[DBG-%ld] ---------------------------------------------\n",
           opt_debug_snl);
  }
}
