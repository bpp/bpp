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

void debug_print_gtree(gtree_t * gtree)
{
  long i;

  printf("Label        Node-index  Child1-index  Child2-index  Parent-index  "
         "Pmat-index    Pop     time\n");
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
  {
    char * label;
    if (gtree->nodes[i]->label)
      label = xstrdup(gtree->nodes[i]->label);
    else
      label = xstrdup("N/A");

    /* shorten label if over 12 characters */
    if (strlen(label) > 12)
    {
      label[9] = '.'; label[10] = '.'; label[11] = '.'; label[12] = 0;
    }

    printf("%-12s", label);
    free(label);
    printf("  %9d", gtree->nodes[i]->node_index);

    if (gtree->nodes[i]->left)
      printf("  %12d", gtree->nodes[i]->left->node_index);
    else
      printf("  %12s", "N/A");

    if (gtree->nodes[i]->right)
      printf("  %12d", gtree->nodes[i]->right->node_index);
    else
      printf("  %12s", "N/A");


    if (gtree->nodes[i]->parent)
      printf("  %12d", gtree->nodes[i]->parent->node_index);
    else
      printf("  %12s", "N/A");

    printf("  %10d", gtree->nodes[i]->pmatrix_index);

    printf(" %5s", gtree->nodes[i]->pop->label);
    printf(" %.12f", gtree->nodes[i]->time);
    printf("\n");
  }
}

void debug_check_relations(stree_t * stree, gtree_t * gtree, long msa_index)
{
  long i,j;

  assert(gtree->root && !gtree->root->parent);

  /* marks */
  for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
  {
    if (!gtree->nodes[i]->parent)
      assert(gtree->nodes[i] == gtree->root);
    assert(!gtree->nodes[i]->mark);
  }

  for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
  {
    gnode_t * x = gtree->nodes[i];

    while (x && !x->mark)
    {
      x->mark = 1;
      
      if (x->parent)
      {
        assert(x->parent->left && x->parent->right);
        assert(x->parent->left == x || x->parent->right == x);
      }
      x = x->parent;
    }
  }

  /* marks */
  for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
  {
    assert(gtree->nodes[i]->mark);
    gtree->nodes[i]->mark = 0;
  }

  /* check leaves */
  for (i = 0; i < gtree->tip_count; ++i)
    assert(gtree->nodes[i]->leaves == 1);
  for (i = gtree->tip_count; i < gtree->tip_count+gtree->inner_count; ++i)
  {
    gnode_t * x = gtree->nodes[i];
    assert(x->leaves == x->left->leaves + x->right->leaves);
  }

  /* check seqin count, coalescences etc */
  long sum = 0;
  for (i = 0; i < stree->tip_count; ++i)
  {
    sum += stree->nodes[i]->seqin_count[msa_index];
  }
  assert(sum == gtree->tip_count);
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode_t * snode = stree->nodes[i];

    long lout = snode->left->seqin_count[msa_index] -
                snode->left->event_count[msa_index];
    long rout = snode->right->seqin_count[msa_index] -
                snode->right->event_count[msa_index];

    assert(snode->seqin_count[msa_index] == lout+rout);
  }
  snode_t * sroot = stree->root;
  assert(sroot->seqin_count[msa_index] - sroot->event_count[msa_index] == 1);

  /* check ages */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
  {
    gnode_t * x = gtree->nodes[i];

    if (!x->parent)
    {
      assert(x->left->time < x->time && x->right->time < x->time);
    }
    else
    {
      assert(x->parent->time > x->time);

      if (x->left)
      {
        assert(x->left->time < x->time && x->right->time < x->time);
      }
    }

    if (x->pop->parent)
      assert(x->time < x->pop->parent->tau);
    if (x->left)
      assert(x->time > x->pop->tau);
  }

  /* check for duplicate pmatrix index */
  /* TODO: Also check for the corresponding index, i.e. pmatindex+edge_count */
  for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
  {
    for (j = i+1; j < gtree->tip_count+gtree->inner_count; ++j)
    {
      if (gtree->nodes[i]->parent && gtree->nodes[j]->parent)
        assert(gtree->nodes[i]->pmatrix_index != gtree->nodes[j]->pmatrix_index);
    }
  }
}

void debug_check_leaves(gtree_t ** gtree)
{
  long i,j;

  for (i = 0; i < opt_locus_count; ++i)
  {
    for (j = 0; j < gtree[i]->tip_count+gtree[i]->inner_count; ++j)
    {
      gnode_t * node = gtree[i]->nodes[j];

      if (node->left)
      {
        /* inner */
        assert(node->leaves == node->left->leaves + node->right->leaves);
      }
      else
      {
        /* tip */
        assert(node->leaves == 1);
      }
    }
  }
}

static void debug_validate_logpg_notheta(stree_t * stree,
                                         gtree_t ** gtree,
                                         locus_t ** locus,
                                         const char * move)
{
  size_t hnodes_total;
  long i,j;
  double logpr, old_logpr;
  const long thread_index = 0;

  hnodes_total = (size_t)(2*stree->hybrid_count);
  old_logpr = stree->notheta_logpr;
  logpr = 0;

  /* temporary storage */
  double * phi_contrib = (double *)xmalloc(hnodes_total *
                                           (size_t)opt_locus_count *
                                           sizeof(double));
  double * old_phi_contrib = (double *)xmalloc(hnodes_total *
                                               (size_t)opt_locus_count *
                                               sizeof(double));
  double * phi_sum = (double *)xmalloc(hnodes_total * sizeof(double));

  double * phi_tmp = phi_contrib;
  double * old_phi_tmp = old_phi_contrib;
  double * phi_sum_tmp = phi_sum;
  for (i = 0; i < stree->hybrid_count; ++i)
  {
    snode_t * x = stree->nodes[stree->tip_count+stree->inner_count+i];

    /* copy old values from mirror hybrid node */
    memcpy(phi_tmp, x->notheta_phi_contrib, opt_locus_count * sizeof(double));
    memcpy(old_phi_tmp,
           x->notheta_old_phi_contrib,
           opt_locus_count * sizeof(double));
    *phi_sum_tmp = x->hphi_sum;

    /* zero-out the arrays */
    memset(x->notheta_phi_contrib,0,opt_locus_count*sizeof(double));
    memset(x->notheta_old_phi_contrib,0,opt_locus_count*sizeof(double));
    x->hphi_sum = 0;

    /* move storage pointers */
    phi_tmp += opt_locus_count;
    old_phi_tmp += opt_locus_count;
    phi_sum_tmp++;

    x = stree->nodes[stree->tip_count+stree->inner_count+i]->hybrid;

    /* copy old values from hybrid node */
    memcpy(phi_tmp, x->notheta_phi_contrib, opt_locus_count * sizeof(double));
    memcpy(old_phi_tmp,
           x->notheta_old_phi_contrib,
           opt_locus_count * sizeof(double));
    *phi_sum_tmp = x->hphi_sum;

    /* zero-out the arrays */
    memset(x->notheta_phi_contrib,0,opt_locus_count*sizeof(double));
    memset(x->notheta_old_phi_contrib,0,opt_locus_count*sizeof(double));
    x->hphi_sum = 0;

    /* move storage pointers */
    phi_tmp += opt_locus_count;
    old_phi_tmp += opt_locus_count;
    phi_sum_tmp++;
  }


  /* recompute population contributions */
  for (i = 0; i < opt_locus_count; ++i)
    logpr = gtree_logprob(stree,locus[i]->heredity[0],i,thread_index);

  /* add constant part */
  logpr += stree->notheta_hfactor+stree->notheta_sfactor;

  phi_sum_tmp = phi_sum;
  phi_tmp = phi_contrib;
  for (i = 0; i < stree->hybrid_count; ++i)
  {
    /* mirror hybrid node */
    snode_t * x = stree->nodes[stree->tip_count+stree->inner_count+i];

    /* validate phi contributions */
    for (j = 0; j < opt_locus_count; ++j)
      if (fabs(x->notheta_phi_contrib[j] - phi_tmp[j]) > 1e-9)
        fatal("[FATAL-%ld] move: %s hybrid: %d mirror: YES    "
              "phi_contrib[%ld]: %f  old_contrib[%ld]: %f    "
              "seqin_count: %d phi: %f\n",
              opt_debug_counter,
              move,
              i,
              j,
              x->notheta_phi_contrib[j],
              j,
              phi_tmp[j],
              x->seqin_count[j],
              x->hphi);


    if (fabs(x->hphi_sum - *phi_sum_tmp) > 1e-9)
    {
      fatal("[FATAL-%ld] move: %s hybrid: %d mirror: YES    "
            "hphi_sum (correct): %f  old hphi_sum (wrong): %f\n",
            opt_debug_counter,
            move,
            x->hphi_sum,
            *phi_sum_tmp);
    }
    ++phi_sum_tmp;
    phi_tmp += opt_locus_count;

    /* normal hybrid node */
    x = stree->nodes[stree->tip_count+stree->inner_count+i]->hybrid;

    /* validate phi contributions */
    for (j = 0; j < opt_locus_count; ++j)
      if (fabs(x->notheta_phi_contrib[j] - phi_tmp[j]) > 1e-9)
        fatal("[FATAL-%ld] move: %s hybrid: %d mirror: NO    "
              "phi_contrib[%ld]: %f  old_contrib[%ld]: %f    "
              "seqin_count: %d phi: %f\n",
              opt_debug_counter,
              move,
              i,
              j,
              x->notheta_phi_contrib[j],
              j,
              phi_tmp[j],
              x->seqin_count[j],
              x->hphi);

    if (fabs(x->hphi_sum - *phi_sum_tmp) > 1e-9)
    {
      fatal("[FATAL-%ld] move: %s hybrid: %d mirror: YES    "
            "hphi_sum (correct): %f  old hphi_sum (wrong): %f\n",
            opt_debug_counter,
            move,
            x->hphi_sum,
            *phi_sum_tmp);
    }
    ++phi_sum_tmp;
    phi_tmp += opt_locus_count;
  }

  if (fabs(logpr - old_logpr) > 1e-5)
    fatal("[FATAL-%ld] move: %s logpr (correct): %f  old logpr (wrong): %f\n",
          opt_debug_counter, move, logpr, old_logpr);

  /* restore */
  phi_tmp = phi_contrib;
  old_phi_tmp = old_phi_contrib;
  phi_sum_tmp = phi_sum;
  for (i = 0; i < stree->hybrid_count; ++i)
  {
    /* normal hybrid node */
    snode_t * x = stree->nodes[stree->tip_count+stree->inner_count+i];

    /* restore elements */
    memcpy(x->notheta_phi_contrib, phi_tmp, opt_locus_count * sizeof(double));
    memcpy(x->notheta_old_phi_contrib,
           old_phi_tmp,
           opt_locus_count * sizeof(double));
    x->hphi_sum = *phi_sum_tmp;

    /* move pointers */
    phi_tmp += opt_locus_count;
    old_phi_tmp += opt_locus_count;
    phi_sum_tmp++;

    /* mirror hybrid node */
    x = stree->nodes[stree->tip_count+stree->inner_count+i]->hybrid;

    /* restore elements */
    memcpy(x->notheta_phi_contrib, phi_tmp, opt_locus_count * sizeof(double));
    memcpy(x->notheta_old_phi_contrib,
           old_phi_tmp,
           opt_locus_count * sizeof(double));
    x->hphi_sum = *phi_sum_tmp;

    /* move pointers */
    phi_tmp += opt_locus_count;
    old_phi_tmp += opt_locus_count;
    phi_sum_tmp++;
  }

  free(phi_contrib);
  free(old_phi_contrib);
  free(phi_sum);
}

void debug_validate_logpg(stree_t * stree,
                          gtree_t ** gtree,
                          locus_t ** locus,
                          const char * move)
{
  if (!opt_est_theta)
    debug_validate_logpg_notheta(stree,gtree,locus,move);
  else
    fatal("Testing of logPG with theta not implemented...");
}

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
            "\n[DBG-%ld] ------ Proposal: SNL Number: %ld ------\n",
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

void debug_gage_stage2(stree_t * stree,
                       gtree_t ** gtree,
                       double logpr_notheta,
                       double r,
                       double lnacceptance)
{
  if (opt_debug_gage == 1)
  {
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
            "[DBG-%ld] Proposal: SNL Number: %ld  lnacceptance: %f  |  %s\n",
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
            "[DBG-%ld] SNL log-acceptance ratio (lnacceptance): %f\n",
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

/* callback function to print the population label for each inner gtree node */
char * cb_gtree_print_pop(const gnode_t * node)
{
  char * s = NULL;
  if (node->left)      /* has a child <=> inner node */
  {
    if (node->pop->label)
      xasprintf(&s, "[%s]", node->pop->label);
    else
      xasprintf(&s, "[%d]", node->pop->node_index);
  }
  else
  {
    /* tip node */
    xasprintf(&s, "%s", node->label);
  }

  return s;
}

/* alternative callback function to print the population label AND age for each
   inner gtree node */
char * cb_gtree_print_pop_and_age(const gnode_t * node)
{
  char * s = NULL;
  if (node->left)      /* has a child <=> inner node */
  {
    if (node->pop->label)
      xasprintf(&s, "[%s]:%f", node->pop->label, node->time);
    else
      xasprintf(&s, "[%d]:%f", node->pop->node_index, node->time);
  }
  else
  {
    /* tip node */
    xasprintf(&s, "%s", node->label);
  }

  return s;
}

void debug_bruce(stree_t * stree,
                 gtree_t ** gtree,
                 const char * move,
                 long iter,
                 FILE * fp_out)
{
  long i;
  char * newick;

  fprintf(fp_out, "%ld - %s - Gene trees:\n", iter, move);
  for (i = 0; i < opt_locus_count; ++i)
  {
    newick = gtree_export_newick(gtree[i]->root,cb_gtree_print_pop);
    fprintf(fp_out, "  %s\n", newick);
    free(newick);
  }
}
