/*
    Copyright (C) 2016-2022 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

void debug_print_migrations(stree_t * stree)
{
  long i;
  for (i = 0; i < opt_migration_count; ++i)
  {
    migspec_t * spec = opt_mig_specs+i;
    unsigned int si = spec->si;
    unsigned int ti = spec->ti;
    printf("  M_%s->%s=%f (%d,%d,%f,%f)",
           stree->nodes[si]->label,
           stree->nodes[ti]->label,
           spec->M,
           spec->si,
           spec->ti,
           spec->alpha,
           spec->beta);
  }
  printf("\n");
}

void debug_print_migmatrix(stree_t * stree)
{
  long i,j;
  long total_nodes = stree->tip_count + stree->inner_count;

  printf("      ");
  for (i = 0; i < total_nodes; ++i)
    printf(" %2s", stree->nodes[i]->label);
  printf("\n");
  for (i = 0; i < total_nodes; ++i)
  {
    printf("   ");
    printf(" %2s", stree->nodes[i]->label);
    for (j = 0; j < total_nodes; ++j)
    {
      if (opt_migration_matrix[i][j] != -1)
        printf(ANSI_COLOR_RED " %2ld" ANSI_COLOR_RESET, opt_migration_matrix[i][j]);
      else
        printf(" %2ld", opt_migration_matrix[i][j]);
    }
    printf("\n");
  }
}


void debug_print_bitmatrix(stree_t * stree)
{
  long i,j;
  long total_nodes = stree->tip_count + stree->inner_count;

  printf("      ");
  for (i = 0; i < total_nodes; ++i)
    printf(" %2s", stree->nodes[i]->label);
  printf("\n");
  for (i = 0; i < total_nodes; ++i)
  {
    printf("   ");
    printf(" %2s", stree->nodes[i]->label);
    for (j = 0; j < total_nodes; ++j)
    {
      if (opt_mig_bitmatrix[i][j] != 0)
        printf(ANSI_COLOR_RED " %2ld" ANSI_COLOR_RESET, opt_mig_bitmatrix[i][j]);
      else
        printf(" %2ld", opt_mig_bitmatrix[i][j]);
    }
    printf("\n");
  }
}

void debug_im_check_sum(stree_t * stree, gtree_t ** gtree)
{
  long i,j,k;
  long migsum;

  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
    for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
      if (opt_mig_bitmatrix[i][j])
      {
        migsum = 0;
        for (k = 0; k < opt_locus_count; ++k)
          migsum += gtree[k]->migcount[i][j];

        if (migsum != stree->migcount_sum[i][j])
        {
          fprintf(stderr, "[DEBUG] Error in function %s\n", __FUNCTION__);
          fatal("migsum_count[%ld][%ld]: %ld   correct: %ld   locus: %ld\n",
                stree->migcount_sum[i][j], migsum, k);
          for (k = 0; k < opt_locus_count; ++k)
            debug_print_gtree(gtree[k]);
        }
        assert(migsum == stree->migcount_sum[i][j]);
      }
}

void debug_print_stree(stree_t * stree)
{
  long i;
  long total_nodes = stree->tip_count+stree->inner_count+stree->hybrid_count;

  printf("Total nodes: %ld\n", total_nodes);

  printf("Label        Node  Child1  Child2  Parent  Linked     Theta       Tau  Leaves  nin[0]  "
         "prop_tau      hphi  htau  has_phi\n");
  for (i = 0; i <total_nodes; ++i)
  {
    char * label;
    if (stree->nodes[i]->label)
      label = xstrdup(stree->nodes[i]->label);
    else
      label = xstrdup("N/A");

    /* shorten label if over 12 characters */
    if (strlen(label) > 12)
    {
      label[9] = '.'; label[10] = '.'; label[11] = '.'; label[12] = 0;
    }

    printf("%-12s", label);
    free(label);
    printf(" %4d", stree->nodes[i]->node_index);

    if (stree->nodes[i]->left)
      printf("  %6d", stree->nodes[i]->left->node_index);
    else
      printf("  %6s", "N/A");
    if (stree->nodes[i]->right)
      printf("  %6d", stree->nodes[i]->right->node_index);
    else
      printf("  %6s", "N/A");
    if (stree->nodes[i]->parent)
      printf("  %6d", stree->nodes[i]->parent->node_index);
    else
      printf("  %6s", "N/A");

    if (stree->nodes[i]->linked_theta)
      printf("  %1s:%2d", "Yes", stree->nodes[i]->linked_theta->node_index);
    else
      printf("  %6s", "No");

    double theta = 0;
    if (stree->nodes[i]->linked_theta)
      theta = stree->nodes[i]->linked_theta->theta;
    else
      theta = stree->nodes[i]->theta;

    if (theta < 0)
      printf("  %8s", "-");
    else
      printf("  %.6f", theta);

    if (stree->nodes[i]->tau)
      printf("  %.6f", stree->nodes[i]->tau);
    else
      printf("  %8s", "-");

    printf("  %6d", stree->nodes[i]->leaves);
    printf("  %6d", stree->nodes[i]->seqin_count[0]);
    printf("  %8d", stree->nodes[i]->prop_tau);
    
    if (opt_msci && stree->nodes[i]->hybrid)
      printf("  %.6f", stree->nodes[i]->hphi);
    else
      printf("  %8s", "-");

    printf("  %4ld", stree->nodes[i]->htau);

    if (opt_msci)
      printf("  %7ld\n", stree->nodes[i]->has_phi);
    else
      printf("\n");
  }
}

void debug_print_gtree(gtree_t * gtree)
{
  long i,j;

  printf("Label        Node-index  Child1-index  Child2-index  Parent-index    Mi  "
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

    if (gtree->nodes[i]->mi && gtree->nodes[i]->mi->count)
      printf("  %4ld", gtree->nodes[i]->mi->count);
    else
      printf("  %4s", "-");

    printf("  %10d", gtree->nodes[i]->pmatrix_index);

    printf(" %5s", gtree->nodes[i]->pop->label);
    printf(" %.12f", gtree->nodes[i]->time);

    if (gtree->root == gtree->nodes[i])
      printf(" (root)");
    printf("\n");
  }
  if (opt_migration)
  {
    long migfound = 0;
    for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
    {
      if (!(gtree->nodes[i]->mi && gtree->nodes[i]->mi->count)) continue;
      migfound++;
    }

    if (migfound)
    {
      printf("\nMigration events:\n");
      printf("Label  Node  #Migs  Migrations (backwards in time)\n");
    }

    for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
    {
      if (!(gtree->nodes[i]->mi && gtree->nodes[i]->mi->count)) continue;

      gnode_t * x = gtree->nodes[i];

      char * label;
      if (x->label)
        label = xstrdup(x->label);
      else
        label = xstrdup("N/A");

      /* shorten label if over 12 characters */
      if (strlen(label) > 6)
      {
        label[4] = '.'; label[5] = '.'; label[6] = 0;
      }

      printf("%-6s", label);
      free(label);
      printf(" %4d", gtree->nodes[i]->node_index);
      printf("  %5ld", x->mi->count);

      printf("  %s -> %s   %f\n",
             x->mi->me[0].source->label,
             x->mi->me[0].target->label,
             x->mi->me[0].time);
      for (j = 1; j < x->mi->count; ++j)
        printf("                    %s -> %s   %f\n",
               x->mi->me[j].source->label,
               x->mi->me[j].target->label,
               x->mi->me[j].time);



    }
  }
}

void debug_print_gtree_detailed(gtree_t * gtree)
{
  long i,j;

  printf("Label        Node-index  Child1-index  Child2-index  Parent-index    Mi  "
         "Pmat-index    Pop            time   Flags\n");
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
  {
    char * label;
    if (gtree->nodes[i]->label)
      label = xstrdup(gtree->nodes[i]->label);
    else
      label = xstrdup("  -");

    /* shorten label if over 12 characters */
    if (strlen(label) > 12)
    {
      label[9] = '.'; label[10] = '.'; label[11] = '.'; label[12] = 0;
    }

    if (gtree->nodes[i]->mark & FLAG_PRUNED)
      printf(ANSI_COLOR_RED);
    printf("%-12s", label);
    free(label);
    char * str_nodeindex = NULL;
    xasprintf(&str_nodeindex,
              gtree->nodes[i] == gtree->root ? "* %d" : "%d",
              gtree->nodes[i]->node_index);
    printf("  %9s", str_nodeindex);
    free(str_nodeindex);

    if (gtree->nodes[i]->left)
      printf("  %12d", gtree->nodes[i]->left->node_index);
    else
      printf("  %12s", "  -");

    if (gtree->nodes[i]->right)
      printf("  %12d", gtree->nodes[i]->right->node_index);
    else
      printf("  %12s", "  -");


    if (gtree->nodes[i]->parent)
      printf("  %12d", gtree->nodes[i]->parent->node_index);
    else
      printf("  %12s", "  -");

    if (gtree->nodes[i]->mi && gtree->nodes[i]->mi->count)
      printf("  %4ld", gtree->nodes[i]->mi->count);
    else
      printf("  %4s", "-");

    printf("  %10d", gtree->nodes[i]->pmatrix_index);

    if (gtree->nodes[i]->pop)
      printf("  %5s", gtree->nodes[i]->pop->label);
    else
      printf("  %5s", "-");

    printf("  %.12f", gtree->nodes[i]->time);
    if (gtree->nodes[i]->mark)
      printf("   ");

    int comma = 0;
    if (gtree->nodes[i]->mark & FLAG_AGE_UPDATE)
    {
      if (comma)
        printf(",");
      printf("AGE_UPDATE");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_POP_UPDATE)
    {
      if (comma)
        printf(",");
      printf("POP_UPDATE");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_BRANCH_UPDATE)
    {
      if (comma)
        printf(",");
      printf("BRANCH_UPDATE");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_MISC)
    {
      if (comma)
        printf(",");
      printf("MISC");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_PRUNED)
    {
      if (comma)
        printf(",");
      printf("PRUNED");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_MIGRATE)
    {
      if (comma)
        printf(",");
      printf("MIGRATE");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_PARTIAL_UPDATE)
    {
      if (comma)
        printf(",");
      printf("PARTIAL_UPDATE");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_RED_LEFT)
    {
      if (comma)
        printf(",");
      printf("RED_LEFT");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_RED_RIGHT)
    {
      if (comma)
        printf(",");
      printf("RED_RIGHT");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_SIMULATE)
    {
      if (comma)
        printf(",");
      printf("SIMULATE");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_PRUNED)
      printf(ANSI_COLOR_RESET);

    printf("\n");
  }
  if (opt_migration)
  {
    long migfound = 0;
    for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
    {
      if (!(gtree->nodes[i]->mi && gtree->nodes[i]->mi->count)) continue;
      migfound++;
    }

    if (migfound)
    {
      printf("\nMigration events:\n");
      printf("Label  Node  #Migs  Migrations (backwards in time)\n");
    }

    for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
    {
      if (!(gtree->nodes[i]->mi && gtree->nodes[i]->mi->count)) continue;

      gnode_t * x = gtree->nodes[i];

      char * label;
      if (x->label)
        label = xstrdup(x->label);
      else
        label = xstrdup("N/A");

      /* shorten label if over 12 characters */
      if (strlen(label) > 6)
      {
        label[4] = '.'; label[5] = '.'; label[6] = 0;
      }

      printf("%-6s", label);
      free(label);
      printf(" %4d", gtree->nodes[i]->node_index);
      printf("  %5ld", x->mi->count);

      printf("  %s -> %s   %f\n", x->mi->me[0].source->label, x->mi->me[0].target->label, x->mi->me[0].time);
      for (j = 1; j < x->mi->count; ++j)
        printf("                    %s -> %s   %f\n", x->mi->me[j].source->label, x->mi->me[j].target->label, x->mi->me[j].time);



    }
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

    if (opt_migration)
    {
      for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
      {
        lout += gtree->migcount[snode->left->node_index][j];
        rout += gtree->migcount[snode->right->node_index][j];

        lout -= gtree->migcount[j][snode->left->node_index];
        rout -= gtree->migcount[j][snode->right->node_index];
      }
    }

    assert(snode->seqin_count[msa_index] == lout+rout);
  }
  snode_t * sroot = stree->root;
  assert(sroot->seqin_count[msa_index] - sroot->event_count[msa_index] == 1);

  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode_t * x = stree->nodes[i];

    for (j = 1; j < x->mb_count; ++j)
      assert(x->migbuffer[j].time != x->migbuffer[j-1].time);
  }

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
static void debug_validate_logpg_theta(stree_t * stree,
                                       gtree_t ** gtree,
                                       locus_t ** locus,
                                       const char * move)
{
  long i,j;
  double logpr, old_logpr;
  const long thread_index = 0;

  logpr = 0;

  assert(!opt_msci);

  for (i = 0; i < opt_locus_count; ++i)
  {
    old_logpr = gtree[i]->logpr;
    if (opt_migration)
      logpr = gtree_logprob_mig(stree,gtree[i],locus[i]->heredity[0],i,thread_index);
    else
      logpr = gtree_logprob(stree,locus[i]->heredity[0],i,thread_index);

    for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
    {
      snode_t * x = stree->nodes[j];

      if (fabs(x->logpr_contrib[i] - x->old_logpr_contrib[i]) > 1e-9)
        fatal("[FATAL-%ld] move: %s locus: %ld snode: %d "
              "logpr_contrib: %f old_contrib: %f\n",
              opt_debug_counter, move, i, j,
              x->logpr_contrib[i], x->old_logpr_contrib[i]);
    }

    if (fabs(logpr - old_logpr) > 1e-9)
      fatal("FATAL-%ld] move: %s locus %ld logpr: %f old_logpr: %f\n",
            opt_debug_counter, move, i, logpr, old_logpr);
  }
}

static void debug_validate_logpg_notheta(stree_t * stree,
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
    debug_validate_logpg_notheta(stree,locus,move);
  else
    debug_validate_logpg_theta(stree,gtree,locus,move);
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

void debug_migration_internals(stree_t * stree, gtree_t * gtree, long msa_index)
{
  long i,j;

  long * migevent_count;
  int * seqin_count;
  long ** migcount;
  snode_t * curpop;
  snode_t * pop;

  /* return; */

  if (!opt_migration) return;

  assert(!gtree->root->parent);

  long total_nodes = stree->tip_count+stree->inner_count;

  migevent_count = (long *)xcalloc((size_t)total_nodes, sizeof(long));
  seqin_count = (int *)xcalloc((size_t)total_nodes, sizeof(int));

  
  migcount = (long **)xcalloc((size_t)total_nodes, sizeof(long *));
  for (i = 0; i < total_nodes; ++i)
    migcount[i] = (long *)xcalloc((size_t)total_nodes, sizeof(long));


  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
  {
    gnode_t * x = gtree->nodes[i];

    if (i < gtree->tip_count)
      seqin_count[x->pop->node_index]++;

    curpop = x->pop;
    if (x->mi)
    {
      for (j = 0; j < x->mi->count; ++j)
      {
        long s = x->mi->me[j].source->node_index;
        long t = x->mi->me[j].target->node_index;

        assert(x->mi->me[j].target->parent);
        assert(x->mi->me[j].target->tau < x->mi->me[j].time);
        assert(x->mi->me[j].target->parent->tau > x->mi->me[j].time);

        migevent_count[s]++;
        migevent_count[t]++;
        migcount[t][s]++;

        if (curpop != x->mi->me[j].source)
        {
          for (pop = curpop->parent; pop != x->mi->me[j].source->parent; pop = pop->parent)
            seqin_count[pop->node_index]++;
        }
        curpop = x->mi->me[j].target;
      }
    }
    snode_t * end = x->parent ? x->parent->pop->parent : NULL;
    for (pop = curpop->parent; pop != end; pop = pop->parent)
      seqin_count[pop->node_index]++;
  }

  /* check */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
  {
    snode_t * snode = stree->nodes[i];
    unsigned int node_index = snode->node_index;

   
    if (snode->seqin_count[msa_index] != seqin_count[node_index])
    {
      printf("\n");
      printf("[DEBUG]: Locus %ld snode %s\n", msa_index, snode->label);
      printf("seqin_count: %d    wrong: %d\n",
             seqin_count[node_index], snode->seqin_count[msa_index]);
      debug_print_gtree(gtree);
      fatal("Exiting");
    }
    //assert(snode->seqin_count[msa_index] == seqin_count[node_index]);
    if (snode->migevent_count[msa_index] != migevent_count[node_index])
    {
      printf("\n");
      printf("[DEBUG]: Locus %ld snode %s\n", msa_index, snode->label);
      printf("migevent_count: %ld    wrong: %ld\n",
             migevent_count[node_index], snode->migevent_count[msa_index]);
      for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
        printf("  node %ld migevent_count %ld  wrong %ld  \n",
               j,migevent_count[j],stree->nodes[j]->migevent_count[msa_index]);
      debug_print_gtree(gtree);
      fatal("Exiting");
    }
    //assert(snode->migevent_count[msa_index] == migevent_count[node_index]);
  }
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    for (j = 0; j < stree->tip_count + stree->inner_count; ++j)
      assert(gtree->migcount[i][j] == migcount[i][j]);

  /* deallocate */
  free(seqin_count);
  free(migevent_count);
  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
    free(migcount[i]);
  free(migcount);
}

void debug_consistency(stree_t * stree, gtree_t ** gtree_list)
{
  long i,j;
  for (i = 0; i < opt_locus_count; ++i)
  {
    gtree_t * gtree = gtree_list[i];
  
    for (j = 0; j < gtree->tip_count; ++j)
    {
      gnode_t * x = gtree->nodes[j];
      if (x->left || x->right)
        fatal("[ERROR] Locus %ld [%d] gnode %ld should be tip but has children",
              i, gtree->tip_count, x->node_index);
    }
    for (j = gtree->tip_count; j < gtree->tip_count + gtree->inner_count; ++j)
    {
      gnode_t * x = gtree->nodes[j];

      if (x->time <= x->pop->tau)
        fatal("[ERROR] Locus %ld [%d] gnode %ld has age %f but its population (%s) starts at %f",
              i, gtree->tip_count, x->node_index, x->time, x->pop->label, x->pop->tau);

      if (x->pop->parent && x->pop->parent->tau <= x->time)
        fatal("[ERROR] Locus %ld [%d] gnode %ld [%s] has age %f but its population's parent (%s) starts at %f",
              i, gtree->tip_count, x->node_index, x->pop->label, x->time, x->pop->parent->label, x->pop->parent->tau);
      
      if (x->parent && x->parent->time <= x->time)
        fatal("[ERROR] Locus %ld [%d] gnode %ld has age %f but its parent has older age %f",
              i, gtree->tip_count, x->node_index, x->time, x->parent->time);

      if (x->left->time >= x->time)
        fatal("[ERROR] Locus %ld [%d] gnode %ld has age %f but its left daughter has older age %f",
              i, gtree->tip_count, x->node_index, x->time, x->left->time);
      if (x->right->time >= x->time)
        fatal("[ERROR] Locus %ld [%d] gnode %ld has age %f but its right daughter has older age %f",
              i, gtree->tip_count, x->node_index, x->time, x->right->time);
    }

    /* chekc if each node is in the correct population */
    for (j = 0; j < gtree->tip_count+gtree->inner_count; ++j)
    {
      gnode_t * x = gtree->nodes[j];
      snode_t * pop = gtree->nodes[j]->pop;

      if (x->time < pop->tau)
      {
        printf("node_index: %d time: %f pop: %s pop->tau: %f\n", x->node_index, x->time, x->pop->label, x->pop->tau);
        debug_print_gtree(gtree);
      }
      assert(x->time >= pop->tau);
      if (pop->parent && pop->parent->tau)
        assert(x->time < pop->parent->tau);
    }
    unsigned int event_count = 0;
    for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
    {
      event_count += stree->nodes[j]->event_count[i];
      dlist_item_t * item = stree->nodes[j]->event[i]->head;
      int pop_event_count = 0;
      while (item)
      {
        gnode_t * x = (gnode_t *)(item->data);
        assert(x->pop == stree->nodes[j]);
        ++pop_event_count;

        item = item->next;
      }
      assert(pop_event_count == stree->nodes[j]->event_count[i]);
    }
    assert(event_count == gtree->inner_count);
  }
}

void debug_write_migs_header(FILE * fp_debug, stree_t * stree)
{
  long j,k;

  fprintf(fp_debug, "Order:");
  for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
    for (k = 0; k < stree->tip_count+stree->inner_count; ++k)
      if (opt_mig_bitmatrix[j][k])
        fprintf(fp_debug, " MC_%s->%s", stree->nodes[j]->label, stree->nodes[k]->label);
  fprintf(fp_debug,"\n");
}
void debug_write_migs(FILE * fp_debug, stree_t * stree, gtree_t ** gtree, long sample)
{
  int print_perlocus = 0;
  int print_sample_number = 0;
  long i,j,k;
  
  long * perlocus  = (long *)xcalloc((size_t)opt_locus_count, sizeof(long));
  long * perthread = (long *)xcalloc((size_t)opt_locus_count, sizeof(long));

  if (print_sample_number)
    fprintf(fp_debug, "%ld\n", sample);
  for (i = 0; i < opt_locus_count; ++i)
  {
    for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
      for (k = 0; k < stree->tip_count+stree->inner_count; ++k)
        if (opt_mig_bitmatrix[j][k])
        {
          if (print_perlocus)
            fprintf(fp_debug," %5ld", gtree[i]->migcount[j][k]);
          perlocus[i] += gtree[i]->migcount[j][k];
        }
    if (print_perlocus)
      fprintf(fp_debug, "  |  %ld\n", perlocus[i]);
  }

  #if 0
  fprintf(fp_debug, "Per-thread migrations:\n");
  #endif
  thread_info_t * ti = threads_ti();
  double mean = 0;
  double stdev = 0;
  for (i = 0; i < opt_threads; ++i)
  {
    thread_info_t * tip = ti+i;
    for (j = tip->locus_first; j < tip->locus_first+tip->locus_count; ++j)
      perthread[i] += perlocus[j];
    if (!i)
      fprintf(fp_debug, "%5ld", perthread[i]);
    else
      fprintf(fp_debug, " %5ld", perthread[i]);
    mean += perthread[i];
  }
  mean /= opt_threads;
  for (i = 0; i < opt_threads; ++i)
    stdev += (perthread[i] - mean)*(perthread[i] - mean);
  stdev = sqrt(stdev / opt_threads);
  
  #if 0
  fprintf(fp_debug, "  | %.2f %.2f\n", mean, stdev);
  #endif

  free(perlocus);
  free(perthread);
  #if 0
  fprintf(fp_debug,"\n\n");
  #else
  fprintf(fp_debug,"\n");
  #endif
}
