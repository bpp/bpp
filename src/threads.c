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

typedef struct qsort_wrapper_s
{
  long index;
  long comp;
} qsort_wrapper_t;

static thread_info_t * ti = NULL;
static pthread_attr_t attr;

static int cb_asc_comp(const void * x, const void * y)
{
  const qsort_wrapper_t * a = *(const qsort_wrapper_t **)x;
  const qsort_wrapper_t * b = *(const qsort_wrapper_t **)y;

  return b->comp - a->comp;
}

static int cb_desc_comp(const void * x, const void * y)
{
  const qsort_wrapper_t * a = *(const qsort_wrapper_t **)x;
  const qsort_wrapper_t * b = *(const qsort_wrapper_t **)y;

  return a->comp - b->comp;
}

#if (defined(__linux__) && !defined(DISABLE_COREPIN))
static void pin_to_core(long t)
{
  cpu_set_t cpuset;

  CPU_ZERO(&cpuset);
  CPU_SET(t,&cpuset);

  if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset))
    fatal("Error while pinning thread to core. "
          "Probably used more threads than available cores?");
}
#endif

static void * threads_worker(void * vp)
{
  long t = (long)vp;
  thread_info_t * tip = ti + t;

#if (defined(__linux__) && !defined(DISABLE_COREPIN))
  if (opt_corepin)
    pin_to_core((opt_threads_start-1)+(t*opt_threads_step));
#endif

  pthread_mutex_lock(&tip->mutex);

  /* loop until signalled to quit */
  while (tip->work >= 0)
  {
    /* wait for work available */
    if (tip->work == 0)
      pthread_cond_wait(&tip->cond, &tip->mutex);

    if (tip->work > 0)
    {
      /* work work! */
      switch (tip->work)
      {
        case THREAD_WORK_GTAGE:
          gtree_propose_ages_parallel(tip->td.locus,
                                      tip->td.gtree,
                                      tip->td.stree,
                                      tip->locus_first,
                                      tip->locus_count,
                                      t,
                                      &tip->td.proposals,
                                      &tip->td.accepted);
          break;
        case THREAD_WORK_GTSPR:
          gtree_propose_spr_parallel(tip->td.locus,
                                     tip->td.gtree,
                                     tip->td.stree,
                                     tip->locus_first,
                                     tip->locus_count,
                                     t,
                                     &tip->td.proposals,
                                     &tip->td.accepted);
          break;
        case THREAD_WORK_TAU:
          propose_tau_update_gtrees(tip->td.locus,
                                    tip->td.gtree,
                                    tip->td.stree,
                                    tip->td.snode,
                                    tip->td.oldage,
                                    tip->td.minage,
                                    tip->td.maxage,
                                    tip->td.minfactor,
                                    tip->td.maxfactor,
                                    tip->locus_first,
                                    tip->locus_count,
                                    tip->td.affected,
                                    tip->td.paffected_count,
                                    &tip->td.count_above,
                                    &tip->td.count_below,
                                    &tip->td.logl_diff,
                                    &tip->td.logpr_diff,
                                    t);
          break;
        case THREAD_WORK_TAU_MIG:
          propose_tau_update_gtrees_mig(tip->td.locus,
                                        tip->td.gtree,
                                        tip->td.stree,
                                        tip->td.snode,
                                        tip->td.oldage,
                                        tip->td.minage,
                                        tip->td.maxage,
                                        tip->td.minfactor,
                                        tip->td.maxfactor,
                                        tip->locus_first,
                                        tip->locus_count,
                                        tip->td.affected,
                                        tip->td.paffected_count,
                                        &tip->td.mig_reject,
                                        &tip->td.count_above,
                                        &tip->td.count_below,
                                        &tip->td.logl_diff,
                                        &tip->td.logpr_diff,
                                        t);
          break;
        case THREAD_WORK_MIXING:
          prop_mixing_update_gtrees(tip->td.locus,
                                    tip->td.gtree,
                                    tip->td.stree,
                                    tip->locus_first,
                                    tip->locus_count,
                                    tip->td.c,
                                    t,
                                    &tip->td.lnacceptance,
                                    NULL);
          break;
        case THREAD_WORK_ALPHA:
          locus_propose_alpha_parallel(tip->td.stree,
                                       tip->td.locus,
                                       tip->td.gtree,
                                       tip->locus_first,
                                       tip->locus_count,
                                       t,
                                       &tip->td.proposals,
                                       &tip->td.accepted);
          break;
        case THREAD_WORK_RATES:
          locus_propose_qrates_parallel(tip->td.stree,
                                        tip->td.locus,
                                        tip->td.gtree,
                                        tip->locus_first,
                                        tip->locus_count,
                                        t,
                                        &tip->td.proposals,
                                        &tip->td.accepted);
          break;
        case THREAD_WORK_FREQS:
          locus_propose_freqs_parallel(tip->td.stree,
                                       tip->td.locus,
                                       tip->td.gtree,
                                       tip->locus_first,
                                       tip->locus_count,
                                       t,
                                       &tip->td.proposals,
                                       &tip->td.accepted);
          break;
        case THREAD_WORK_BRATE:
          prop_branch_rates_parallel(tip->td.gtree,
                                     tip->td.stree,
                                     tip->td.locus,
                                     tip->locus_first,
                                     tip->locus_count,
                                     t,
                                     &tip->td.proposals,
                                     &tip->td.accepted);
          break;
        default:
          fatal("Unknown work function assigned to thread worker %ld", t);
                             
      }
        
      tip->work = 0;
      pthread_cond_signal(&tip->cond);
    }
  }
  pthread_mutex_unlock(&tip->mutex);

  pthread_exit(NULL);
}

void threads_pin_master()
{
  #if (defined(__linux__) && !defined(DISABLE_COREPIN))
  if (opt_corepin)
    pin_to_core(opt_threads_start-1);
  #endif
}

thread_info_t * threads_ti()
{
  return ti;
}

/* this is only used when resuming from a checkpoint */
void threads_set_ti(thread_info_t * tip)
{
  ti = tip;
}

static void load_balance_none(msa_t ** msa_list)
{
  long t;
  long loci_per_thread = opt_locus_count / opt_threads;
  long loci_remaining = opt_locus_count % opt_threads;
  long loci_start = 0;

  /* static load allocation */
  /* TODO: We discussed with Ziheng better ways to allocate loci on threads. One
     idea was to use some weights on the site count of each locus, i.e.
     sites^{2/3} as we need to account for the MSC density computation time and
     not only on the phylogenetic likelihood.

     For now we only distrubute an equal number of loci to each thread and
     cyclically distribute the overflow */

  for (t = 0; t < opt_threads; ++t)
  {
    thread_info_t * tip = ti + t;
    tip->work = 0;

    /* allocate loci for thread t */
    tip->locus_first = loci_start;
    tip->locus_count = loci_per_thread + (loci_remaining > 0 ? 1 : 0);

    if (loci_remaining)
      --loci_remaining;
    loci_start += tip->locus_count;
  }
}

static long * load_balance_zigzag(msa_t ** msa_list)
{
  long i,t;
  long core, increment;
  long * assign;
  long * shuffle_indices;
  qsort_wrapper_t ** loadi;
  msa_t ** reorder_msa;

  loadi = (qsort_wrapper_t **)xmalloc((size_t)opt_locus_count *
                                      sizeof(qsort_wrapper_t *));
  for (i = 0; i < opt_locus_count; ++i)
  {
    loadi[i] = (qsort_wrapper_t *)xmalloc(sizeof(qsort_wrapper_t));
    loadi[i]->index = i;
    loadi[i]->comp  = msa_list[i]->count * msa_list[i]->length;  /* load */
  }
  qsort(loadi,opt_locus_count,sizeof(loadi),cb_asc_comp);

  /* zig-zag distribution */
  core = 0;
  increment = 1;
  assign = (long *)xmalloc((size_t)opt_locus_count * sizeof(long));
  for (i = 0; i < opt_locus_count; ++i)
  {
    assign[i] = core;
    core += increment;

    if (core == opt_threads)
    {
      increment = -1;
      core = opt_threads - 1;
    }
    else if (core == -1)
    {
      increment = 1;
      core = 0;
    }
  }

  /* sort loci by assigned thread index */
  for (i = 0; i < opt_locus_count; ++i)
    loadi[i]->comp = assign[i];
  qsort(loadi,opt_locus_count,sizeof(loadi),cb_desc_comp);

  /* reorder loci and gene tree arrays according to thread assignments */
  reorder_msa = (msa_t **)xmalloc((size_t)opt_locus_count * sizeof(msa_t *));
  for (i = 0; i < opt_locus_count; ++i)
    reorder_msa[i] = msa_list[loadi[i]->index];

  shuffle_indices = (long *)xmalloc((size_t)opt_locus_count * sizeof(long));

  for (i = 0; i < opt_locus_count; ++i)
  {
    /* keep the shuffling order for other arrays that need to be reordered */
    shuffle_indices[i] = loadi[i]->index;

    /* reoder msa_list */
    msa_list[i] = reorder_msa[i];
  }

  free(reorder_msa);
  for (i = 0; i < opt_locus_count; ++i)
    free(loadi[i]);
  free(loadi);

  /* get number of loci assigned to each thread */
  long * lc_thread = (long *)xcalloc((size_t)opt_threads, sizeof(long));
  for (i = 0; i < opt_locus_count; ++i)
    lc_thread[assign[i]]++;
  free(assign);

  long loci_start = 0;

  /* init and create worker threads */
  for (t = 0; t < opt_threads; ++t)
  {
    thread_info_t * tip = ti + t;
    tip->work = 0;

    /* allocate loci for thread t */
    tip->locus_first = loci_start;
    tip->locus_count = lc_thread[t];

    loci_start += tip->locus_count;
  }
  free(lc_thread);
  return shuffle_indices;
}

long * threads_load_balance(msa_t ** msa_list)
{
  long * shuffle_indices = NULL;

  /* allocate memory for thread info */
  ti = (thread_info_t *)xmalloc((size_t)opt_threads * sizeof(thread_info_t));

  if (opt_load_balance == BPP_LB_ZIGZAG)
    shuffle_indices = load_balance_zigzag(msa_list);
  else
    load_balance_none(msa_list);

  return shuffle_indices;
}

void threads_lb_stats(locus_t ** locus, FILE * fp_out)
{
  int ls_digits = 0;
  int le_digits = 0;
  int p_digits  = 0;
  int s_digits  = 0;
  int l_digits  = 0;
  int t_digits  = 0;
  long i,t,n;
  long lindex;

  long * patterns = (long *)xcalloc((size_t)opt_threads, sizeof(long));
  long * seqs = (long *)xcalloc((size_t)opt_threads, sizeof(long));
  long * load = (long *)xcalloc((size_t)opt_threads, sizeof(long));

  /* max number of digits for threads */
  t_digits = (int)(floor(log10(opt_threads)+1));

  for (t = 0; t < opt_threads; ++t)
  {
    thread_info_t * tip = ti + t;

    for (i = 0; i < tip->locus_count; ++i)
    {
      lindex    = tip->locus_first+i;
      patterns[t] += locus[lindex]->sites;
      seqs[t]     += locus[lindex]->tips;
      load[t]     += locus[lindex]->sites * locus[lindex]->tips;
    }

    /* locus start max digits */
    n = (int)(floor(log10(tip->locus_first+1)+1));
    if (n > ls_digits)
      ls_digits = n;

    /* locus end max digits */
    n = (int)(floor(log10(tip->locus_first+1+tip->locus_count)+1));
    if (n > le_digits)
      le_digits = n;

    /* patterns max digits */
    n = (int)(floor(log10(patterns[t])+1));
    if (n > p_digits)
      p_digits = n;

    /* sequences max digits */
    n = (int)(floor(log10(seqs[t])+1));
    if (n > s_digits)
      s_digits = n;

    /* load max digits */
    n = (int)(floor(log10(seqs[t])+1));
    if (n > l_digits)
      l_digits = n;
  }

  fprintf(stdout, "\nDistributing workload to threads:\n");
  fprintf(fp_out, "\nDistributing workload to threads:\n");
  for (t = 0; t < opt_threads; ++t)
  {
    thread_info_t * tip = ti + t;

    fprintf(stdout,
            " Thread %*ld : loci [%*ld - %*ld), Patterns/Seqs/Load : %*ld / %*ld / %*ld\n",
            t_digits, t,
            ls_digits, tip->locus_first+1,
            le_digits, tip->locus_first+1+tip->locus_count,
            p_digits, patterns[t],
            s_digits, seqs[t],
            l_digits, load[t]);
    fprintf(fp_out,
            " Thread %*ld : loci [%*ld - %*ld), Patterns/Seqs/Load : %*ld / %*ld / %*ld\n",
            t_digits, t,
            ls_digits, tip->locus_first+1,
            le_digits, tip->locus_first+1+tip->locus_count,
            p_digits, patterns[t],
            s_digits, seqs[t],
            l_digits, load[t]);
  }

  free(patterns);
  free(seqs);
  free(load);
}

void threads_init()
{
  long t;

  if (opt_threads > opt_locus_count)
    fatal("The number of threads cannot be greater than the number of loci");
  assert(opt_threads <= opt_locus_count);

#if (defined(__linux__) && !defined(DISABLE_COREPIN))
  if (opt_corepin)
    pin_to_core(opt_threads_start-1);
#endif

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  if (!ti)
    fatal("Internal error - call load balance routine");

  /* init and create worker threads */
  for (t = 0; t < opt_threads; ++t)
  {
    thread_info_t * tip = ti + t;
    tip->work = 0;
    tip->td.locus = NULL;
    tip->td.gtree = NULL;
    tip->td.stree = NULL;

    pthread_mutex_init(&tip->mutex, NULL);
    pthread_cond_init(&tip->cond, NULL);
    if (pthread_create(&tip->thread, &attr, threads_worker, (void *)(long)t))
      fatal("Cannot create thread");
  }
}

void threads_wakeup(int work_type, thread_data_t * data)
{
  long t; 

  /* dynamic load distribution */
  /* Currently we do not do any dynamic load distribution. We assign a static
     workload at initialization, which then never changes */


  for (t = 0; t < opt_threads; ++t)
  {
    thread_info_t * tip = ti + t;

    pthread_mutex_lock(&tip->mutex);

    memcpy(&tip->td,data,sizeof(thread_data_t));

    tip->work = work_type;

    pthread_cond_signal(&tip->cond);
    pthread_mutex_unlock(&tip->mutex);
  }

  /* wait for threads to finish their work */
  for (t = 0; t < opt_threads; ++t)
  {
    thread_info_t * tip = ti+t;

    pthread_mutex_lock(&tip->mutex);
    while (tip->work > 0)
      pthread_cond_wait(&tip->cond,&tip->mutex);

    pthread_mutex_unlock(&tip->mutex);
  }

  if (work_type == THREAD_WORK_GTAGE ||
      work_type == THREAD_WORK_GTSPR ||
      work_type == THREAD_WORK_ALPHA ||
      work_type == THREAD_WORK_RATES ||
      work_type == THREAD_WORK_FREQS ||
      work_type == THREAD_WORK_BRATE)
  {
    long proposals = 0;
    long accepted = 0;
    for (t = 0; t < opt_threads; ++t)
    {
      thread_info_t * tip = ti+t;
      proposals += tip->td.proposals;
      accepted  += tip->td.accepted;
    }

    data->proposals = proposals;
    data->accepted = accepted;
  }
  else if (work_type == THREAD_WORK_TAU)
  {
    data->count_above = 0;
    data->count_below = 0;
    data->logl_diff = 0;
    data->logpr_diff = 0;
    for (t = 0; t < opt_threads; ++t)
    {
      thread_info_t * tip = ti+t;

      data->count_above += tip->td.count_above;
      data->count_below += tip->td.count_below;
      data->logl_diff   += tip->td.logl_diff;
      data->logpr_diff  += tip->td.logpr_diff;
    }
  }
  else if (work_type == THREAD_WORK_TAU_MIG)
  {
    data->count_above = 0;
    data->count_below = 0;
    data->logl_diff = 0;
    data->logpr_diff = 0;
    data->mig_reject = 0;
    for (t = 0; t < opt_threads; ++t)
    {
      thread_info_t * tip = ti+t;

      if (tip->td.mig_reject)
      {
        data->mig_reject = 1;
        break;
      }

      data->count_above += tip->td.count_above;
      data->count_below += tip->td.count_below;
      data->logl_diff   += tip->td.logl_diff;
      data->logpr_diff  += tip->td.logpr_diff;
    }
  }
  else if (work_type == THREAD_WORK_MIXING)
  {
    data->lnacceptance = 0;
    for (t = 0; t < opt_threads; ++t)
    {
      thread_info_t * tip = ti+t;
      data->lnacceptance += tip->td.lnacceptance;
    }
  }
  else
    assert(0);
}

void threads_exit()
{
  long t;

  for (t = 0; t < opt_threads; ++t)
  {
    thread_info_t * tip = ti + t;

    /* tell worker to quit */
    pthread_mutex_lock(&tip->mutex);
    tip->work = -1;
    pthread_cond_signal(&tip->cond);
    pthread_mutex_unlock(&tip->mutex);

    /* wait for worker to quit */
    if (pthread_join(tip->thread, 0))
      fatal("Cannot join thread");

    pthread_cond_destroy(&tip->cond);
    pthread_mutex_destroy(&tip->mutex);
  }

  free(ti);
  pthread_attr_destroy(&attr);
}
