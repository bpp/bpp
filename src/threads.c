/*
    Copyright (C) 2016-2025 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

/* Dynamic load balancing work queue with chunking */
typedef struct work_queue_s
{
  volatile long next_locus;      /* next locus index to process */
  volatile long completed_count; /* number of completed chunks */
  long total_loci;               /* total number of loci */
  long chunk_size;               /* number of loci per chunk */
  long total_chunks;             /* total number of chunks */
  pthread_mutex_t mutex;         /* mutex for queue access */
  pthread_cond_t done_cond;      /* condition for completion */
  pthread_barrier_t barrier;     /* barrier for synchronized start */
  int barrier_initialized;       /* flag to track barrier state */
} work_queue_t;

static work_queue_t wq = {0, 0, 0, 1, 0, PTHREAD_MUTEX_INITIALIZER, PTHREAD_COND_INITIALIZER, {0}, 0};

/* Get next chunk from work queue, returns chunk start index and count.
   Returns -1 in *chunk_start if queue is empty. */
static void workqueue_get_chunk(long * chunk_start, long * chunk_count)
{
  pthread_mutex_lock(&wq.mutex);
  if (wq.next_locus < wq.total_loci)
  {
    *chunk_start = wq.next_locus;
    *chunk_count = wq.chunk_size;

    /* Don't exceed total loci */
    if (wq.next_locus + wq.chunk_size > wq.total_loci)
      *chunk_count = wq.total_loci - wq.next_locus;

    wq.next_locus += *chunk_count;
  }
  else
  {
    *chunk_start = -1;
    *chunk_count = 0;
  }
  pthread_mutex_unlock(&wq.mutex);
}

/* Signal that a chunk has been completed */
static void workqueue_complete_chunk(void)
{
  pthread_mutex_lock(&wq.mutex);
  wq.completed_count++;
  if (wq.completed_count >= wq.total_chunks)
  {
    pthread_cond_broadcast(&wq.done_cond);
  }
  pthread_mutex_unlock(&wq.mutex);
}

/* Initialize work queue for a new parallel phase.
   Chunk size is calculated to give ~4 chunks per thread for good balancing
   while reducing mutex overhead. Minimum chunk size is 1. */
static void workqueue_init(long total_loci)
{
  long chunks_per_thread = 4;  /* target chunks per thread for balancing */
  long target_chunks;

  pthread_mutex_lock(&wq.mutex);
  wq.next_locus = 0;
  wq.completed_count = 0;
  wq.total_loci = total_loci;

  /* Calculate chunk size: aim for ~4 chunks per thread */
  target_chunks = opt_threads * chunks_per_thread;
  wq.chunk_size = (total_loci + target_chunks - 1) / target_chunks;
  if (wq.chunk_size < 1) wq.chunk_size = 1;

  /* Calculate actual number of chunks */
  wq.total_chunks = (total_loci + wq.chunk_size - 1) / wq.chunk_size;

  /* Destroy old barrier if it exists, then create new one */
  if (wq.barrier_initialized)
    pthread_barrier_destroy(&wq.barrier);

  pthread_barrier_init(&wq.barrier, NULL, (unsigned)opt_threads);
  wq.barrier_initialized = 1;

  pthread_mutex_unlock(&wq.mutex);
}

/* Wait at barrier for all threads to synchronize before starting work */
static void workqueue_barrier_wait(void)
{
  pthread_barrier_wait(&wq.barrier);
}

/* Wait for all work to complete */
static void workqueue_wait(void)
{
  pthread_mutex_lock(&wq.mutex);
  while (wq.completed_count < wq.total_chunks)
  {
    pthread_cond_wait(&wq.done_cond, &wq.mutex);
  }
  pthread_mutex_unlock(&wq.mutex);
}

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
  struct timespec ts_start, ts_end;
  int work_type;
  long local_proposals, local_accepted;
  unsigned int local_count_above, local_count_below;
  double local_logl_diff, local_logpr_diff, local_lnacceptance;
  long local_mig_reject;

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
      /* start timing (for static mode) */
      work_type = tip->work;

      /* Dynamic load balancing mode: pull chunks from shared queue */
      if (opt_load_balance == BPP_LB_DYNAMIC)
      {
        long chunk_start, chunk_count;

        /* Initialize thread-local accumulators */
        tip->td.proposals = 0;
        tip->td.accepted = 0;
        tip->td.count_above = 0;
        tip->td.count_below = 0;
        tip->td.logl_diff = 0;
        tip->td.logpr_diff = 0;
        tip->td.lnacceptance = 0;
        tip->td.mig_reject = 0;

        /* Wait for all threads to be ready before starting */
        workqueue_barrier_wait();

        /* Start timing AFTER barrier - don't count synchronization wait */
        clock_gettime(CLOCK_MONOTONIC, &ts_start);

        /* Pull chunks from queue until empty */
        workqueue_get_chunk(&chunk_start, &chunk_count);
        while (chunk_start >= 0)
        {
          /* Reset per-chunk accumulators */
          local_proposals = 0;
          local_accepted = 0;
          local_count_above = 0;
          local_count_below = 0;
          local_logl_diff = 0;
          local_logpr_diff = 0;
          local_lnacceptance = 0;
          local_mig_reject = 0;

          /* Process chunk based on work type */
          switch (work_type)
          {
            case THREAD_WORK_GTAGE:
              gtree_propose_ages_parallel(tip->td.locus,
                                          tip->td.gtree,
                                          tip->td.stree,
                                          chunk_start,
                                          chunk_count,
                                          t,
                                          &local_proposals,
                                          &local_accepted);
              tip->td.proposals += local_proposals;
              tip->td.accepted += local_accepted;
              break;
            case THREAD_WORK_GTSPR:
              gtree_propose_spr_parallel(tip->td.locus,
                                         tip->td.gtree,
                                         tip->td.stree,
                                         chunk_start,
                                         chunk_count,
                                         t,
                                         &local_proposals,
                                         &local_accepted);
              tip->td.proposals += local_proposals;
              tip->td.accepted += local_accepted;
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
                                        chunk_start,
                                        chunk_count,
                                        tip->td.affected,
                                        tip->td.paffected_count,
                                        &local_count_above,
                                        &local_count_below,
                                        &local_logl_diff,
                                        &local_logpr_diff,
                                        t);
              tip->td.count_above += local_count_above;
              tip->td.count_below += local_count_below;
              tip->td.logl_diff += local_logl_diff;
              tip->td.logpr_diff += local_logpr_diff;
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
                                            chunk_start,
                                            chunk_count,
                                            tip->td.affected,
                                            tip->td.paffected_count,
                                            &local_mig_reject,
                                            &local_count_above,
                                            &local_count_below,
                                            &local_logl_diff,
                                            &local_logpr_diff,
                                            t);
              if (local_mig_reject)
                tip->td.mig_reject = 1;
              tip->td.count_above += local_count_above;
              tip->td.count_below += local_count_below;
              tip->td.logl_diff += local_logl_diff;
              tip->td.logpr_diff += local_logpr_diff;
              break;
            case THREAD_WORK_MIXING:
              prop_mixing_update_gtrees(tip->td.locus,
                                        tip->td.gtree,
                                        tip->td.stree,
                                        chunk_start,
                                        chunk_count,
                                        tip->td.c,
                                        t,
                                        &local_lnacceptance);
              tip->td.lnacceptance += local_lnacceptance;
              break;
            case THREAD_WORK_ALPHA:
              locus_propose_alpha_parallel(tip->td.stree,
                                           tip->td.locus,
                                           tip->td.gtree,
                                           chunk_start,
                                           chunk_count,
                                           t,
                                           &local_proposals,
                                           &local_accepted);
              tip->td.proposals += local_proposals;
              tip->td.accepted += local_accepted;
              break;
            case THREAD_WORK_RATES:
              locus_propose_qrates_parallel(tip->td.stree,
                                            tip->td.locus,
                                            tip->td.gtree,
                                            chunk_start,
                                            chunk_count,
                                            t,
                                            &local_proposals,
                                            &local_accepted);
              tip->td.proposals += local_proposals;
              tip->td.accepted += local_accepted;
              break;
            case THREAD_WORK_FREQS:
              locus_propose_freqs_parallel(tip->td.stree,
                                           tip->td.locus,
                                           tip->td.gtree,
                                           chunk_start,
                                           chunk_count,
                                           t,
                                           &local_proposals,
                                           &local_accepted);
              tip->td.proposals += local_proposals;
              tip->td.accepted += local_accepted;
              break;
            case THREAD_WORK_BRATE:
              prop_branch_rates_parallel(tip->td.gtree,
                                         tip->td.stree,
                                         tip->td.locus,
                                         chunk_start,
                                         chunk_count,
                                         t,
                                         &local_proposals,
                                         &local_accepted);
              tip->td.proposals += local_proposals;
              tip->td.accepted += local_accepted;
              break;
            default:
              fatal("Unknown work function assigned to thread worker %ld", t);
          }

          /* Signal this chunk is complete and get next */
          workqueue_complete_chunk();
          workqueue_get_chunk(&chunk_start, &chunk_count);
        }
      }
      else
      {
        /* Static load balancing mode: process assigned range */
        clock_gettime(CLOCK_MONOTONIC, &ts_start);

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
                                      &tip->td.lnacceptance);
                                      #if 0
                                      NULL);
                                      #endif
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
      }

      /* end timing and record */
      clock_gettime(CLOCK_MONOTONIC, &ts_end);
      tip->elapsed_us = (ts_end.tv_sec - ts_start.tv_sec) * 1000000.0 +
                        (ts_end.tv_nsec - ts_start.tv_nsec) / 1000.0;
      tip->total_time_us[work_type] += tip->elapsed_us;
      tip->call_count[work_type]++;

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

/*
 * Weighted load balancing with improved metric for MSC/MSC-M models.
 *
 * The load metric accounts for:
 * 1. Phylogenetic likelihood: scales with sites × inner_nodes
 * 2. MSC density computation: scales approximately with sites^(2/3)
 * 3. Tree complexity: scales with tips (which determines inner_nodes = tips-1)
 * 4. Migration factor: adds weight for MSC-M models
 *
 * Formula: load = sites^(2/3) × (tips - 1) × mig_factor
 *
 * The 2/3 exponent comes from the observation that MSC density computation
 * doesn't scale linearly with sites (as noted in TODO comment above).
 * For MSC-M, migration events add significant overhead.
 */
static long * load_balance_weighted(msa_t ** msa_list)
{
  long i, t;
  long core, increment;
  long * assign;
  long * shuffle_indices;
  qsort_wrapper_t ** loadi;
  msa_t ** reorder_msa;
  double mig_factor;

  /* migration factor: MSC-M has ~1.5x more work due to migration proposals */
  mig_factor = opt_migration ? 1.5 : 1.0;

  loadi = (qsort_wrapper_t **)xmalloc((size_t)opt_locus_count *
                                      sizeof(qsort_wrapper_t *));
  for (i = 0; i < opt_locus_count; ++i)
  {
    loadi[i] = (qsort_wrapper_t *)xmalloc(sizeof(qsort_wrapper_t));
    loadi[i]->index = i;

    /* weighted load metric:
     * - sites^(2/3) accounts for MSC density scaling (not linear with sites)
     * - (tips - 1) = inner_count approximates tree complexity
     * - mig_factor adds weight for migration models
     */
    double sites = (double)msa_list[i]->length;
    double tips = (double)msa_list[i]->count;
    double inner_count = tips - 1.0;
    if (inner_count < 1.0) inner_count = 1.0;

    double load = pow(sites, 2.0/3.0) * inner_count * mig_factor;
    loadi[i]->comp = (long)(load + 0.5);  /* round to nearest integer */
  }
  qsort(loadi, opt_locus_count, sizeof(loadi), cb_asc_comp);

  /* zig-zag distribution (same as load_balance_zigzag) */
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
  qsort(loadi, opt_locus_count, sizeof(loadi), cb_desc_comp);

  /* reorder loci and gene tree arrays according to thread assignments */
  reorder_msa = (msa_t **)xmalloc((size_t)opt_locus_count * sizeof(msa_t *));
  for (i = 0; i < opt_locus_count; ++i)
    reorder_msa[i] = msa_list[loadi[i]->index];

  shuffle_indices = (long *)xmalloc((size_t)opt_locus_count * sizeof(long));

  for (i = 0; i < opt_locus_count; ++i)
  {
    /* keep the shuffling order for other arrays that need to be reordered */
    shuffle_indices[i] = loadi[i]->index;

    /* reorder msa_list */
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

  /* init thread info */
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

  if (opt_load_balance == BPP_LB_WEIGHTED)
    shuffle_indices = load_balance_weighted(msa_list);
  else if (opt_load_balance == BPP_LB_ZIGZAG)
    shuffle_indices = load_balance_zigzag(msa_list);
  else if (opt_load_balance == BPP_LB_DYNAMIC)
  {
    /* Dynamic mode: threads pull from work queue, static assignment not used.
       Initialize thread info with placeholder values. */
    load_balance_none(msa_list);
  }
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
  int w_digits  = 0;
  int t_digits  = 0;
  long i,t,n;
  long lindex;
  double mig_factor;

  long * patterns = (long *)xcalloc((size_t)opt_threads, sizeof(long));
  long * seqs = (long *)xcalloc((size_t)opt_threads, sizeof(long));
  long * load = (long *)xcalloc((size_t)opt_threads, sizeof(long));
  long * wload = (long *)xcalloc((size_t)opt_threads, sizeof(long));

  /* migration factor for weighted load */
  mig_factor = opt_migration ? 1.5 : 1.0;

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

      /* weighted load: sites^(2/3) * (tips-1) * mig_factor */
      double sites = (double)locus[lindex]->sites;
      double inner = (double)(locus[lindex]->tips - 1);
      if (inner < 1.0) inner = 1.0;
      wload[t] += (long)(pow(sites, 2.0/3.0) * inner * mig_factor + 0.5);
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
    n = (int)(floor(log10(load[t])+1));
    if (n > l_digits)
      l_digits = n;

    /* weighted load max digits */
    n = (int)(floor(log10(wload[t])+1));
    if (n > w_digits)
      w_digits = n;
  }

  /* print load balance mode */
  const char * lb_mode = "none";
  if (opt_load_balance == BPP_LB_ZIGZAG) lb_mode = "zigzag";
  else if (opt_load_balance == BPP_LB_WEIGHTED) lb_mode = "weighted";
  else if (opt_load_balance == BPP_LB_DYNAMIC) lb_mode = "dynamic";

  fprintf(stdout, "\nDistributing workload to threads (mode: %s):\n", lb_mode);
  fprintf(fp_out, "\nDistributing workload to threads (mode: %s):\n", lb_mode);
  for (t = 0; t < opt_threads; ++t)
  {
    thread_info_t * tip = ti + t;

    fprintf(stdout,
            " Thread %*ld : loci [%*ld - %*ld), Patt/Seqs/Load/WLoad : %*ld / %*ld / %*ld / %*ld\n",
            t_digits, t,
            ls_digits, tip->locus_first+1,
            le_digits, tip->locus_first+1+tip->locus_count,
            p_digits, patterns[t],
            s_digits, seqs[t],
            l_digits, load[t],
            w_digits, wload[t]);
    fprintf(fp_out,
            " Thread %*ld : loci [%*ld - %*ld), Patt/Seqs/Load/WLoad : %*ld / %*ld / %*ld / %*ld\n",
            t_digits, t,
            ls_digits, tip->locus_first+1,
            le_digits, tip->locus_first+1+tip->locus_count,
            p_digits, patterns[t],
            s_digits, seqs[t],
            l_digits, load[t],
            w_digits, wload[t]);
  }

  free(wload);

  free(patterns);
  free(seqs);
  free(load);
}

void threads_init()
{
  long t, i;

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

    /* initialize timing instrumentation */
    tip->elapsed_us = 0;
    for (i = 0; i < THREAD_WORK_COUNT; ++i)
    {
      tip->total_time_us[i] = 0;
      tip->call_count[i] = 0;
    }

    pthread_mutex_init(&tip->mutex, NULL);
    pthread_cond_init(&tip->cond, NULL);
    if (pthread_create(&tip->thread, &attr, threads_worker, (void *)(long)t))
      fatal("Cannot create thread");
  }
}

void threads_wakeup(int work_type, thread_data_t * data)
{
  long t;

  /* Initialize work queue for dynamic load balancing */
  if (opt_load_balance == BPP_LB_DYNAMIC)
  {
    workqueue_init(opt_locus_count);
  }

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

void threads_reset_timing()
{
  long t, i;

  if (!ti) return;

  for (t = 0; t < opt_threads; ++t)
  {
    thread_info_t * tip = ti + t;
    tip->elapsed_us = 0;
    for (i = 0; i < THREAD_WORK_COUNT; ++i)
    {
      tip->total_time_us[i] = 0;
      tip->call_count[i] = 0;
    }
  }
}

static const char * work_type_name(int work_type)
{
  switch (work_type)
  {
    case THREAD_WORK_GTAGE:   return "GTAGE";
    case THREAD_WORK_GTSPR:   return "GTSPR";
    case THREAD_WORK_TAU:     return "TAU";
    case THREAD_WORK_TAU_MIG: return "TAU_MIG";
    case THREAD_WORK_MIXING:  return "MIXING";
    case THREAD_WORK_ALPHA:   return "ALPHA";
    case THREAD_WORK_RATES:   return "RATES";
    case THREAD_WORK_FREQS:   return "FREQS";
    case THREAD_WORK_BRATE:   return "BRATE";
    default:                  return "UNKNOWN";
  }
}

void threads_print_timing_summary(FILE * fp_out, int detailed)
{
  long t, i;
  double total_time[THREAD_WORK_COUNT] = {0};
  double min_time[THREAD_WORK_COUNT];
  double max_time[THREAD_WORK_COUNT];
  double avg_time[THREAD_WORK_COUNT];
  long   total_calls[THREAD_WORK_COUNT] = {0};
  double grand_total = 0;
  double overall_min_time = 0, overall_max_time = 0, overall_total_time = 0;

  if (!ti || opt_threads <= 1) return;

  /* initialize min/max */
  for (i = 0; i < THREAD_WORK_COUNT; ++i)
  {
    min_time[i] = 1e30;
    max_time[i] = 0;
  }

  /* aggregate statistics */
  for (t = 0; t < opt_threads; ++t)
  {
    thread_info_t * tip = ti + t;
    double thread_total = 0;

    for (i = 1; i < THREAD_WORK_COUNT; ++i)  /* skip 0 = no work */
    {
      if (tip->call_count[i] > 0)
      {
        total_time[i] += tip->total_time_us[i];
        total_calls[i] += tip->call_count[i];
        if (tip->total_time_us[i] < min_time[i])
          min_time[i] = tip->total_time_us[i];
        if (tip->total_time_us[i] > max_time[i])
          max_time[i] = tip->total_time_us[i];
        thread_total += tip->total_time_us[i];
      }
    }
    grand_total += thread_total;
    overall_total_time += thread_total;
    if (t == 0 || thread_total < overall_min_time)
      overall_min_time = thread_total;
    if (t == 0 || thread_total > overall_max_time)
      overall_max_time = thread_total;
  }

  if (grand_total == 0) return;

  /* compute averages */
  for (i = 1; i < THREAD_WORK_COUNT; ++i)
  {
    if (total_calls[i] > 0)
      avg_time[i] = total_time[i] / opt_threads;
    else
      avg_time[i] = 0;
  }

  /* print summary */
  fprintf(stdout, "\n--- Thread Load Balancing Statistics ---\n");
  fprintf(fp_out, "\n--- Thread Load Balancing Statistics ---\n");

  fprintf(stdout, "Threads: %ld\n\n", opt_threads);
  fprintf(fp_out, "Threads: %ld\n\n", opt_threads);

  /* overall imbalance */
  double overall_avg = overall_total_time / opt_threads;
  double overall_imbalance = (overall_avg > 0) ? overall_max_time / overall_avg : 1.0;
  double overall_efficiency = (overall_max_time > 0) ?
                              (overall_total_time / (opt_threads * overall_max_time)) * 100.0 : 100.0;

  fprintf(stdout, "Overall Thread Time (seconds):\n");
  fprintf(stdout, "  Min: %.3f  Max: %.3f  Avg: %.3f\n",
          overall_min_time / 1e6, overall_max_time / 1e6, overall_avg / 1e6);
  fprintf(stdout, "  Imbalance ratio: %.3f (1.0 = perfect)\n", overall_imbalance);
  fprintf(stdout, "  Parallel efficiency: %.1f%%\n\n", overall_efficiency);

  fprintf(fp_out, "Overall Thread Time (seconds):\n");
  fprintf(fp_out, "  Min: %.3f  Max: %.3f  Avg: %.3f\n",
          overall_min_time / 1e6, overall_max_time / 1e6, overall_avg / 1e6);
  fprintf(fp_out, "  Imbalance ratio: %.3f (1.0 = perfect)\n", overall_imbalance);
  fprintf(fp_out, "  Parallel efficiency: %.1f%%\n\n", overall_efficiency);

  /* per-work-type summary */
  fprintf(stdout, "Per-proposal-type timing (seconds):\n");
  fprintf(stdout, "  %-8s %10s %10s %10s %10s %10s\n",
          "Type", "Total", "Min", "Max", "Avg", "Imbal");
  fprintf(fp_out, "Per-proposal-type timing (seconds):\n");
  fprintf(fp_out, "  %-8s %10s %10s %10s %10s %10s\n",
          "Type", "Total", "Min", "Max", "Avg", "Imbal");

  for (i = 1; i < THREAD_WORK_COUNT; ++i)
  {
    if (total_calls[i] > 0)
    {
      double imbalance = (avg_time[i] > 0) ? max_time[i] / avg_time[i] : 1.0;
      fprintf(stdout, "  %-8s %10.3f %10.3f %10.3f %10.3f %10.3f\n",
              work_type_name(i),
              total_time[i] / 1e6,
              min_time[i] / 1e6,
              max_time[i] / 1e6,
              avg_time[i] / 1e6,
              imbalance);
      fprintf(fp_out, "  %-8s %10.3f %10.3f %10.3f %10.3f %10.3f\n",
              work_type_name(i),
              total_time[i] / 1e6,
              min_time[i] / 1e6,
              max_time[i] / 1e6,
              avg_time[i] / 1e6,
              imbalance);
    }
  }

  /* detailed per-thread breakdown */
  if (detailed)
  {
    fprintf(stdout, "\nPer-thread breakdown (seconds):\n");
    fprintf(fp_out, "\nPer-thread breakdown (seconds):\n");

    /* header */
    fprintf(stdout, "  Thread");
    fprintf(fp_out, "  Thread");
    for (i = 1; i < THREAD_WORK_COUNT; ++i)
    {
      if (total_calls[i] > 0)
      {
        fprintf(stdout, " %8s", work_type_name(i));
        fprintf(fp_out, " %8s", work_type_name(i));
      }
    }
    fprintf(stdout, "    Total\n");
    fprintf(fp_out, "    Total\n");

    /* per-thread data */
    for (t = 0; t < opt_threads; ++t)
    {
      thread_info_t * tip = ti + t;
      double thread_total = 0;

      fprintf(stdout, "  %6ld", t);
      fprintf(fp_out, "  %6ld", t);

      for (i = 1; i < THREAD_WORK_COUNT; ++i)
      {
        if (total_calls[i] > 0)
        {
          fprintf(stdout, " %8.3f", tip->total_time_us[i] / 1e6);
          fprintf(fp_out, " %8.3f", tip->total_time_us[i] / 1e6);
          thread_total += tip->total_time_us[i];
        }
      }
      fprintf(stdout, " %8.3f\n", thread_total / 1e6);
      fprintf(fp_out, " %8.3f\n", thread_total / 1e6);
    }
  }

  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");
}
