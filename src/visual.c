/*
    Copyright (C) 2016-2024 Tomas Flouri, Bruce Rannala and Ziheng Yang

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
#include "pdfgen.h"

typedef struct coord_s
{
  double x;
  double y;
} coord_t;

//static double doc_height = PDF_MM_TO_POINT(150.0f);
static double doc_height = PDF_A4_WIDTH;
//static double doc_width = PDF_MM_TO_POINT(198.0f);
static double doc_width = PDF_A4_HEIGHT;
static double margin_top = 40;
static double margin_bottom = 60;
static double margin_left = 50;
static double margin_right = 90;
static double padding_top = 20;
static double padding_bottom = 20;
static double padding_left = 0;
static double padding_right = 0;


static coord_t * coord = NULL;

#if 0
static int cb_trav_full(snode_t * x)
{
  if (!x->left)
    return 0;

  return 1;
}
#endif

static void pdf_add_vert_arrow_dn(struct pdf_doc *pdf, struct pdf_object *page,
                                  float x, float y, float width, float height,
                                  float border_width, uint32_t colour)
{
    double ah_ratio = 0.1;  /* arrow head ratio */
    double sh_ratio = 0.4;  /* shaft width ratio */
    float arrow_x[] = {x + width*(sh_ratio / 2),        /* 1 */
                       x + width*(sh_ratio / 2),        /* 2 */
                       x + width / 2,                   /* 3 */
                       x,                               /* 4 */
                       x - width / 2,                   /* 5 */
                       x - width*(sh_ratio / 2),        /* 6 */
                       x - width*(sh_ratio / 2),        /* 7 */
                       x + width*(sh_ratio / 2)};       /* 8 */

    float arrow_y[] = {y,                               /* 1 */
                       y - height*(1-ah_ratio),         /* 2 */
                       y - height*(1-ah_ratio),         /* 3 */
                       y - height,                      /* 4 */
                       y - height*(1-ah_ratio),         /* 5 */
                       y - height*(1-ah_ratio),         /* 6 */
                       y,                               /* 7 */
                       y};                              /* 8 */
    pdf_add_polygon(pdf, page, arrow_x, arrow_y,
                    sizeof(arrow_x) / sizeof(arrow_x[0]), border_width,
                    colour);
}
static void pdf_add_vert_arrow_up(struct pdf_doc *pdf, struct pdf_object *page,
                                  float x, float y, float width, float height,
                                  float border_width, uint32_t colour)
{
    double ah_ratio = 0.1;  /* arrow head ratio */
    double sh_ratio = 0.4;  /* shaft width ratio */
    float arrow_x[] = {x + width*(sh_ratio / 2),        /* 1 */
                       x + width*(sh_ratio / 2),        /* 2 */
                       x + width / 2,                   /* 3 */
                       x,                               /* 4 */
                       x - width / 2,                   /* 5 */
                       x - width*(sh_ratio / 2),        /* 6 */
                       x - width*(sh_ratio / 2),        /* 7 */
                       x + width*(sh_ratio / 2)};       /* 8 */

    float arrow_y[] = {y,                               /* 1 */
                       y + height*(1-ah_ratio),         /* 2 */
                       y + height*(1-ah_ratio),         /* 3 */
                       y + height,                      /* 4 */
                       y + height*(1-ah_ratio),         /* 5 */
                       y + height*(1-ah_ratio),         /* 6 */
                       y,                               /* 7 */
                       y};                              /* 8 */
    pdf_add_polygon(pdf, page, arrow_x, arrow_y,
                    sizeof(arrow_x) / sizeof(arrow_x[0]), border_width,
                    colour);
}


static int cb_trav_full(snode_t * x)
{
  if (!x->left && !x->right && !x->hybrid)
    return 0;

  if (x->hybrid && x->hybrid->mark[0])
    return 0;

  return 1;
}

static void set_solid_branch(const stree_t * stree, snode_t * x)
{
  snode_t * h;
  const long thread_index_zero = 0;

  assert(x->hybrid);

  if (node_is_bidirection(x))
  {
    /* set hybrid node branch as solid */
    h = (x->left && x->right) ? x : x->hybrid;
    assert(!h->hybrid->left && !h->hybrid->right);

    h->mark[0] = 1;
    h->hybrid->mark[0] = 0;

    return;
  }

  /* already solved ? */
  if ((x->mark[0] && !x->hybrid->mark[0]) || (!x->mark[0] && x->hybrid->mark[0]))
    return;

  /* check if parallel */
  if (x->parent == x->hybrid->parent)
  {
    h = (x->left || x->right) ? x : x->hybrid;
    assert(h->left);
    assert(!h->hybrid->left && !h->hybrid->right);

    h->mark[0] = 1;
    h->hybrid->mark[0] = 0;

    return;
  }

  /* all other cases */

  /* if the parent of the hybridization node is hybrid, sort it out first */
  if (x->parent->hybrid)
    set_solid_branch(stree, x->parent);
  if (x->hybrid->parent->hybrid)
    set_solid_branch(stree, x->hybrid->parent);


  snode_t * p1 = x->parent;
  snode_t * p2 = x->hybrid->parent;

  while (p1->hybrid)
  {
    assert((p1->mark[0] && !p1->hybrid->mark[0]) || (!p1->mark[0] && p1->hybrid->mark[0]));
    p1 = p1->mark[0] ? p1->parent : p1->hybrid->parent;
  }

  while (p2->hybrid)
  {
    assert((p2->mark[0] && !p2->hybrid->mark[0]) || (!p2->mark[0] && p2->hybrid->mark[0]));
    p2 = p2->mark[0] ? p2->parent : p2->hybrid->parent;
  }

  assert(p1 && p2);
  if (stree->pptable[p1->node_index][p2->node_index] == 1)
  {
    x->mark[0] = 0;
    x->hybrid->mark[0] = 1;
  }
  else if (stree->pptable[p2->node_index][p1->node_index] == 1)
  {
    x->mark[0] = 1;
    x->hybrid->mark[0] = 0;
  }
  else
  {
    if (x->hphi > x->hybrid->hphi)
    {
      x->mark[0] = 1;
      x->hybrid->mark[0] = 0;
    }
    else if (x->hphi < x->hybrid->hphi)
    {
      x->mark[0] = 0;
      x->hybrid->mark[0] = 1;
    }
    else
    {
      /* randomly pick one */
      if (legacy_rndu(thread_index_zero) < 0.5)
      {
        x->mark[0] = 1;
        x->hybrid->mark[0] = 0;
      }
      else
      {
        x->mark[0] = 0;
        x->hybrid->mark[0] = 1;
      }
    }
  }
}

static void network_traverse_postorder(snode_t * node,
                                       int (*cbtrav)(snode_t *),
                                       unsigned int * index,
                                       snode_t ** outbuffer)
{
  printf("Visiting node: %s (%d)\n", node->label, node->node_index);
  if (node->hybrid)
    printf("    hybrid node: htau: %ld hyb->htau: %ld\n", node->htau, node->hybrid->htau);
  if (!node->left && !node->right && !node->hybrid)
  {
    if (cbtrav(node))
    {
      outbuffer[*index] = node;
      (*index)++;
    }
    printf("Skipped1 node %s %d(parent: %s %d)\n", node->label, node->node_index, node->parent->label, node->parent->node_index);
    return;
  }
  if (!cbtrav(node))
  {
    printf("Skipped2 node %s\n", node->label);
    return;
  }

  if (node->hybrid && !node->htau && node->hybrid->htau)
  {
    printf("Skipping node %s (%d)\n", node->label, node->node_index);
    return;
  }

  if (node->hybrid)
    node->mark[0] = node->hybrid->mark[0] = 1;

  if (node->hybrid)
  {
    if (node->left) // && (!node->left->hybrid || !node->left->mark))
      network_traverse_postorder(node->left, cbtrav, index, outbuffer);
    else if (node->hybrid->left)
      network_traverse_postorder(node->hybrid->left, cbtrav, index, outbuffer);
  }
  else
  {
    if (node->left) // && (!node->left->hybrid || !node->left->mark))
      network_traverse_postorder(node->left, cbtrav, index, outbuffer);
  }

  if (node->hybrid)
  {
    if (node->right) // && (!node->right->hybrid || !node->right->mark))
      network_traverse_postorder(node->right, cbtrav, index, outbuffer);
    else if (node->hybrid->right)
      network_traverse_postorder(node->hybrid->right, cbtrav, index, outbuffer);
  }
  else
  {
    if (node->right) // && (!node->left->hybrid || !node->left->mark))
      network_traverse_postorder(node->right, cbtrav, index, outbuffer);
  }

  outbuffer[*index] = node;
  (*index)++;
  if (node->hybrid)
    printf(ANSI_COLOR_RED
           "Stored node: %s %d [%ld (%s) %ld (%s)]\n" ANSI_COLOR_RESET,
           node->label, node->node_index, node->htau, node->parent->label,
           node->hybrid->htau, node->hybrid->parent->label);
  else
    printf(ANSI_COLOR_RED "Stored node: %s\n" ANSI_COLOR_RESET, node->label);
}
static void plot_tree_with_grid(struct pdf_doc * pdf,
                                const stree_t * stree,
                                double linewidth)
{
  unsigned int i,node_count = 0;
  unsigned int li,ri,xi;
  unsigned int tip_count = 0;

  long xtics = 6;
  double xtic_xshift = 3;
  double xtic_yshift = -3;
  double xtic_fontsize = 8;
  double xtic_angle = -5*BPP_PI/12;
  double ci_xshift = 0;
  double ci_yshift = 5;
  double ci_lw = 1;
  double fontsize_label = 8;
  double arrow_offset = 0.05;
  double arrow_width = 8;

  double canvas_width = doc_width - margin_left - margin_right;
  double canvas_height = doc_height - margin_bottom - margin_top;
  double treecanvas_height = canvas_height - (padding_top+padding_bottom);
  double treecanvas_width = canvas_width - (padding_left+padding_right);

  nodepinfo_t * rootni = (nodepinfo_t *)(stree->root->data);
  double maxage_ext = 0.02;
  double maxage = MAX(rootni->age,rootni->hi)*(1+maxage_ext);
  double xscaler = treecanvas_width / maxage;
  double tipspace = treecanvas_height / (stree->tip_count-1);
  double cap_size = 6;
  double timeline_lw = linewidth/2.;

  char * label = NULL;
  snode_t * left;
  snode_t * right;
  snode_t * x;

  /* TODO: Rescale timeline if left or right padding */
  assert(padding_left == 0);
  assert(padding_right == 0);

  /* draw background */
  pdf_add_filled_rectangle(pdf,
                           NULL,
                           margin_left,
                           margin_bottom,
                           canvas_width,
                           canvas_height,
                           0,
                           PDF_RGB(0xeb,0xeb,0xeb),
                           PDF_BLACK);

  #if 0
  /* draw canvas border */
  pdf_add_rectangle(pdf,
                    NULL,
                    margin_left,
                    margin_bottom,
                    canvas_width,
                    canvas_height,
                    linewidth,
                    PDF_BLACK);
  #endif

  /* draw timeline */
  pdf_add_line(pdf,
               NULL,
               margin_left,
               margin_bottom,
               margin_left+canvas_width,
               margin_bottom,
               timeline_lw,
               PDF_BLACK);

  /* draw left and right cap as independent line segments until I learn how
     to draw a line with a cap in postscript (I expect to die first) */

  /* draw left cap */
  pdf_add_line(pdf,
               NULL,
               margin_left,
               margin_bottom-cap_size/2.,
               margin_left,
               margin_bottom+cap_size/2.,
               timeline_lw,
               PDF_BLACK);

  /* draw right cap */
  pdf_add_line(pdf,
               NULL,
               margin_left+canvas_width,
               margin_bottom-cap_size/2.,
               margin_left+canvas_width,
               margin_bottom+cap_size/2.,
               timeline_lw,
               PDF_BLACK);

  /* show number of xtics */
  double xtic_space = maxage / (xtics+1);
  for (i = 0; i < xtics; ++i)
  {
    /* convert x position canvas coordinate */
    double x = xtic_space*(i+1)*xscaler;

    /* draw xtic marker */
    pdf_add_line(pdf,
                 NULL,
                 margin_left+canvas_width-x,
                 margin_bottom-cap_size/2.,
                 margin_left+canvas_width-x,
                 margin_bottom+cap_size/2.,
                 timeline_lw,
                 PDF_BLACK);

    /* draw xtic value */
    xasprintf(&label,"%.6f",xtic_space*(i+1));
    pdf_add_text_rotate(pdf,
                        NULL,
                        label,
                        xtic_fontsize,
                        margin_left+canvas_width-x+xtic_xshift,
                        margin_bottom+xtic_yshift,
                        xtic_angle,
                        PDF_BLACK);
    free(label);
  }

  /* draw major xtic gridlines */
  for (i = 0; i < xtics; ++i)
  {
    /* convert x position canvas coordinate */
    double x = xtic_space*(i+1)*xscaler;

    pdf_add_line(pdf,
                 NULL,
                 margin_left+canvas_width-x,
                 margin_bottom+cap_size/2.,
                 margin_left+canvas_width-x,
                 margin_bottom+canvas_height,
                 timeline_lw,
                 PDF_RGB(0xf6,0xf6,0xf6));
  }

  /* get nodes in postorder */
  stree_traverse(stree->root,
                 TREE_TRAVERSE_POSTORDER,
                 cb_trav_full,
                 stree->td,
                 &node_count);

  char trunclabel[256];
  size_t truncsize = 7;
  for (i  = 0; i < node_count; ++i)
  {
    /* get trio of nodes and indices */
    x = stree->td[i];
    left = x->left;
    right = x->right;

    li = left->node_index;
    ri = right->node_index;
    xi = x->node_index;

    /* if left or right child is a tip, calculate its coordinates */
    if (li < stree->tip_count)
    {
      if (strlen(stree->nodes[li]->label) > 7)
      {
        memcpy(trunclabel,stree->nodes[li]->label,truncsize);
        trunclabel[truncsize+0] = '.';
        trunclabel[truncsize+1] = '.';
        trunclabel[truncsize+2] = '.';
        trunclabel[truncsize+3] = '\0';
      }
      else
      {
        strcpy(trunclabel,stree->nodes[li]->label);
      }

      coord[li].x = doc_width - margin_right;
      coord[li].y = margin_bottom+padding_bottom + tipspace * tip_count++;

      pdf_add_text_rotate(pdf, NULL, trunclabel, fontsize_label,
                          coord[li].x+5, coord[li].y-5, 0, PDF_BLACK);
    }
    if (ri < stree->tip_count)
    {
      if (strlen(stree->nodes[ri]->label) > 7)
      {
        memcpy(trunclabel,stree->nodes[ri]->label,truncsize);
        trunclabel[truncsize+0] = '.';
        trunclabel[truncsize+1] = '.';
        trunclabel[truncsize+2] = '.';
        trunclabel[truncsize+3] = '\0';
      }
      else
      {
        strcpy(trunclabel,stree->nodes[ri]->label);
      }

      coord[ri].x = doc_width - margin_right;
      coord[ri].y = margin_bottom+padding_bottom + tipspace * tip_count++;

      pdf_add_text_rotate(pdf, NULL, trunclabel, fontsize_label,
                          coord[ri].x+5, coord[ri].y-5, 0, PDF_BLACK);
    }

    /* calculate xy coordinates for node x */
    coord[xi].x = margin_left + canvas_width - x->tau * xscaler;
    if (coord[ri].y > coord[li].y)
      coord[xi].y = coord[li].y + (coord[ri].y - coord[li].y) / 2.0f;
    else
      coord[xi].y = coord[ri].y + (coord[li].y - coord[ri].y) / 2.0f;

    /* draw the three line segments around an inner node */
    pdf_add_line(pdf, NULL, coord[li].x, coord[li].y,
                 coord[xi].x, coord[li].y, linewidth, PDF_BLACK);        
    pdf_add_line(pdf, NULL, coord[ri].x, coord[ri].y,
                 coord[xi].x, coord[ri].y, linewidth, PDF_BLACK);        
    pdf_add_line(pdf, NULL, coord[xi].x, coord[li].y,
                 coord[xi].x, coord[ri].y, linewidth, PDF_BLACK);        

    /* print inner node label (node index) */
    xasprintf(&label, "%d", x->node_index+1);
    if (x->parent)
      pdf_add_text(pdf, NULL, label, fontsize_label,
                   coord[xi].x+5, coord[xi].y-5, PDF_BLACK);
    else
      pdf_add_text(pdf, NULL, label, fontsize_label,
                   coord[xi].x-13, coord[xi].y-5, PDF_BLACK);
    free(label);
  }

  /* draw 95% HPD interval */
  for (i  = 0; i < node_count; ++i)
  {
    /* get trio of nodes and indices */
    x = stree->td[i];
    xi = x->node_index;

    /* draw 95% HPD interval */
    nodepinfo_t * ni = (nodepinfo_t *)(x->data);
    double lox = margin_left + canvas_width - ni->lo*xscaler;
    double hix = margin_left + canvas_width - ni->hi*xscaler;
    pdf_add_line(pdf,
                 NULL,
                 lox+ci_xshift,
                 coord[xi].y+ci_yshift,
                 hix+ci_xshift,
                 coord[xi].y+ci_yshift,
                 ci_lw,
                 PDF_BLUE);
    pdf_add_line(pdf,
                 NULL,
                 lox+ci_xshift,
                 coord[xi].y+ci_yshift+cap_size/2.,
                 lox+ci_xshift,
                 coord[xi].y+ci_yshift+cap_size/2.,
                 ci_lw,
                 PDF_BLUE);

  }

  /* print max value on timeline */
  xasprintf(&label, "%.6f", maxage );
  pdf_add_text_rotate(pdf,
                      NULL,
                      label,
                      xtic_fontsize,
                      margin_left+xtic_xshift,
                      margin_bottom+xtic_yshift,
                      xtic_angle,
                      PDF_BLACK);
  free(label);


  /* print 0 (present) on timeline */
  xasprintf(&label, "0");
  pdf_add_text_rotate(pdf,
                      NULL,
                      label,
                      xtic_fontsize,
                      margin_left+canvas_width+xtic_xshift,
                      margin_bottom+xtic_yshift,
                      0,
                      PDF_BLACK);
  free(label);

  /* show all migration events */
  for (i = 0; i < opt_migration_count; ++i)
  {
    migspec_t * spec = opt_mig_specs+i;

    snode_t * src = stree->nodes[spec->si];
    snode_t * dst = stree->nodes[spec->ti];

    assert(src->parent);
    assert(dst->parent);

    /* get migration window */
    double s1t = coord[src->node_index].x;
    double s1p = coord[src->parent->node_index].x;
    double d1t = coord[dst->node_index].x;
    double d1p = coord[dst->parent->node_index].x;

    double lo = MIN(s1t,d1t);
    double hi = MAX(s1p,d1p);

    double ax = hi+(lo-hi)/2;

    double ay = coord[src->node_index].y;
    double aheight = fabs(ay-coord[dst->node_index].y);
    #if 0
    fprintf(stdout, "%s->%s (y:%f y:%f)  x: %f y: %f awith: %f\n",
            src->label, dst->label, coord[src->node_index].y,
            coord[dst->node_index].y, ax, ay, awidth);
    #endif
    if (coord[src->node_index].y > coord[dst->node_index].y)
      pdf_add_vert_arrow_dn(pdf,
                            NULL,
                            ax,
                            ay-arrow_offset*aheight,
                            arrow_width,
                            aheight*(1-2*arrow_offset),
                            1,
                            PDF_RGB(0xff, 0x80, 0x40));
    else
      pdf_add_vert_arrow_up(pdf,
                            NULL,
                            ax,
                            ay+arrow_offset*aheight,
                            arrow_width,
                            aheight*(1-2*arrow_offset),
                            1,
                            PDF_RGB(0xff, 0x80, 0x40));
  }
}

static void plot_tree_with_timeline(struct pdf_doc * pdf,
                                    const stree_t * stree,
                                    double linewidth,
                                    double fontsize,
                                    double timeline_height)
{
  unsigned int i,node_count = 0;
  unsigned int li,ri,xi;
  unsigned int tip_count = 0;

  long xtics = 5;
  double xtic_xshift = 3;
  double xtic_yshift = -3;
  double xtic_fontsize = 8;
  double xtic_angle = -5*BPP_PI/12;
  double ci_xshift = 0;
  double ci_yshift = 5;
  double ci_lw = 2;
  double tl_projection_width = 1;
  double timeline_y = margin_bottom+timeline_height/4;

  double canvas_width = doc_width - margin_left - margin_right;
  double canvas_height = doc_height - margin_bottom - margin_top;
  double treecanvas_height = canvas_height - timeline_height;

  nodepinfo_t * rootni = (nodepinfo_t *)(stree->root->data);
  double maxage = MAX(rootni->age,rootni->hi);
  double xscaler = canvas_width / maxage;
  double tipspace = treecanvas_height / (stree->tip_count-1);
  double cap_size = 6;
  double timeline_lw = linewidth/2.;

  char * label = NULL;
  snode_t * left;
  snode_t * right;
  snode_t * x;
  #if 0
  /* draw canvas border */
  pdf_add_rectangle(pdf,
                    NULL,
                    margin_left,
                    margin_bottom,
                    canvas_width,
                    canvas_height,
                    linewidth,
                    PDF_BLACK);
  #endif

  /* draw timeline */
  pdf_add_line(pdf,
               NULL,
               margin_left,
               timeline_y,
               margin_left+canvas_width,
               timeline_y,
               timeline_lw,
               PDF_GREEN);

  /* draw left and right cap as independent line segments until I learn how
     to draw a line with a cap in postscript (I expect to die first) */

  /* draw left cap */
  pdf_add_line(pdf,
               NULL,
               margin_left,
               timeline_y-cap_size/2.,
               margin_left,
               timeline_y+cap_size/2.,
               timeline_lw,
               PDF_GREEN);

  /* draw right cap */
  pdf_add_line(pdf,
               NULL,
               margin_left+canvas_width,
               timeline_y-cap_size/2.,
               margin_left+canvas_width,
               timeline_y+cap_size/2.,
               timeline_lw,
               PDF_GREEN);

  /* get nodes in postorder */
  stree_traverse(stree->root,
                 TREE_TRAVERSE_POSTORDER,
                 cb_trav_full,
                 stree->td,
                 &node_count);

  for (i  = 0; i < node_count; ++i)
  {
    /* get trio of nodes and indices */
    x = stree->td[i];
    left = x->left;
    right = x->right;

    li = left->node_index;
    ri = right->node_index;
    xi = x->node_index;

    /* if left or right child is a tip, calculate its coordinates */
    if (li < stree->tip_count)
    {
      coord[li].x = doc_width - margin_right;
      coord[li].y = margin_bottom+timeline_height + tipspace * tip_count++;

      //pdf_add_text(pdf, NULL, stree->nodes[li]->label, fontsize,
      //             coord[li].x+5, coord[li].y-5, PDF_BLACK);
      pdf_add_text_rotate(pdf, NULL, stree->nodes[li]->label, fontsize,
                   coord[li].x+5, coord[li].y-5, 0, PDF_BLACK);
    }
    if (ri < stree->tip_count)
    {
      coord[ri].x = doc_width - margin_right;
      coord[ri].y = margin_bottom+timeline_height + tipspace * tip_count++;

      //pdf_add_text(pdf, NULL, stree->nodes[ri]->label, fontsize,
      //             coord[ri].x+5, coord[ri].y-5, PDF_BLACK);
      pdf_add_text_rotate(pdf, NULL, stree->nodes[ri]->label, fontsize,
                          coord[ri].x+5, coord[ri].y-5, 0, PDF_BLACK);
    }

    /* calculate xy coordinates for node x */
    coord[xi].x = margin_left + canvas_width - x->tau * xscaler;
    if (coord[ri].y > coord[li].y)
      coord[xi].y = coord[li].y + (coord[ri].y - coord[li].y) / 2.0f;
    else
      coord[xi].y = coord[ri].y + (coord[li].y + coord[ri].y) / 2.0f;

    /* draw the three line segments around an inner node */
    pdf_add_line(pdf, NULL, coord[li].x, coord[li].y,
                 coord[xi].x, coord[li].y, linewidth, PDF_BLACK);        
    pdf_add_line(pdf, NULL, coord[ri].x, coord[ri].y,
                 coord[xi].x, coord[ri].y, linewidth, PDF_BLACK);        
    pdf_add_line(pdf, NULL, coord[xi].x, coord[li].y,
                 coord[xi].x, coord[ri].y, linewidth, PDF_BLACK);        

    /* print inner node label (node index) */
    xasprintf(&label, "%d", x->node_index+1);
    if (x->parent)
      pdf_add_text(pdf, NULL, label, fontsize,
                   coord[xi].x+5, coord[xi].y-5, PDF_BLACK);
    else
      pdf_add_text(pdf, NULL, label, fontsize,
                   coord[xi].x-13, coord[xi].y-5, PDF_BLACK);
    free(label);

    /* add line from current node to timeline */
    pdf_add_line(pdf,
                 NULL,
                 coord[xi].x,
                 coord[xi].y,
                 coord[xi].x,
                 timeline_y,
                 tl_projection_width,
                 PDF_BLUE);

    /* draw 95% HPD interval */
    nodepinfo_t * ni = (nodepinfo_t *)(x->data);
    double lox = margin_left + canvas_width - ni->lo*xscaler;
    double hix = margin_left + canvas_width - ni->hi*xscaler;
    pdf_add_line(pdf,
                 NULL,
                 lox+ci_xshift,
                 coord[xi].y+ci_yshift,
                 hix+ci_xshift,
                 coord[xi].y+ci_yshift,
                 ci_lw,
                 PDF_BLUE);
    pdf_add_line(pdf,
                 NULL,
                 lox+ci_xshift,
                 coord[xi].y+ci_yshift+cap_size/2.,
                 lox+ci_xshift,
                 coord[xi].y+ci_yshift+cap_size/2.,
                 ci_lw,
                 PDF_BLUE);
  }

  /* print max value on timeline */
  xasprintf(&label, "%.6f", maxage );
  pdf_add_text_rotate(pdf,
                      NULL,
                      label,
                      xtic_fontsize,
                      margin_left+xtic_xshift,
                      timeline_y+xtic_yshift,
                      xtic_angle,
                      PDF_BLACK);
  free(label);


  /* print 0 (present) on timeline */
  xasprintf(&label, "0");
  pdf_add_text_rotate(pdf,
                      NULL,
                      label,
                      xtic_fontsize,
                      margin_left+canvas_width+xtic_xshift,
                      timeline_y+xtic_yshift,
                      0,
                      PDF_BLACK);
  free(label);

  /* show all migration events */
  for (i = 0; i < opt_migration_count; ++i)
  {
    migspec_t * spec = opt_mig_specs+i;

    snode_t * src = stree->nodes[spec->si];
    snode_t * dst = stree->nodes[spec->ti];
    
    assert(src->parent);
    assert(dst->parent);
  }



  /* show number of xtics */
  double xtic_space = maxage / (xtics+1);
  for (i = 0; i < xtics; ++i)
  {
    /* convert x position canvas coordinate */
    double x = xtic_space*(i+1)*xscaler;

    /* draw xtic marker */
    pdf_add_line(pdf,
                 NULL,
                 margin_left+canvas_width-x,
                 timeline_y-cap_size/2.,
                 margin_left+canvas_width-x,
                 timeline_y+cap_size/2.,
                 timeline_lw,
                 PDF_GREEN);

    /* draw xtic value */
    xasprintf(&label,"%.6f",xtic_space*(i+1));
    pdf_add_text_rotate(pdf,
                        NULL,
                        label,
                        xtic_fontsize,
                        margin_left+canvas_width-x+xtic_xshift,
                        timeline_y+xtic_yshift,
                        xtic_angle,
                        PDF_BLACK);
    free(label);
  }
}

static void plot_tree(struct pdf_doc * pdf,
                      const stree_t * stree,
                      double linewidth,
                      double fontsize)
{
  unsigned int i,node_count = 0;
  unsigned int li,ri,xi;
  unsigned int tip_count = 0;
  snode_t * left;
  snode_t * right;
  snode_t * x;

  double canvas_width = doc_width - margin_left - margin_right;
  double canvas_height = doc_height - margin_bottom - margin_top;
  double xscaler = canvas_width / stree->root->tau;
  double tipspace = canvas_height / (stree->tip_count-1);

  stree_traverse(stree->root,
                 TREE_TRAVERSE_POSTORDER,
                 cb_trav_full,
                 stree->td,
                 &node_count);

  for (i  = 0; i < node_count; ++i)
  {
    x = stree->td[i];
    left = x->left;
    right = x->right;

    li = left->node_index;
    ri = right->node_index;
    xi = x->node_index;

    /* if left or right child is a tip, calculate its coordinates */
    if (li < stree->tip_count)
    {
      coord[li].x = doc_width - margin_right;
      coord[li].y = margin_bottom + tipspace * tip_count++;

      pdf_add_text(pdf, NULL, stree->nodes[li]->label, fontsize,
                   coord[li].x+5, coord[li].y-5, PDF_BLACK);
    }
    if (ri < stree->tip_count)
    {
      coord[ri].x = doc_width - margin_right;
      coord[ri].y = margin_bottom + tipspace * tip_count++;

      pdf_add_text(pdf, NULL, stree->nodes[ri]->label, fontsize,
                   coord[ri].x+5, coord[ri].y-5, PDF_BLACK);
    }

    /* calculate coordinates of node x */
    coord[xi].x = margin_left + canvas_width - x->tau * xscaler;
    if (coord[ri].y > coord[li].y)
      coord[xi].y = coord[li].y + (coord[ri].y - coord[li].y) / 2.0f;
    else
      coord[xi].y = coord[ri].y + (coord[li].y + coord[ri].y) / 2.0f;

    pdf_add_line(pdf, NULL, coord[li].x, coord[li].y,
                 coord[xi].x, coord[li].y, linewidth, PDF_BLACK);        
    pdf_add_line(pdf, NULL, coord[ri].x, coord[ri].y,
                 coord[xi].x, coord[ri].y, linewidth, PDF_BLACK);        
    pdf_add_line(pdf, NULL, coord[xi].x, coord[li].y,
                 coord[xi].x, coord[ri].y, linewidth, PDF_BLACK);        

    char * nodelabel = NULL;
    xasprintf(&nodelabel, "%d", x->node_index+1);
    if (x->parent)
      pdf_add_text(pdf, NULL, nodelabel, fontsize,
                   coord[xi].x+5, coord[xi].y-5, PDF_BLACK);
    else
      pdf_add_text(pdf, NULL, nodelabel, fontsize,
                   coord[xi].x-13, coord[xi].y-5, PDF_BLACK);
    free(nodelabel);
  }
}

static void plot_network(struct pdf_doc * pdf, const stree_t * stree)
{
  unsigned int i,node_count = 0;
  unsigned int li,ri,xi;
  unsigned int tip_count = 0;

  long xtics = 6;
  double xtic_xshift = 3;
  double xtic_yshift = -3;
  double xtic_fontsize = 8;
  double xtic_angle = -5*BPP_PI/12;
  double ci_xshift = 0;
  double ci_yshift = 5;
  double ci_lw = 1;
  double fontsize_label = 8;
  //double arrow_offset = 0.05;
  //double arrow_width = 8;
  double linewidth = 3;

  double canvas_width = doc_width - margin_left - margin_right;
  double canvas_height = doc_height - margin_bottom - margin_top;
  double treecanvas_height = canvas_height - (padding_top+padding_bottom);
  double treecanvas_width = canvas_width - (padding_left+padding_right);

  nodepinfo_t * rootni = (nodepinfo_t *)(stree->root->data);
  double maxage_ext = 0.02;
  double maxage = MAX(rootni->age,rootni->hi)*(1+maxage_ext);
  double xscaler = treecanvas_width / maxage;
  double tipspace = treecanvas_height / (stree->tip_count-1);
  double cap_size = 6;
  double timeline_lw = linewidth/2.;

  char * label = NULL;
  snode_t * left;
  snode_t * right;
  snode_t * x;

  double fontsize = 8;

  /* draw background */
  pdf_add_filled_rectangle(pdf,
                           NULL,
                           margin_left,
                           margin_bottom,
                           canvas_width,
                           canvas_height,
                           0,
                           PDF_RGB(0xeb,0xeb,0xeb),
                           PDF_BLACK);

  network_traverse_postorder(stree->root,
                             cb_trav_full,
                             &node_count,
                             stree->td);
  printf("Nodes in buffer: %d\n", node_count);
  for (i = 0; i < node_count; ++i)
  {
    snode_t * x = stree->td[i];
    if (x->hybrid)
      printf("%s (%d) [%ld (%s) %ld (%s)]\n", x->label, x->node_index, x->htau, x->parent->label, x->hybrid->htau, x->hybrid->parent->label);
    else
      printf("%s (%d)\n", x->label, x->node_index);
  }
  
  for (i = 0; i < stree->hybrid_count; ++i)
    set_solid_branch(stree,stree->nodes[stree->tip_count+stree->inner_count+i]);

  for (i = 0; i < node_count; ++i)
  {
    x = stree->td[i];

    if (x->hybrid && node_is_bidirection(x))
      continue;

    left = x->left;
    right = x->right;

    if (x->hybrid && !x->left && x->hybrid->left)
      left = x->hybrid->left;


    /* if left or right child is a tip, calculate its coordinates */
    if (left && left->node_index < stree->tip_count)
    {
      printf("Processing left tip %s\n", left->label);
      li = left->node_index;
      coord[li].x = doc_width - margin_right;
      coord[li].y = margin_bottom + tipspace * tip_count++;

      pdf_add_text(pdf, NULL, stree->nodes[li]->label, fontsize,
                   coord[li].x+5, coord[li].y-5, PDF_BLACK);
    }
    if (right && right->node_index < stree->tip_count)
    {
      printf("Processing right tip %s\n", right->label);
      ri = right->node_index;
      coord[ri].x = doc_width - margin_right;
      coord[ri].y = margin_bottom + tipspace * tip_count++;

      pdf_add_text(pdf, NULL, stree->nodes[ri]->label, fontsize,
                   coord[ri].x+5, coord[ri].y-5, PDF_BLACK);
    }

    /* calculate coordinates of node x */
    if (x->left && x->right && x->left->hybrid && !x->left->mark[0] && x->right->hybrid && !x->right->mark[0])
    {
      assert(0);
      /* TODO: Create a point between the middle of the two hybrid nodes */
    }
    else if (x->left && x->right && x->right->hybrid && !x->right->mark[0])
    {
      printf("++ Processing node %s\n", x->label);
      li = left->node_index;
      xi = x->node_index;

      coord[xi].x = margin_left + canvas_width - x->tau * xscaler;
      coord[xi].y = coord[li].y;

      pdf_add_line(pdf, NULL, coord[li].x, coord[li].y,
                   coord[xi].x, coord[xi].y, linewidth, PDF_BLACK);        

      char * nodelabel = NULL;
      xasprintf(&nodelabel, "%d", x->node_index+1);
      if (x->parent)
        pdf_add_text(pdf, NULL, nodelabel, fontsize,
                     coord[xi].x+5, coord[xi].y-5, PDF_BLACK);
      else
        pdf_add_text(pdf, NULL, nodelabel, fontsize,
                     coord[xi].x-13, coord[xi].y-5, PDF_BLACK);
      free(nodelabel);
    }
    else if (x->left && x->right && x->left->hybrid && !x->left->mark[0])
    {
      printf("-- Processing node %s\n", x->label);
      ri = right->node_index;
      xi = x->node_index;

      coord[xi].x = margin_left + canvas_width - x->tau * xscaler;
      coord[xi].y = coord[ri].y;

      pdf_add_line(pdf, NULL, coord[ri].x, coord[ri].y,
                   coord[xi].x, coord[xi].y, linewidth, PDF_BLACK);        

      char * nodelabel = NULL;
      xasprintf(&nodelabel, "%d", x->node_index+1);
      if (x->parent)
        pdf_add_text(pdf, NULL, nodelabel, fontsize,
                     coord[xi].x+5, coord[xi].y-5, PDF_BLACK);
      else
        pdf_add_text(pdf, NULL, nodelabel, fontsize,
                     coord[xi].x-13, coord[xi].y-5, PDF_BLACK);
      free(nodelabel);
    }
    else if ((x->left && x->right) && 
             ((!x->left->hybrid && !x->right->hybrid) ||
              (x->left->hybrid && x->left->mark[0]) ||
              (x->right->hybrid && x->right->mark[0])))
    {
      li = left->node_index;
      ri = right->node_index;
      xi = x->node_index;

      coord[xi].x = margin_left + canvas_width - x->tau * xscaler;
      if (coord[ri].y > coord[li].y)
        coord[xi].y = coord[li].y + (coord[ri].y - coord[li].y) / 2.0f;
      else
        coord[xi].y = coord[ri].y + (coord[li].y + coord[ri].y) / 2.0f;

      printf("Processing node: %s [ %s (%d)  %s (%d) ]     :   (%f,%f)  [ (%f,%f)   (%f,%f)  ]\n", x->label, x->left->label, x->left->node_index+1, x->right->label, x->right->node_index+1, coord[xi].x, coord[xi].y, coord[li].x, coord[li].y, coord[ri].x, coord[ri].y);
      pdf_add_line(pdf, NULL, coord[li].x, coord[li].y,
                   coord[xi].x, coord[li].y, linewidth, PDF_BLACK);        
      pdf_add_line(pdf, NULL, coord[ri].x, coord[ri].y,
                   coord[xi].x, coord[ri].y, linewidth, PDF_BLACK);        
      pdf_add_line(pdf, NULL, coord[xi].x, coord[li].y,
                   coord[xi].x, coord[ri].y, linewidth, PDF_BLACK);        

      char * nodelabel = NULL;
      xasprintf(&nodelabel, "%d", x->node_index+1);
      if (x->parent)
        pdf_add_text(pdf, NULL, nodelabel, fontsize,
                     coord[xi].x+5, coord[xi].y-5, PDF_BLACK);
      else
        pdf_add_text(pdf, NULL, nodelabel, fontsize,
                     coord[xi].x-13, coord[xi].y-5, PDF_BLACK);
      free(nodelabel);
    }
    else if (x->hybrid)
    {
      /* TODO: Workout bidirection */
      assert(!node_is_bidirection(x));
      if (!x->left && !x->right)
        x = x->hybrid;
      assert(x->left && !x->right);

      left = x->left;

      li = left->node_index;
      xi = x->node_index;

      coord[xi].x = margin_left + canvas_width - x->tau * xscaler;
      coord[xi].y = coord[li].y;
      /* also set hybrid */
      coord[x->hybrid->node_index].x = coord[xi].x;
      coord[x->hybrid->node_index].y = coord[xi].y;

      printf("Processing hybrid: %s [ %ld (%s) %ld (%s) ]\n", x->label, x->htau, x->parent->label, x->hybrid->htau, x->hybrid->parent->label);
      if (x->left)
        printf("     Left: %s   [ %f   %f]\n", x->left->label, coord[x->left->node_index].x,coord[x->left->node_index].y);
      if (x->right)
        printf("     Right: %s   [ %f   %f]\n", x->right->label, coord[x->right->node_index].x,coord[x->right->node_index].y);

      pdf_add_line(pdf, NULL, coord[li].x, coord[li].y,
                   coord[xi].x, coord[xi].y, linewidth, PDF_BLACK);        
      /* TF: 23.10.2024 */ 
      //pdf_add_vert_arrow_up(pdf,
      //                      NULL,
      //                      coord[li].x,
      //                      coord[li].y,
      //                      5, //coord[xi].x,
      //                      coord[xi].y-coord[li].y,
      //                      1,
      //                      PDF_RGB(0xff, 0x80, 0x40));

      char * nodelabel= NULL;
      xasprintf(&nodelabel, "%d", x->node_index+1);
      if (x->parent)
        pdf_add_text(pdf, NULL, nodelabel, fontsize,
                     coord[xi].x+5, coord[xi].y-5, PDF_RED);
      else
        pdf_add_text(pdf, NULL, nodelabel, fontsize,
                     coord[xi].x-13, coord[xi].y-5, PDF_RED);
      free(nodelabel);
    }
    else
    {
      printf("x node: %s\n", x->label);
      assert(x->left && x->right);
      if (x->left->hybrid)
        printf("left->hybrid (%s) mark: %d\n", x->left->label, x->left->mark[0]);
      if (x->right->hybrid)
        printf("right->hybrid (%s) mark: %d\n", x->right->label, x->right->mark[0]);
      assert(x->hybrid);
      assert(0);
      /* many cases here */
    }
  }
  for (i = 0; i < node_count; ++i)
  {
    snode_t * x = stree->td[i];
    if (!x->hybrid) continue;

    x = x->hybrid;
    snode_t * p = x->parent;
    unsigned int pi = p->node_index;
    xi = x->node_index;

    pdf_add_line(pdf, NULL, coord[pi].x, coord[pi].y,
                 coord[xi].x, coord[xi].y, linewidth, PDF_BLACK);
  }
}

void stree_export_pdf(const stree_t * stree)
{
  unsigned int i,total_nodes;
  char * outfile = NULL;
  struct pdf_doc * pdf;

  struct pdf_info info = {
      .creator = "bpp" "4.7.1",
      .producer = "bpp" "4.7.1",
      .title = "Binary species tree",
      .author = "Tomas Flouri",
      .subject = "",
      .date = "Today"
  };

  total_nodes = stree->tip_count + stree->inner_count + stree->hybrid_count;

  coord = (coord_t *)xmalloc((size_t)total_nodes * sizeof(coord_t));

  pdf = pdf_create(doc_width, doc_height, &info);
  pdf_append_page(pdf);
  //pdf_set_font(pdf, "Times-Roman");
  pdf_set_font(pdf, "Helvetica");

  for (i = 0; i < total_nodes; ++i)
    stree->nodes[i]->mark[0] = 0;

  if (!opt_msci)
  {
    plot_tree_with_grid(pdf,stree,3);
  }
  else
  {
    plot_network(pdf,stree);
    //exit(0);
  }

  for (i = 0; i < total_nodes; ++i)
    stree->nodes[i]->mark[0] = 0;
  
  /* time signature */
  struct tm * lt = NULL;
  char buffer[256];
  time_t t = time(NULL);
  lt = localtime(&t);
  assert(lt);

  /* strftime(buffer, 256, "%a %b %d %T %Y", lt); */
  strftime(buffer, 256, "%c", lt);

  float buffer_size = 0;
  pdf_get_font_text_width(pdf,"Courier-Bold", buffer, 8, &buffer_size);

  pdf_set_font(pdf, "Courier-Bold");
  pdf_add_text(pdf, NULL, buffer, 8,
               5, 5, PDF_BLACK);
  pdf_add_text(pdf, NULL, cmdline, 8,
               5+buffer_size+5, 5, PDF_BLACK);
  pdf_add_text(pdf,
               NULL, 
               "Created with: bpp " PROG_VERSION,
               8,
               5, doc_height-15, PDF_BLACK);
  
  pdf_get_font_text_width(pdf, "Courier-Bold", PVER_SHA1, 8, &buffer_size);
  pdf_add_text(pdf,
               NULL, 
               PVER_SHA1,
               8,
               doc_width-buffer_size-10, doc_height-15, PDF_BLACK);

  xasprintf(&outfile, "%s.pdf", opt_jobname);
  pdf_save(pdf,outfile);
  free(outfile);
  pdf_destroy(pdf);

  free(coord);
}

