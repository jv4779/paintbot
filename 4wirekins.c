/********************************************************************
* Description:               _
*   Kinematics for 4-4 wire RPR parallel robot with 2-DOF
*
*
* Author: 
* License: GPL Version 2
* System: Linux
*    
* Copyright (c) 2004 All rights reserved.
*
* Last change:
********************************************************************/

/*
-- A2 *__                                ______* A4
 ^       \__L2                    __L4__/
 |           \___    Eb    ______/
 |               *--------*
Ly             Eh|        |
 |            ___*--------*______
 |        ___/                   \______
 v     __/ L1                       L3  \______
-- A1 *                                        * A3
      |<-------------- Lx -------------------->|

A1 = {0, 0}
A2 = {0, Ly}
A3 = {Lx, 0}
A4 = {Lx, Ly}

 */

#include "kinematics.h"             /* these decls */

/* ident tag */
#ifndef __GNUC__
#ifndef __attribute__
#define __attribute__(x)
#endif
#endif

#ifdef MAIN
#define hal_float_t double
#define rtapi_print printf
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#else
#include "rtapi_math.h"
#include "hal.h"
#endif
struct haldata {
    hal_float_t *lx, *ly, *eb, *eh;
} *haldata = 0;

#define Lx (*(haldata->lx))
#define Ly (*(haldata->ly))
#define Eb (*(haldata->eb))
#define Eh (*(haldata->eh))

#define sq(x) ((x)*(x))

int circle_circle_intersection(double x0, double y0, double r0,
                               double x1, double y1, double r1,
                               double *xi, double *yi,
                               double *xi_prime, double *yi_prime);

/*
  forward kinematics takes three strut lengths and computes Dx, Dy, and Dz
  pos->tran.x,y,z, respectively. The forward flag is used to resolve
  D above/below the xy plane. The inverse flags are not set since there
  are no ambiguities going from world to joint coordinates.

*/
int kinematicsForward(const double * joints,
                      EmcPose * pos,
                      const KINEMATICS_FORWARD_FLAGS * fflags,
                      KINEMATICS_INVERSE_FLAGS * iflags)
{
#define L1 (joints[0])
#define L2 (joints[1])
#define L3 (joints[2])
#define L4 (joints[3])
#define Dx (pos->tran.x)
#define Dy (pos->tran.y)
  double xi_1_3,yi_1_3,xi_2_4,yi_2_4;
  double halfEb = Eb/2.0, halfEh = Eh/2.0;

  //rtapi_print("kinematicsForward 1=%f 2=%f 3=%f 4=%f\n",L1,L2,L3,L4);

  if ( !circle_circle_intersection(
		halfEb, halfEh, L1,
        	Lx - halfEb, halfEh, L3,
        	&xi_1_3, &yi_1_3, 0, 0) )
  {
    return -1;
  }
  if ( !circle_circle_intersection(
		halfEb, Ly - halfEh, L2,
        	Lx - halfEb, Ly - halfEh, L4,
        	0, 0, &xi_2_4, &yi_2_4) )
  {
    return -1;
  }

  Dx = (xi_1_3 + xi_2_4)/2.0;
  Dy = (yi_1_3 + yi_2_4)/2.0;

  //rtapi_print("kinematicsForward result x=%f y=%f\n",Dx,Dy);

  pos->tran.z = 0.0;
  pos->a = 0.0;
  pos->b = 0.0;
  pos->c = 0.0;

  return 0;

#undef L1
#undef L2
#undef L3
#undef L4
#undef Dx
#undef Dy
}

/*
L1 = Sqrt[Abs[0.5 b - x]^2 + Abs[0.5 h - y]^2]

L2 = Sqrt[Abs[0.5 b - x]^2 + Abs[-0.5 h + Ly - y]^2]

L3 = Sqrt[Abs[-0.5 b + Lx - x]^2 + Abs[0.5 h - y]^2]

L4 = Sqrt[Abs[-0.5 b + Lx - x]^2 + Abs[-0.5 h + Ly - y]^2]
*/

int kinematicsInverse(const EmcPose * pos,
                      double * joints,
                      const KINEMATICS_INVERSE_FLAGS * iflags,
                      KINEMATICS_FORWARD_FLAGS * fflags)
{
#define L1 (joints[0])
#define L2 (joints[1])
#define L3 (joints[2])
#define L4 (joints[3])
#define Dx (pos->tran.x)
#define Dy (pos->tran.y)

  double halfEb = Eb/2.0, halfEh = Eh/2.0;

  //rtapi_print("kinematicsInverse x=%f y=%f\n",Dx,Dy);

  if ( Dx < halfEb || Dy < halfEh || Dx > Lx - halfEb || Dy > Ly - halfEh ) {
    return -1;
  }

  L1 = sqrt(sq(fabs(halfEb - Dx)) + sq(fabs(halfEh - Dy)));

  L2 = sqrt(sq(fabs(halfEb - Dx)) + sq(fabs(-halfEh + Ly - Dy)));

  L3 = sqrt(sq(fabs(-halfEb + Lx - Dx)) + sq(fabs(halfEh - Dy)));

  L4 = sqrt(sq(fabs(-halfEb + Lx - Dx)) + sq(fabs(-halfEh + Ly - Dy)));

  //rtapi_print("kinematicsInverse result 1=%f 2=%f 3=%f 4=%f\n",L1,L2,L3,L4);

  return 0;

#undef AD
#undef BD
#undef CD
#undef DD
#undef Dx
#undef Dy
}

KINEMATICS_TYPE kinematicsType()
{
  return KINEMATICS_BOTH;
}

#ifdef MAIN

#include <stdio.h>
#include <string.h>

/*
  Interactive testing of kins.
  Compile: gcc -g -DRTAPI -DMAIN -I/usr/include/linuxcnc 4wirekins.c -lm
  Syntax: a.out <Lx> <Ly> <Eb> <Eh>
*/
int main(int argc, char *argv[])
{
#ifndef BUFFERLEN
#define BUFFERLEN 256
#endif
  char buffer[BUFFERLEN];
  char cmd[BUFFERLEN];
  EmcPose pos, vel;
  double joints[3], jointvels[3];
  char inverse;
  char flags;
  KINEMATICS_FORWARD_FLAGS fflags;

  haldata = malloc(sizeof(haldata));
  haldata->lx = malloc(sizeof(double));
  haldata->ly = malloc(sizeof(double));
  haldata->eb = malloc(sizeof(double));
  haldata->eh = malloc(sizeof(double));

  inverse = 0;			/* forwards, by default */
  fflags = 0;			/* above xy plane, by default */
  if (argc != 5 ||
      1 != sscanf(argv[1], "%lf", &Lx) ||
      1 != sscanf(argv[2], "%lf", &Ly) ||
      1 != sscanf(argv[3], "%lf", &Eb) ||
      1 != sscanf(argv[4], "%lf", &Eh)) {
    fprintf(stderr, "syntax: %s Lx Ly Eb Eh\n", argv[0]);
    return 1;
  }

  while (! feof(stdin)) {
    if (inverse) {
	printf("inv> ");
    }
    else {
	printf("fwd> ");
    }
    fflush(stdout);

    if (NULL == fgets(buffer, BUFFERLEN, stdin)) {
      break;
    }
    if (1 != sscanf(buffer, "%s", cmd)) {
      continue;
    }

    if (! strcmp(cmd, "quit")) {
      break;
    }
    if (! strcmp(cmd, "i")) {
      inverse = 1;
      continue;
    }
    if (! strcmp(cmd, "f")) {
      inverse = 0;
      continue;
    }

    if (inverse) {		/* inverse kins */
      if (2 != sscanf(buffer, "%lf %lf", 
		      &pos.tran.x,
		      &pos.tran.y)) {
	printf("need X Y\n");
	continue;
      }
      if (0 != kinematicsInverse(&pos, joints, NULL, &fflags)) {
	printf("inverse kin error\n");
      }
      else {
	printf("%f\t%f\t%f\t%f\n", joints[0], joints[1], joints[2], joints[3]);
	if (0 != kinematicsForward(joints, &pos, &fflags, NULL)) {
	  printf("forward kin error\n");
	}
	else {
	  printf("%f\t%f\n", pos.tran.x, pos.tran.y);
	}
      }
    }
    else {			/* forward kins */
      if (4 != sscanf(buffer, "%lf %lf %lf %lf", 
			&joints[0],
			&joints[1],
			&joints[2],
			&joints[3])) {
	printf("need joints 0 1 2 3\n");
	continue;
      }
      if (0 != kinematicsForward(joints, &pos, &fflags, NULL)) {
	printf("forward kin error\n");
      }
      else {
	printf("%f\t%f\n", pos.tran.x, pos.tran.y);
	if (0 != kinematicsInverse(&pos, joints, NULL, &fflags)) {
	  printf("inverse kin error\n");
	}
	else {
	  printf("%f\t%f\t%f\t%f\n", joints[0], joints[1], joints[2], joints[3]);
	}
      }
    }
  } /* end while (! feof(stdin)) */

  return 0;
}

#else

#include "rtapi.h"		/* RTAPI realtime OS API */
#include "rtapi_app.h"		/* RTAPI realtime module decls */
#include "hal.h"

EXPORT_SYMBOL(kinematicsType);
EXPORT_SYMBOL(kinematicsForward);
EXPORT_SYMBOL(kinematicsInverse);

MODULE_LICENSE("GPL");



int comp_id;
int rtapi_app_main(void) {
    int res = 0;

    comp_id = hal_init("4wirekins");
    if(comp_id < 0) return comp_id;

    haldata = hal_malloc(sizeof(struct haldata));
    if(!haldata) goto error;

    if((res = hal_pin_float_new("4wirekins.Lx", HAL_IO, &(haldata->lx), comp_id)) < 0) goto error;
    if((res = hal_pin_float_new("4wirekins.Ly", HAL_IO, &(haldata->ly), comp_id)) < 0) goto error;
    if((res = hal_pin_float_new("4wirekins.Eb", HAL_IO, &(haldata->eb), comp_id)) < 0) goto error;
    if((res = hal_pin_float_new("4wirekins.Eh", HAL_IO, &(haldata->eh), comp_id)) < 0) goto error;

    Lx = Ly = 1.0;
    Eb = Eh = 0.0;
    hal_ready(comp_id);
    return 0;

error:
    hal_exit(comp_id);
    return res;
}

void rtapi_app_exit(void) { hal_exit(comp_id); }

#endif /* MAIN */

/* circle_circle_intersection() *
 * Determine the points where 2 circles in a common plane intersect.
 *
 * int circle_circle_intersection(
 *                                // center and radius of 1st circle
 *                                double x0, double y0, double r0,
 *                                // center and radius of 2nd circle
 *                                double x1, double y1, double r1,
 *                                // 1st intersection point
 *                                double *xi, double *yi,              
 *                                // 2nd intersection point
 *                                double *xi_prime, double *yi_prime)
 *
 * This is a public domain work. 3/26/2005 Tim Voght
 *
 */

int circle_circle_intersection(double x0, double y0, double r0,
                               double x1, double y1, double r1,
                               double *xi, double *yi,
                               double *xi_prime, double *yi_prime)
{
  double a, dx, dy, d, h, rx, ry;
  double x2, y2;

  /* dx and dy are the vertical and horizontal distances between
   * the circle centers.
   */
  dx = x1 - x0;
  dy = y1 - y0;

  /* Determine the straight-line distance between the centers. */
  d = sqrt((dy*dy) + (dx*dx));

  /* Check for solvability. */
  if (d > (r0 + r1))
  {
    /* no solution. circles do not intersect. */
    return 0;
  }
  if (d < fabs(r0 - r1))
  {
    /* no solution. one circle is contained in the other */
    return 0;
  }

  /* 'point 2' is the point where the line through the circle
   * intersection points crosses the line between the circle
   * centers.  
   */

  /* Determine the distance from point 0 to point 2. */
  a = ((r0*r0) - (r1*r1) + (d*d)) / (2.0 * d) ;

  /* Determine the coordinates of point 2. */
  x2 = x0 + (dx * a/d);
  y2 = y0 + (dy * a/d);

  /* Determine the distance from point 2 to either of the
   * intersection points.
   */
  h = sqrt((r0*r0) - (a*a));

  /* Now determine the offsets of the intersection points from
   * point 2.
   */
  rx = -dy * (h/d);
  ry = dx * (h/d);

  /* Determine the absolute intersection points. */
  if ( xi )
    *xi = x2 + rx;
  if ( xi_prime )
    *xi_prime = x2 - rx;
  if ( yi )
    *yi = y2 + ry;
  if ( yi_prime )
    *yi_prime = y2 - ry;

  return 1;
}

