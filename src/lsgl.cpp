/* Routines for linear multiple output using sparse group lasso regression.
 Intended for use with R.
 Copyright (C) 2014 Martin Vincent

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>
 */

//Uncomment to turn on debuging
//#undef NDEBUG

//Configuration
//Debugging
#ifndef NDEBUG
#define SGL_DEBUG
#endif

//Runtime checking for numerical problems
#define SGL_RUNTIME_CHECKS

//Check dimension of input objects
#define SGL_DIM_CHECKS

//Converges checks
#define SGL_CONVERGENCE_CHECK

//Exception handling
#define SGL_CATCH_EXCEPTIONS

// print information abt convergence
//#define SGL_DEBUG_INFO_QUADRATIC

//Should the timers be activated (only needed for profiling the code)
//#define DO_TIMING

//Show entering and leving selected functions
//#define FUNC_ENTER


//Sgl optimizer
#include <sgl.h>
#include "pkg_c_config.h"

// Objectives
#include "frobenius_norm.h"
#include "frobenius_norm_weighted.h"

/**********************************
 *
 *  lsgl X dense Y dense module
 *
 *********************************/

// Module name
#define MODULE_NAME lsgl_xd_yd

#define OBJECTIVE frobenius

#include <sgl/RInterface/sgl_lambda_seq.h>
#include <sgl/RInterface/sgl_fit.h>

#define PREDICTOR sgl::LinearPredictor < sgl::matrix , sgl::LinearResponse >

#include <sgl/RInterface/sgl_predict.h>
#include <sgl/RInterface/sgl_subsampling.h>

/**********************************
 *
 *  lsgl weighted X dense Y dense module
 *
 *********************************/

// Reset macros
#undef MODULE_NAME
#undef OBJECTIVE
#undef PREDICTOR

// Module name
#define MODULE_NAME lsgl_w_xd_yd

#define OBJECTIVE frobenius_w

#include <sgl/RInterface/sgl_lambda_seq.h>
#include <sgl/RInterface/sgl_fit.h>

#define PREDICTOR sgl::LinearPredictor < sgl::matrix , sgl::LinearResponse >

#include <sgl/RInterface/sgl_predict.h>
#include <sgl/RInterface/sgl_subsampling.h>

/*********************************
 *
 *  lsgl X sparse Y dense module
 *
 *********************************/
// Reset macros
#undef MODULE_NAME
#undef OBJECTIVE
#undef PREDICTOR

// Module name
#define MODULE_NAME lsgl_xs_yd

#define OBJECTIVE frobenius_spx

#include <sgl/RInterface/sgl_lambda_seq.h>
#include <sgl/RInterface/sgl_fit.h>

#define PREDICTOR sgl::LinearPredictor < sgl::sparse_matrix , sgl::LinearResponse >

#include <sgl/RInterface/sgl_predict.h>
#include <sgl/RInterface/sgl_subsampling.h>

/*********************************
 *
 *  lsgl X dense Y sparse module
 *
 *********************************/
// Reset macros
#undef MODULE_NAME
#undef OBJECTIVE
#undef PREDICTOR

// Module name
#define MODULE_NAME lsgl_xd_ys

#define OBJECTIVE frobenius_spy

#include <sgl/RInterface/sgl_lambda_seq.h>
#include <sgl/RInterface/sgl_fit.h>

#define PREDICTOR sgl::LinearPredictor < sgl::matrix , sgl::LinearResponse >

#include <sgl/RInterface/sgl_predict.h>
#include <sgl/RInterface/sgl_subsampling.h>

/*********************************
 *
 *  lsgl X dense Y sparse module
 *
 *********************************/
// Reset macros
#undef MODULE_NAME
#undef OBJECTIVE
#undef PREDICTOR

// Module name
#define MODULE_NAME lsgl_xs_ys

#define OBJECTIVE frobenius_spx_spy

#include <sgl/RInterface/sgl_lambda_seq.h>
#include <sgl/RInterface/sgl_fit.h>

#define PREDICTOR sgl::LinearPredictor < sgl::sparse_matrix , sgl::LinearResponse >

#include <sgl/RInterface/sgl_predict.h>
#include <sgl/RInterface/sgl_subsampling.h>


/* **********************************
 *
 *  Registration of methods
 *
 ***********************************/

#include <R_ext/Rdynload.h>

static const R_CallMethodDef sglCallMethods[] = {

  SGL_LAMBDA(lsgl_xd_yd), SGL_LAMBDA(lsgl_xs_yd),
  SGL_LAMBDA(lsgl_xd_ys), SGL_LAMBDA(lsgl_xs_ys),
  SGL_LAMBDA(lsgl_w_xd_yd),

  SGL_FIT(lsgl_xd_yd), SGL_FIT(lsgl_xs_yd),
  SGL_FIT(lsgl_xd_ys), SGL_FIT(lsgl_xs_ys),
  SGL_FIT(lsgl_w_xd_yd),

  SGL_PREDICT(lsgl_xd_yd), SGL_PREDICT(lsgl_xs_yd),
  SGL_PREDICT(lsgl_xd_ys), SGL_PREDICT(lsgl_xs_ys),

  SGL_SUBSAMPLING(lsgl_xd_yd), SGL_SUBSAMPLING(lsgl_xs_yd),
  SGL_SUBSAMPLING(lsgl_xd_ys), SGL_SUBSAMPLING(lsgl_xs_ys),
  SGL_SUBSAMPLING(lsgl_w_xd_yd),

  {"r_pkg_c_config", (DL_FUNC) &r_pkg_c_config, 0},

  {NULL, NULL, 0}

};

extern "C" {
  void R_init_lsgl(DllInfo *info);
}

void R_init_lsgl(DllInfo *info)
{
  // Register the .Call routines.
  R_registerRoutines(info, NULL, sglCallMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
