/*
 * EssexAidProcessVFrameSwitchable.c
 *
 * Embedded MATLAB Coder code generation for function 'EssexAidProcessVFrameSwitchable'
 *
 * C source code generated on: Thu Jul 14 16:03:12 2011
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "EssexAidProcessVFrameSwitchable.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void m_abs(real_T eml_x_data[6912], int32_T eml_x_sizes[2], real_T
 eml_y_data[6912], int32_T eml_y_sizes[2]);
static void m_b_eml_scalexp_alloc(real_T eml_varargin_1_data[6912], int32_T
 eml_varargin_1_sizes[2], real_T eml_z_data[6912], int32_T eml_z_sizes[2]);
static void m_b_filter(real_T eml_b[2], real_T eml_a[3], real_T
 eml_x_data[6912], int32_T eml_x_sizes[2], real_T eml_zi_data[190], int32_T
 eml_zi_sizes[1], real_T eml_y_data[6912], int32_T eml_y_sizes[2], real_T
 eml_zf[2]);
static void m_b_max(real_T eml_varargin_1_data[6912], int32_T
 eml_varargin_1_sizes[2], real_T eml_maxval_data[6912], int32_T
 eml_maxval_sizes[2]);
static void m_b_power(real_T eml_a_data[6912], int32_T eml_a_sizes[2], real_T
 eml_b, real_T eml_y_data[6912], int32_T eml_y_sizes[2]);
static void m_c_eml_scalexp_alloc(real_T eml_varargin_1_data[6912], int32_T
 eml_varargin_1_sizes[2], real_T eml_z_data[6912], int32_T eml_z_sizes[2]);
static void m_c_filter(real_T eml_b[5], real_T eml_a[5], real_T
 eml_x_data[6912], int32_T eml_x_sizes[2], real_T eml_zi_data[3000], int32_T
 eml_zi_sizes[1], real_T eml_y_data[6912], int32_T eml_y_sizes[2], real_T
 eml_zf[4]);
static void m_c_power(real_T eml_b_data[6912], int32_T eml_b_sizes[2], real_T
 eml_y_data[6912], int32_T eml_y_sizes[2]);
static int32_T m_compute_nones(boolean_T eml_x_data[6912], int32_T
 eml_x_sizes[2], int32_T eml_n);
static void m_d_eml_scalexp_alloc(real_T eml_varargin_2_data[6912], int32_T
 eml_varargin_2_sizes[2], real_T eml_z_data[6912], int32_T eml_z_sizes[2]);
static void m_eml_li_find(boolean_T eml_x_data[6912], int32_T eml_x_sizes[2],
 int32_T eml_y_data[6912], int32_T eml_y_sizes[2]);
static void m_eml_scalexp_alloc(real_T eml_varargin_1_data[6912], int32_T
 eml_varargin_1_sizes[2], real_T eml_z_data[6912], int32_T eml_z_sizes[2]);
static void m_filter(real_T eml_b[2], real_T eml_a[2], real_T eml_x_data[6912],
 int32_T eml_x_sizes[2], real_T eml_zi, real_T eml_y_data[6912], int32_T
 eml_y_sizes[2], real_T *eml_zf);
static void m_find(boolean_T eml_x_data[6912], int32_T eml_x_sizes[2], real_T
 eml_i_data[6912], int32_T eml_i_sizes[2]);
static boolean_T m_isfinite(real_T eml_x);
static real_T m_max(real_T eml_varargin_1_data[6912], int32_T
 eml_varargin_1_sizes[2]);
static void m_power(real_T eml_a_data[6912], int32_T eml_a_sizes[2], real_T
 eml_y_data[6912], int32_T eml_y_sizes[2]);
static void m_refp1_DRNL_brokenstick_nl(real_T eml_x_data[6912], int32_T
 eml_x_sizes[2], real_T eml_b, real_T eml_c);
static void m_refp1_log10(real_T eml_x_data[6912], int32_T eml_x_sizes[2]);
static void m_refp1_sign(real_T eml_x_data[6912], int32_T eml_x_sizes[2]);
static void m_refp1_sqrt(real_T eml_x_data[6912], int32_T eml_x_sizes[2]);
static void m_sum(real_T eml_x_data[13824], int32_T eml_x_sizes[2], real_T
 eml_y_data[6912], int32_T eml_y_sizes[2]);

/* Function Definitions */
static void m_abs(real_T eml_x_data[6912], int32_T eml_x_sizes[2], real_T
 eml_y_data[6912], int32_T eml_y_sizes[2])
{
  int32_T eml_loop_ub;
  int16_T eml_iv9[2];
  int16_T eml_iv10[2];
  int16_T eml_tmp_sizes[2];
  int32_T eml_k;
  for(eml_loop_ub = 0; eml_loop_ub < 2; eml_loop_ub++) {
    eml_iv9[eml_loop_ub] = (int16_T)eml_x_sizes[eml_loop_ub];
  }
  eml_iv10[0] = 1;
  eml_iv10[1] = eml_iv9[1];
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_iv10[1];
  for(eml_loop_ub = 0; eml_loop_ub < 2; eml_loop_ub++) {
    eml_y_sizes[eml_loop_ub] = eml_tmp_sizes[eml_loop_ub];
  }
  eml_loop_ub = eml_x_sizes[1];
  for(eml_k = 1; eml_k <= eml_loop_ub; eml_k++) {
    eml_y_data[eml_k - 1] = fabs(eml_x_data[eml_k - 1]);
  }
}
static void m_b_eml_scalexp_alloc(real_T eml_varargin_1_data[6912], int32_T
 eml_varargin_1_sizes[2], real_T eml_z_data[6912],
 int32_T eml_z_sizes[2])
{
  int32_T eml_i1;
  int16_T eml_iv12[2];
  int16_T eml_iv13[2];
  int16_T eml_tmp_sizes[2];
  if(!(eml_varargin_1_sizes[1] == 1)) {
    for(eml_i1 = 0; eml_i1 < 2; eml_i1++) {
      eml_iv12[eml_i1] = (int16_T)eml_varargin_1_sizes[eml_i1];
    }
    eml_iv13[0] = 1;
    eml_iv13[1] = eml_iv12[1];
    eml_tmp_sizes[0] = 1;
    eml_tmp_sizes[1] = eml_iv13[1];
    for(eml_i1 = 0; eml_i1 < 2; eml_i1++) {
      eml_z_sizes[eml_i1] = eml_tmp_sizes[eml_i1];
    }
  } else {
    for(eml_i1 = 0; eml_i1 < 2; eml_i1++) {
      eml_z_sizes[eml_i1] = 1;
    }
  }
}
static void m_b_filter(real_T eml_b[2], real_T eml_a[3], real_T
 eml_x_data[6912], int32_T eml_x_sizes[2], real_T eml_zi_data[190],
 int32_T eml_zi_sizes[1], real_T eml_y_data[6912], int32_T eml_y_sizes[2],
 real_T eml_zf[2])
{
  real_T eml_a1;
  int32_T eml_k;
  int16_T eml_iv5[2];
  int16_T eml_iv6[2];
  int16_T eml_tmp_sizes[2];
  int32_T eml_nx;
  real_T eml_dbuffer[3];
  int32_T eml_j;
  eml_a1 = eml_a[0];
  if((!m_isfinite(eml_a1)) || (eml_a1 == 0.0) || (!(eml_a1 != 1.0))) {
  } else {
    for(eml_k = 0; eml_k < 2; eml_k++) {
      eml_b[eml_k] /= eml_a1;
    }
    for(eml_k = 2; eml_k < 4; eml_k++) {
      eml_a[eml_k - 1] /= eml_a1;
    }
    eml_a[0] = 1.0;
  }
  for(eml_k = 0; eml_k < 2; eml_k++) {
    eml_iv5[eml_k] = (int16_T)eml_x_sizes[eml_k];
  }
  eml_iv6[0] = 1;
  eml_iv6[1] = eml_iv5[1];
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_iv6[1];
  eml_nx = eml_x_sizes[1];
  for(eml_k = 0; eml_k < 2; eml_k++) {
    eml_y_sizes[eml_k] = eml_tmp_sizes[eml_k];
    eml_dbuffer[eml_k + 1] = eml_zi_data[eml_k];
  }
  for(eml_j = 1; eml_j <= eml_nx; eml_j++) {
    for(eml_k = 0; eml_k < 2; eml_k++) {
      eml_dbuffer[eml_k] = eml_dbuffer[eml_k + 1];
    }
    eml_dbuffer[2] = 0.0;
    for(eml_k = 0; eml_k < 2; eml_k++) {
      eml_dbuffer[eml_k] += eml_x_data[eml_j - 1] * eml_b[eml_k];
    }
    for(eml_k = 2; eml_k < 4; eml_k++) {
      eml_dbuffer[eml_k - 1] -= eml_dbuffer[0] * eml_a[eml_k - 1];
    }
    eml_y_data[eml_j - 1] = eml_dbuffer[0];
  }
  for(eml_k = 0; eml_k < 2; eml_k++) {
    eml_zf[eml_k] = eml_dbuffer[eml_k + 1];
  }
}
static void m_b_max(real_T eml_varargin_1_data[6912], int32_T
 eml_varargin_1_sizes[2], real_T eml_maxval_data[6912], int32_T
 eml_maxval_sizes[2])
{
  int32_T eml_loop_ub;
  int32_T eml_k;
  m_c_eml_scalexp_alloc(eml_varargin_1_data, eml_varargin_1_sizes,
   eml_maxval_data, eml_maxval_sizes);
  eml_loop_ub = eml_maxval_sizes[1];
  for(eml_k = 1; eml_k <= eml_loop_ub; eml_k++) {
    eml_maxval_data[eml_k - 1] = rt_MAXd_snf(eml_varargin_1_data[eml_k - 1],
     0.0);
  }
}
static void m_b_power(real_T eml_a_data[6912], int32_T eml_a_sizes[2], real_T
 eml_b, real_T eml_y_data[6912], int32_T eml_y_sizes[2]
 )
{
  int32_T eml_loop_ub;
  int32_T eml_k;
  m_b_eml_scalexp_alloc(eml_a_data, eml_a_sizes, eml_y_data, eml_y_sizes);
  eml_loop_ub = eml_y_sizes[1];
  for(eml_k = 1; eml_k <= eml_loop_ub; eml_k++) {
    eml_y_data[eml_k - 1] = rt_pow_snf(eml_a_data[eml_k - 1], eml_b);
  }
}
static void m_c_eml_scalexp_alloc(real_T eml_varargin_1_data[6912], int32_T
 eml_varargin_1_sizes[2], real_T eml_z_data[6912],
 int32_T eml_z_sizes[2])
{
  int32_T eml_i2;
  int16_T eml_iv15[2];
  int16_T eml_iv16[2];
  int16_T eml_tmp_sizes[2];
  if(!(eml_varargin_1_sizes[1] == 1)) {
    for(eml_i2 = 0; eml_i2 < 2; eml_i2++) {
      eml_iv15[eml_i2] = (int16_T)eml_varargin_1_sizes[eml_i2];
    }
    eml_iv16[0] = 1;
    eml_iv16[1] = eml_iv15[1];
    eml_tmp_sizes[0] = 1;
    eml_tmp_sizes[1] = eml_iv16[1];
    for(eml_i2 = 0; eml_i2 < 2; eml_i2++) {
      eml_z_sizes[eml_i2] = eml_tmp_sizes[eml_i2];
    }
  } else {
    for(eml_i2 = 0; eml_i2 < 2; eml_i2++) {
      eml_z_sizes[eml_i2] = 1;
    }
  }
}
static void m_c_filter(real_T eml_b[5], real_T eml_a[5], real_T
 eml_x_data[6912], int32_T eml_x_sizes[2], real_T eml_zi_data[3000],
 int32_T eml_zi_sizes[1], real_T eml_y_data[6912], int32_T eml_y_sizes[2],
 real_T eml_zf[4])
{
  real_T eml_a1;
  int32_T eml_k;
  int16_T eml_iv7[2];
  int16_T eml_iv8[2];
  int16_T eml_tmp_sizes[2];
  int32_T eml_nx;
  real_T eml_dbuffer[5];
  int32_T eml_j;
  eml_a1 = eml_a[0];
  if((!m_isfinite(eml_a1)) || (eml_a1 == 0.0) || (!(eml_a1 != 1.0))) {
  } else {
    for(eml_k = 0; eml_k < 5; eml_k++) {
      eml_b[eml_k] /= eml_a1;
    }
    for(eml_k = 2; eml_k < 6; eml_k++) {
      eml_a[eml_k - 1] /= eml_a1;
    }
    eml_a[0] = 1.0;
  }
  for(eml_k = 0; eml_k < 2; eml_k++) {
    eml_iv7[eml_k] = (int16_T)eml_x_sizes[eml_k];
  }
  eml_iv8[0] = 1;
  eml_iv8[1] = eml_iv7[1];
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_iv8[1];
  for(eml_k = 0; eml_k < 2; eml_k++) {
    eml_y_sizes[eml_k] = eml_tmp_sizes[eml_k];
  }
  eml_nx = eml_x_sizes[1];
  for(eml_k = 0; eml_k < 4; eml_k++) {
    eml_dbuffer[eml_k + 1] = eml_zi_data[eml_k];
  }
  for(eml_j = 1; eml_j <= eml_nx; eml_j++) {
    for(eml_k = 0; eml_k < 4; eml_k++) {
      eml_dbuffer[eml_k] = eml_dbuffer[eml_k + 1];
    }
    eml_dbuffer[4] = 0.0;
    for(eml_k = 0; eml_k < 5; eml_k++) {
      eml_dbuffer[eml_k] += eml_x_data[eml_j - 1] * eml_b[eml_k];
    }
    for(eml_k = 2; eml_k < 6; eml_k++) {
      eml_dbuffer[eml_k - 1] -= eml_dbuffer[0] * eml_a[eml_k - 1];
    }
    eml_y_data[eml_j - 1] = eml_dbuffer[0];
  }
  for(eml_k = 0; eml_k < 4; eml_k++) {
    eml_zf[eml_k] = eml_dbuffer[eml_k + 1];
  }
}
static void m_c_power(real_T eml_b_data[6912], int32_T eml_b_sizes[2], real_T
 eml_y_data[6912], int32_T eml_y_sizes[2])
{
  int32_T eml_loop_ub;
  int32_T eml_k;
  m_d_eml_scalexp_alloc(eml_b_data, eml_b_sizes, eml_y_data, eml_y_sizes);
  eml_loop_ub = eml_y_sizes[1];
  for(eml_k = 1; eml_k <= eml_loop_ub; eml_k++) {
    eml_y_data[eml_k - 1] = rt_pow_snf(10.0, eml_b_data[eml_k - 1]);
  }
}
static int32_T m_compute_nones(boolean_T eml_x_data[6912], int32_T
 eml_x_sizes[2], int32_T eml_n)
{
  int32_T eml_k;
  int32_T eml_i;
  eml_k = 0;
  for(eml_i = 1; eml_i <= eml_n; eml_i++) {
    if(eml_x_data[eml_i - 1]) {
      eml_k++;
    }
  }
  return eml_k;
}
static void m_d_eml_scalexp_alloc(real_T eml_varargin_2_data[6912], int32_T
 eml_varargin_2_sizes[2], real_T eml_z_data[6912],
 int32_T eml_z_sizes[2])
{
  int32_T eml_i3;
  int16_T eml_iv17[2];
  int16_T eml_iv18[2];
  int16_T eml_tmp_sizes[2];
  if(!(eml_varargin_2_sizes[1] == 1)) {
    for(eml_i3 = 0; eml_i3 < 2; eml_i3++) {
      eml_iv17[eml_i3] = (int16_T)eml_varargin_2_sizes[eml_i3];
    }
    eml_iv18[0] = 1;
    eml_iv18[1] = eml_iv17[1];
    eml_tmp_sizes[0] = 1;
    eml_tmp_sizes[1] = eml_iv18[1];
    for(eml_i3 = 0; eml_i3 < 2; eml_i3++) {
      eml_z_sizes[eml_i3] = eml_tmp_sizes[eml_i3];
    }
  } else {
    for(eml_i3 = 0; eml_i3 < 2; eml_i3++) {
      eml_z_sizes[eml_i3] = 1;
    }
  }
}
static void m_eml_li_find(boolean_T eml_x_data[6912], int32_T eml_x_sizes[2],
 int32_T eml_y_data[6912], int32_T eml_y_sizes[2])
{
  int32_T eml_n;
  int32_T eml_iv4[2];
  int32_T eml_tmp_sizes[2];
  int32_T eml_j;
  int32_T eml_i;
  eml_n = eml_x_sizes[1];
  eml_iv4[0] = 1;
  eml_iv4[1] = m_compute_nones(eml_x_data, eml_x_sizes, eml_n);
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_iv4[1];
  for(eml_j = 0; eml_j < 2; eml_j++) {
    eml_y_sizes[eml_j] = eml_tmp_sizes[eml_j];
  }
  eml_j = 1;
  for(eml_i = 1; eml_i <= eml_n; eml_i++) {
    if(eml_x_data[eml_i - 1]) {
      eml_y_data[eml_j - 1] = eml_i;
      eml_j++;
    }
  }
}
static void m_eml_scalexp_alloc(real_T eml_varargin_1_data[6912], int32_T
 eml_varargin_1_sizes[2], real_T eml_z_data[6912], int32_T
 eml_z_sizes[2])
{
  int32_T eml_i0;
  int16_T eml_iv0[2];
  int16_T eml_iv1[2];
  int16_T eml_tmp_sizes[2];
  if(!(eml_varargin_1_sizes[1] == 1)) {
    for(eml_i0 = 0; eml_i0 < 2; eml_i0++) {
      eml_iv0[eml_i0] = (int16_T)eml_varargin_1_sizes[eml_i0];
    }
    eml_iv1[0] = 1;
    eml_iv1[1] = eml_iv0[1];
    eml_tmp_sizes[0] = 1;
    eml_tmp_sizes[1] = eml_iv1[1];
    for(eml_i0 = 0; eml_i0 < 2; eml_i0++) {
      eml_z_sizes[eml_i0] = eml_tmp_sizes[eml_i0];
    }
  } else {
    for(eml_i0 = 0; eml_i0 < 2; eml_i0++) {
      eml_z_sizes[eml_i0] = 1;
    }
  }
}
static void m_filter(real_T eml_b[2], real_T eml_a[2], real_T eml_x_data[6912],
 int32_T eml_x_sizes[2], real_T eml_zi, real_T
 eml_y_data[6912], int32_T eml_y_sizes[2], real_T *eml_zf)
{
  real_T eml_a1;
  int32_T eml_k;
  int16_T eml_iv2[2];
  int16_T eml_iv3[2];
  int16_T eml_tmp_sizes[2];
  real_T eml_dbuffer[2];
  int32_T eml_j;
  int32_T eml_b_k;
  eml_a1 = eml_a[0];
  if((!m_isfinite(eml_a1)) || (eml_a1 == 0.0) || (!(eml_a1 != 1.0))) {
  } else {
    for(eml_k = 0; eml_k < 2; eml_k++) {
      eml_b[eml_k] /= eml_a1;
    }
    eml_a[1] /= eml_a1;
    eml_a[0] = 1.0;
  }
  for(eml_k = 0; eml_k < 2; eml_k++) {
    eml_iv2[eml_k] = (int16_T)eml_x_sizes[eml_k];
  }
  eml_iv3[0] = 1;
  eml_iv3[1] = eml_iv2[1];
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_iv3[1];
  for(eml_k = 0; eml_k < 2; eml_k++) {
    eml_y_sizes[eml_k] = eml_tmp_sizes[eml_k];
  }
  eml_k = eml_x_sizes[1];
  eml_dbuffer[1] = eml_zi;
  for(eml_j = 1; eml_j <= eml_k; eml_j++) {
    eml_dbuffer[0] = eml_dbuffer[1];
    eml_dbuffer[1] = 0.0;
    for(eml_b_k = 0; eml_b_k < 2; eml_b_k++) {
      eml_dbuffer[eml_b_k] += eml_x_data[eml_j - 1] * eml_b[eml_b_k];
    }
    eml_dbuffer[1] -= eml_dbuffer[0] * eml_a[1];
    eml_y_data[eml_j - 1] = eml_dbuffer[0];
  }
  *eml_zf = eml_dbuffer[1];
}
static void m_find(boolean_T eml_x_data[6912], int32_T eml_x_sizes[2], real_T
 eml_i_data[6912], int32_T eml_i_sizes[2])
{
  int32_T eml_nx;
  int32_T eml_idx;
  int32_T eml_iv11[2];
  int32_T eml_tmp_sizes[2];
  int32_T eml_ii;
  boolean_T eml_exitg1;
  boolean_T eml_guard1 = FALSE;
  int32_T eml_hoistedExpr_sizes[1];
  static uint32_T eml_hoistedExpr_data[6912];
  int32_T eml_b_tmp_sizes[2];
  static uint32_T eml_tmp_data[6912];
  int32_T eml_b_i_sizes[2];
  static int32_T eml_b_i_data[6912];
  eml_nx = eml_x_sizes[1];
  eml_idx = 0;
  eml_iv11[0] = 1;
  eml_iv11[1] = eml_nx;
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_iv11[1];
  for(eml_ii = 0; eml_ii < 2; eml_ii++) {
    eml_i_sizes[eml_ii] = eml_tmp_sizes[eml_ii];
  }
  eml_ii = 1;
  eml_exitg1 = 0U;
  while((eml_exitg1 == 0U) && (eml_ii <= eml_nx)) {
    eml_guard1 = FALSE;
    if(eml_x_data[eml_ii - 1]) {
      eml_idx++;
      eml_i_data[eml_idx - 1] = (real_T)eml_ii;
      if(eml_idx >= eml_nx) {
        eml_exitg1 = 1U;
      } else {
        eml_guard1 = TRUE;
      }
    } else {
      eml_guard1 = TRUE;
    }
    if(eml_guard1 == TRUE) {
      eml_ii++;
    }
  }
  if(eml_nx == 1) {
    if(eml_idx == 0) {
      eml_i_sizes[0] = 1;
      eml_i_sizes[1] = 0;
      /* Empty Loop. */
    }
  } else {
    if(1 > eml_idx) {
      eml_idx = 0;
    }
    eml_hoistedExpr_sizes[0] = eml_idx;
    eml_ii = eml_idx - 1;
    for(eml_nx = 0; eml_nx <= eml_ii; eml_nx++) {
      eml_hoistedExpr_data[eml_nx] = (uint32_T)(1 + eml_nx);
    }
    eml_b_tmp_sizes[0] = 1;
    eml_iv11[0] = 1;
    eml_iv11[1] = eml_hoistedExpr_sizes[0];
    eml_b_tmp_sizes[1] = eml_iv11[1];
    eml_ii = eml_hoistedExpr_sizes[0] - 1;
    for(eml_nx = 0; eml_nx <= eml_ii; eml_nx++) {
      eml_tmp_data[eml_nx] = eml_hoistedExpr_data[eml_nx];
    }
    eml_iv11[0] = 1;
    eml_iv11[1] = eml_b_tmp_sizes[1];
    eml_b_i_sizes[0] = 1;
    eml_b_i_sizes[1] = eml_iv11[1];
    eml_nx = eml_iv11[1] - 1;
    for(eml_ii = 0; eml_ii <= eml_nx; eml_ii++) {
      for(eml_idx = 0; eml_idx <= 0; eml_idx = 1) {
        eml_b_i_data[eml_ii] = (int32_T)eml_i_data[(int32_T)eml_tmp_data[eml_ii]
          - 1];
      }
    }
    eml_i_sizes[0] = 1;
    eml_i_sizes[1] = eml_b_i_sizes[1];
    eml_ii = eml_b_i_sizes[1] - 1;
    for(eml_nx = 0; eml_nx <= eml_ii; eml_nx++) {
      eml_i_data[eml_nx] = (real_T)eml_b_i_data[eml_nx];
    }
  }
}
static boolean_T m_isfinite(real_T eml_x)
{
  return (!rtIsInf(eml_x)) && (!rtIsNaN(eml_x));
}
static real_T m_max(real_T eml_varargin_1_data[6912], int32_T
 eml_varargin_1_sizes[2])
{
  real_T eml_maxval;
  int32_T eml_n;
  int32_T eml_itmp;
  int32_T eml_ix;
  boolean_T eml_guard1 = FALSE;
  boolean_T eml_searchingForNonNaN;
  int32_T eml_k;
  boolean_T eml_exitg1;
  eml_n = eml_varargin_1_sizes[1];
  eml_maxval = eml_varargin_1_data[0];
  eml_itmp = 1;
  if(eml_n == 1) {
  } else {
    eml_ix = 1;
    eml_guard1 = FALSE;
    if(rtIsNaN(eml_maxval)) {
      eml_searchingForNonNaN = TRUE;
      eml_k = 2;
      eml_exitg1 = 0U;
      while((eml_exitg1 == 0U) && (eml_k <= eml_n)) {
        eml_ix++;
        if(!rtIsNaN(eml_varargin_1_data[eml_ix - 1])) {
          eml_maxval = eml_varargin_1_data[eml_ix - 1];
          eml_itmp = eml_k;
          eml_searchingForNonNaN = FALSE;
          eml_exitg1 = 1U;
        } else {
          eml_k++;
        }
      }
      if(eml_searchingForNonNaN) {
      } else {
        eml_guard1 = TRUE;
      }
    } else {
      eml_guard1 = TRUE;
    }
    if(eml_guard1 == TRUE) {
      for(eml_itmp++; eml_itmp <= eml_n; eml_itmp++) {
        eml_ix++;
        if(eml_varargin_1_data[eml_ix - 1] > eml_maxval) {
          eml_maxval = eml_varargin_1_data[eml_ix - 1];
        }
      }
    }
  }
  return eml_maxval;
}
static void m_power(real_T eml_a_data[6912], int32_T eml_a_sizes[2], real_T
 eml_y_data[6912], int32_T eml_y_sizes[2])
{
  int32_T eml_loop_ub;
  int32_T eml_k;
  m_eml_scalexp_alloc(eml_a_data, eml_a_sizes, eml_y_data, eml_y_sizes);
  eml_loop_ub = eml_y_sizes[1];
  for(eml_k = 1; eml_k <= eml_loop_ub; eml_k++) {
    eml_y_data[eml_k - 1] = rt_pow_snf(eml_a_data[eml_k - 1], 2.0);
  }
}
static void m_refp1_DRNL_brokenstick_nl(real_T eml_x_data[6912], int32_T
 eml_x_sizes[2], real_T eml_b, real_T eml_c)
{
  int32_T eml_abs_x_sizes[2];
  static real_T eml_abs_x_data[6912];
  real_T eml_compressionThreshold;
  int32_T eml_b_abs_x_sizes[2];
  int32_T eml_loop_ub;
  int32_T eml_i4;
  boolean_T eml_b_abs_x_data[6912];
  int32_T eml_idx_sizes[2];
  static real_T eml_idx_data[6912];
  int32_T eml_c_abs_x_sizes[2];
  static real_T eml_c_abs_x_data[6912];
  int32_T eml_hoistedExpr_sizes[2];
  int32_T eml_i5;
  static real_T eml_hoistedExpr_data[6912];
  /* nick modified broken stick function */
  /*  y = sign(x).* min(a*abs_x,  b*abs_x .^ c); */
  /*  This function could be replaced by a lookup table */
  m_abs(eml_x_data, eml_x_sizes, eml_abs_x_data, eml_abs_x_sizes);
  /*  linear (low amplitude) response */
  eml_x_sizes[0] = 1;
  /*  compressed high amplitude */
  eml_compressionThreshold = rt_pow_snf(10.0, 1.0 / (1.0 - eml_c) *
   log10(eml_b));
  /*  only values outside the compression threshold */
  /*   need be subject to compression */
  eml_b_abs_x_sizes[0] = 1;
  eml_b_abs_x_sizes[1] = eml_abs_x_sizes[1];
  eml_loop_ub = eml_abs_x_sizes[0] * eml_abs_x_sizes[1] - 1;
  for(eml_i4 = 0; eml_i4 <= eml_loop_ub; eml_i4++) {
    eml_b_abs_x_data[eml_i4] = (eml_abs_x_data[eml_i4] >
      eml_compressionThreshold);
  }
  m_find(eml_b_abs_x_data, eml_b_abs_x_sizes, eml_idx_data, eml_idx_sizes);
  eml_c_abs_x_sizes[0] = 1;
  eml_b_abs_x_sizes[0] = 1;
  eml_b_abs_x_sizes[1] = eml_idx_sizes[1];
  eml_c_abs_x_sizes[1] = eml_b_abs_x_sizes[1];
  eml_loop_ub = eml_idx_sizes[1] - 1;
  for(eml_i4 = 0; eml_i4 <= eml_loop_ub; eml_i4++) {
    eml_c_abs_x_data[eml_i4] = eml_abs_x_data[(int32_T)eml_idx_data[eml_i4] - 1];
  }
  m_b_power(eml_c_abs_x_data, eml_c_abs_x_sizes, eml_c, eml_abs_x_data,
   eml_abs_x_sizes);
  eml_b_abs_x_sizes[0] = 1;
  eml_b_abs_x_sizes[1] = eml_idx_sizes[1];
  eml_hoistedExpr_sizes[0] = 1;
  eml_hoistedExpr_sizes[1] = eml_b_abs_x_sizes[1];
  eml_loop_ub = eml_b_abs_x_sizes[1] - 1;
  for(eml_i4 = 0; eml_i4 <= eml_loop_ub; eml_i4++) {
    for(eml_i5 = 0; eml_i5 <= 0; eml_i5 = 1) {
      eml_hoistedExpr_data[eml_i4] = eml_x_data[(int32_T)eml_idx_data[eml_i4] -
        1];
    }
  }
  m_refp1_sign(eml_hoistedExpr_data, eml_hoistedExpr_sizes);
  eml_loop_ub = eml_hoistedExpr_sizes[1] - 1;
  for(eml_i4 = 0; eml_i4 <= eml_loop_ub; eml_i4++) {
    eml_x_data[(int32_T)eml_idx_data[eml_i4] - 1] =
      eml_hoistedExpr_data[eml_hoistedExpr_sizes[0] * eml_i4] * (eml_b *
      eml_abs_x_data[eml_abs_x_sizes[0] * eml_i4]);
  }
  /* of DRNL_brokenstick_nl */
}
static void m_refp1_log10(real_T eml_x_data[6912], int32_T eml_x_sizes[2])
{
  int32_T eml_b_x_sizes[2];
  int32_T eml_loop_ub;
  int32_T eml_k;
  eml_b_x_sizes[0] = 1;
  eml_b_x_sizes[1] = eml_x_sizes[1];
  eml_loop_ub = eml_b_x_sizes[1];
  for(eml_k = 1; eml_k <= eml_loop_ub; eml_k++) {
    eml_x_data[eml_k - 1] = log10(eml_x_data[eml_k - 1]);
  }
}
static void m_refp1_sign(real_T eml_x_data[6912], int32_T eml_x_sizes[2])
{
  int32_T eml_b_x_sizes[2];
  int32_T eml_loop_ub;
  int32_T eml_k;
  real_T eml_x;
  eml_b_x_sizes[0] = 1;
  eml_b_x_sizes[1] = eml_x_sizes[1];
  eml_loop_ub = eml_b_x_sizes[1];
  for(eml_k = 1; eml_k <= eml_loop_ub; eml_k++) {
    eml_x = eml_x_data[eml_k - 1];
    if(rtIsNaN(eml_x)) {
      eml_x = rtNaN;
    } else if(eml_x > 0.0) {
      eml_x = 1.0;
    } else if(eml_x < 0.0) {
      eml_x = -1.0;
    } else {
      eml_x = 0.0;
    }
    eml_x_data[eml_k - 1] = eml_x;
  }
}
static void m_refp1_sqrt(real_T eml_x_data[6912], int32_T eml_x_sizes[2])
{
  int32_T eml_b_x_sizes[2];
  int32_T eml_loop_ub;
  int32_T eml_k;
  eml_b_x_sizes[0] = 1;
  eml_b_x_sizes[1] = eml_x_sizes[1];
  eml_loop_ub = eml_b_x_sizes[1];
  for(eml_k = 1; eml_k <= eml_loop_ub; eml_k++) {
    eml_x_data[eml_k - 1] = sqrt(eml_x_data[eml_k - 1]);
  }
}
static void m_sum(real_T eml_x_data[13824], int32_T eml_x_sizes[2], real_T
 eml_y_data[6912], int32_T eml_y_sizes[2])
{
  int32_T eml_b;
  int16_T eml_sz[2];
  int16_T eml_iv14[2];
  int16_T eml_tmp_sizes[2];
  int32_T eml_ix;
  int32_T eml_iy;
  int32_T eml_i;
  int32_T eml_ixstart;
  for(eml_b = 0; eml_b < 2; eml_b++) {
    eml_sz[eml_b] = (int16_T)eml_x_sizes[eml_b];
  }
  eml_sz[0] = 1;
  eml_iv14[0] = 1;
  eml_iv14[1] = eml_sz[1];
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_iv14[1];
  for(eml_b = 0; eml_b < 2; eml_b++) {
    eml_y_sizes[eml_b] = eml_tmp_sizes[eml_b];
  }
  if(eml_x_sizes[1] == 0) {
    eml_b = eml_y_sizes[1] - 1;
    for(eml_ix = 0; eml_ix <= eml_b; eml_ix++) {
      eml_iy = eml_y_sizes[0] - 1;
      for(eml_i = 0; eml_i <= eml_iy; eml_i++) {
        eml_y_data[eml_i + eml_y_sizes[0] * eml_ix] = 0.0;
      }
    }
  } else {
    eml_b = eml_x_sizes[1];
    eml_ix = 0;
    eml_iy = 0;
    for(eml_i = 1; eml_i <= eml_b; eml_i++) {
      eml_ixstart = eml_ix + 1;
      eml_ix = eml_ixstart + 1;
      eml_iy++;
      eml_y_data[eml_iy - 1] = eml_x_data[eml_ixstart - 1] + eml_x_data[eml_ix -
        1];
    }
  }
}
void EssexAidProcessVFrameSwitchable(real_T eml_frameBufferL[6912], real_T
 eml_frameBufferR[6912], real_T eml_filterStatesL[3000],
 real_T eml_filterStatesR[3000], const real_T eml_filterCoeffs[5000], real_T
 eml_numChannels, real_T eml_numSamples, real_T
 eml_ARampL[6912], real_T eml_ARampR[6912], real_T eml_ARthresholdPa, real_T
 eml_filterOrder, const real_T eml_DRNLb[11], const
 real_T eml_DRNLc[11], const real_T eml_MOCthreshold_dBOP[11], real_T
 eml_MOCfactor, real_T eml_peakIPL[11], real_T eml_peakOPL[
 11], real_T eml_rmsIPL[11], real_T eml_rmsOPL[11], real_T eml_peakIPR[11],
 real_T eml_peakOPR[11], real_T eml_rmsIPR[11], real_T
 eml_rmsOPR[11], real_T eml_MOCend[11], real_T eml_MOCcontrol[76032], const
 real_T eml_mainGain[11], boolean_T eml_useGTF)
{
  int32_T eml_loop_ub;
  int32_T eml_hoistedExpr_sizes[1];
  int32_T eml_b_loop_ub;
  static real_T eml_hoistedExpr_data[6912];
  int32_T eml_b_hoistedExpr_sizes[2];
  int32_T eml_tmp_sizes[2];
  static real_T eml_b_hoistedExpr_data[6912];
  int32_T eml_stapesVelL_sizes[2];
  int32_T eml_y;
  static real_T eml_stapesVelL_data[6912];
  int32_T eml_c_hoistedExpr_sizes[1];
  static real_T eml_c_hoistedExpr_data[6912];
  int32_T eml_stapesVelR_sizes[2];
  static real_T eml_stapesVelR_data[6912];
  int32_T eml_d_hoistedExpr_sizes[1];
  static real_T eml_yR_data[6912];
  int32_T eml_frameBufferL_sizes[2];
  int16_T eml_iv19[2];
  int32_T eml_e_hoistedExpr_sizes[1];
  static real_T eml_yL_data[6912];
  int16_T eml_iv20[2];
  int32_T eml_f_hoistedExpr_sizes[1];
  static real_T eml_d_hoistedExpr_data[6912];
  int32_T eml_b_frameBufferL_sizes[2];
  int32_T eml_g_hoistedExpr_sizes[1];
  static real_T eml_b_yR_data[6912];
  real_T eml_b_filterCoeffs[2];
  real_T eml_c_filterCoeffs[2];
  real_T eml_d_filterCoeffs[2];
  real_T eml_e_filterCoeffs[2];
  real_T eml_b_MOCthreshold_dBOP;
  int32_T eml_yR_sizes[2];
  int32_T eml_h_hoistedExpr_sizes[1];
  static real_T eml_e_hoistedExpr_data[6912];
  int32_T eml_i_hoistedExpr_sizes[1];
  static real_T eml_f_hoistedExpr_data[6912];
  int32_T eml_j_hoistedExpr_sizes[1];
  static real_T eml_g_hoistedExpr_data[6912];
  int32_T eml_k_hoistedExpr_sizes[1];
  static real_T eml_h_hoistedExpr_data[6912];
  int32_T eml_l_hoistedExpr_sizes[1];
  static real_T eml_i_hoistedExpr_data[6912];
  int32_T eml_yL_sizes[2];
  boolean_T eml_b_yL_data[6912];
  static int32_T eml_tmp_data[6912];
  int32_T eml_b_tmp_sizes[2];
  static int32_T eml_b_tmp_data[6912];
  int32_T eml_m_hoistedExpr_sizes[1];
  static real_T eml_j_hoistedExpr_data[6912];
  int32_T eml_b_yL_sizes[2];
  int32_T eml_filterCount;
  int32_T eml_b_yR_sizes[2];
  static real_T eml_c_yR_data[6912];
  int32_T eml_n_hoistedExpr_sizes[1];
  static real_T eml_k_hoistedExpr_data[6912];
  real_T eml_c_loop_ub;
  real_T eml_d0;
  int32_T eml_o_hoistedExpr_sizes[1];
  real_T eml_l_hoistedExpr_data[182];
  int32_T eml_c_tmp_sizes[2];
  real_T eml_c_tmp_data[182];
  real_T eml_f_filterCoeffs[3];
  int32_T eml_p_hoistedExpr_sizes[2];
  real_T eml_m_hoistedExpr_data[182];
  int32_T eml_filterStatesL_sizes[1];
  real_T eml_filterStatesL_data[190];
  int32_T eml_c_yR_sizes[2];
  int32_T eml_q_hoistedExpr_sizes[1];
  real_T eml_n_hoistedExpr_data[182];
  int32_T eml_d_tmp_sizes[2];
  int32_T eml_d_tmp_data[182];
  int32_T eml_r_hoistedExpr_sizes[1];
  real_T eml_o_hoistedExpr_data[182];
  int32_T eml_s_hoistedExpr_sizes[2];
  real_T eml_p_hoistedExpr_data[182];
  int32_T eml_filterStatesR_sizes[1];
  real_T eml_filterStatesR_data[190];
  int32_T eml_d_yR_sizes[2];
  int32_T eml_t_hoistedExpr_sizes[1];
  real_T eml_q_hoistedExpr_data[182];
  real_T eml_nn;
  int32_T eml_u_hoistedExpr_sizes[1];
  real_T eml_r_hoistedExpr_data[180];
  int32_T eml_e_tmp_sizes[2];
  real_T eml_e_tmp_data[180];
  real_T eml_g_filterCoeffs[5];
  real_T eml_h_filterCoeffs[5];
  int32_T eml_v_hoistedExpr_sizes[2];
  real_T eml_s_hoistedExpr_data[180];
  int32_T eml_b_filterStatesL_sizes[1];
  static real_T eml_b_filterStatesL_data[3000];
  int32_T eml_e_yR_sizes[2];
  real_T eml_dv0[4];
  int32_T eml_w_hoistedExpr_sizes[1];
  real_T eml_t_hoistedExpr_data[180];
  int32_T eml_f_tmp_sizes[2];
  int32_T eml_f_tmp_data[180];
  int32_T eml_x_hoistedExpr_sizes[1];
  real_T eml_u_hoistedExpr_data[180];
  int32_T eml_y_hoistedExpr_sizes[2];
  real_T eml_v_hoistedExpr_data[180];
  int32_T eml_b_filterStatesR_sizes[1];
  real_T eml_b_filterStatesR_data[3000];
  int32_T eml_f_yR_sizes[2];
  int32_T eml_ab_hoistedExpr_sizes[1];
  real_T eml_w_hoistedExpr_data[180];
  int32_T eml_bb_hoistedExpr_sizes[1];
  real_T eml_x_hoistedExpr_data[190];
  int32_T eml_g_tmp_sizes[2];
  real_T eml_g_tmp_data[190];
  int32_T eml_cb_hoistedExpr_sizes[2];
  real_T eml_y_hoistedExpr_data[190];
  int32_T eml_c_filterStatesL_sizes[1];
  real_T eml_c_filterStatesL_data[190];
  int32_T eml_g_yR_sizes[2];
  int32_T eml_db_hoistedExpr_sizes[1];
  real_T eml_ab_hoistedExpr_data[190];
  int32_T eml_h_tmp_sizes[2];
  int32_T eml_h_tmp_data[190];
  int32_T eml_eb_hoistedExpr_sizes[1];
  real_T eml_bb_hoistedExpr_data[190];
  int32_T eml_fb_hoistedExpr_sizes[2];
  real_T eml_cb_hoistedExpr_data[190];
  int32_T eml_c_filterStatesR_sizes[1];
  real_T eml_c_filterStatesR_data[190];
  int32_T eml_h_yR_sizes[2];
  int32_T eml_gb_hoistedExpr_sizes[1];
  real_T eml_db_hoistedExpr_data[190];
  int32_T eml_hb_hoistedExpr_sizes[1];
  static real_T eml_eb_hoistedExpr_data[3000];
  int32_T eml_i_tmp_sizes[2];
  static real_T eml_i_tmp_data[3000];
  int32_T eml_ib_hoistedExpr_sizes[2];
  static real_T eml_fb_hoistedExpr_data[3000];
  int32_T eml_d_filterStatesL_sizes[1];
  static real_T eml_d_filterStatesL_data[3000];
  int32_T eml_i_yR_sizes[2];
  int32_T eml_jb_hoistedExpr_sizes[1];
  static real_T eml_gb_hoistedExpr_data[3000];
  int32_T eml_j_tmp_sizes[2];
  int32_T eml_j_tmp_data[3000];
  int32_T eml_kb_hoistedExpr_sizes[1];
  real_T eml_hb_hoistedExpr_data[3000];
  int32_T eml_lb_hoistedExpr_sizes[2];
  real_T eml_ib_hoistedExpr_data[3000];
  int32_T eml_d_filterStatesR_sizes[1];
  real_T eml_d_filterStatesR_data[3000];
  int32_T eml_j_yR_sizes[2];
  int32_T eml_mb_hoistedExpr_sizes[1];
  real_T eml_jb_hoistedExpr_data[3000];
  int32_T eml_nb_hoistedExpr_sizes[2];
  static real_T eml_kb_hoistedExpr_data[13824];
  int32_T eml_ob_hoistedExpr_sizes[2];
  int32_T eml_c_yL_sizes[2];
  int32_T eml_pb_hoistedExpr_sizes[1];
  static real_T eml_lb_hoistedExpr_data[6912];
  int32_T eml_qb_hoistedExpr_sizes[2];
  int32_T eml_rb_hoistedExpr_sizes[1];
  static real_T eml_mb_hoistedExpr_data[6912];
  int32_T eml_sb_hoistedExpr_sizes[1];
  static real_T eml_nb_hoistedExpr_data[6912];
  int32_T eml_c_frameBufferL_sizes[2];
  int32_T eml_tb_hoistedExpr_sizes[1];
  static real_T eml_ob_hoistedExpr_data[6912];
  int32_T eml_ub_hoistedExpr_sizes[1];
  static real_T eml_pb_hoistedExpr_data[6912];
  int32_T eml_frameBufferR_sizes[2];
  /* ESSEXAIDPROCESSFRAME Essex aid algorithm in frame processing mode */
  /*    This code will look a bit odd to most Matlab programmers. This is */
  /*    because the intended target is a C function that will be called on a */
  /*    sub-millisecond basis. The bizzare enumerations assist the */
  /*    pass-by-reference functionality that allows this function to fly in */
  /*    real-time. This function works on a need to know basis, eliminating any */
  /*    unnecessary data copying or parameter calculation. */
  /*  eml.varsize('frameBuffer', 6192); */
  /* Fake enumeration - must be kept up to date with juce enum */
  /*  enumC_BPb1 = 8; */
  /*  enumC_BPa1 = 13; */
  /*  enumC_BPb2 = 18; */
  /*  enumC_BPa2 = 23; */
  /*  enumC_BPb3 = 28; */
  /*  enumC_BPa3 = 33; */
  /*  enumC_BPb4 = 38; */
  /*  enumC_BPa4 = 43; */
  /*  enumS_MOC1  = 1; */
  /*  enumS_BPin_1_1 = 2; */
  /*  enumS_BPin_2_1 = 6; */
  /*  enumS_BPout_1_1 = 10; */
  /*  enumS_BPout_2_1 = 14; */
  /*  */
  /*  enumS_MOC2 = 18; */
  /*  enumS_BPin_1_2 = 19; */
  /*  enumS_BPin_2_2 = 23; */
  /*  enumS_BPout_1_2 = 27; */
  /*  enumS_BPout_2_2 = 31; */
  /*  ... */
  /*  rmsLev[0] = iunput RMS from AR smoothed response */
  /* Initial gain */
  /*  frameBuffer(1:numSamples) = frameBuffer(1:numSamples)*ipScalar; */
  /* % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  Place where conversion to velocity once lived */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  tympanic membrane response in meters */
  if(1.0 > eml_numSamples) {
    eml_loop_ub = 0;
  } else {
    eml_loop_ub = (int32_T)eml_numSamples;
  }
  eml_hoistedExpr_sizes[0] = eml_loop_ub;
  eml_loop_ub--;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_hoistedExpr_data[eml_b_loop_ub] = 1.0 + (real_T)eml_b_loop_ub;
  }
  eml_b_hoistedExpr_sizes[0] = 1;
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_hoistedExpr_sizes[0];
  eml_b_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
  eml_loop_ub = eml_hoistedExpr_sizes[0] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_b_hoistedExpr_data[eml_b_loop_ub] = eml_hoistedExpr_data[eml_b_loop_ub];
  }
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_b_hoistedExpr_sizes[1];
  eml_stapesVelL_sizes[0] = 1;
  eml_stapesVelL_sizes[1] = eml_tmp_sizes[1];
  eml_loop_ub = eml_tmp_sizes[1] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    for(eml_y = 0; eml_y <= 0; eml_y = 1) {
      eml_stapesVelL_data[eml_b_loop_ub] =
        eml_frameBufferL[(int32_T)eml_b_hoistedExpr_data[eml_b_loop_ub] - 1];
    }
  }
  if(1.0 > eml_numSamples) {
    eml_loop_ub = 0;
  } else {
    eml_loop_ub = (int32_T)eml_numSamples;
  }
  eml_c_hoistedExpr_sizes[0] = eml_loop_ub;
  eml_loop_ub--;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_c_hoistedExpr_data[eml_b_loop_ub] = 1.0 + (real_T)eml_b_loop_ub;
  }
  eml_b_hoistedExpr_sizes[0] = 1;
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_c_hoistedExpr_sizes[0];
  eml_b_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
  eml_loop_ub = eml_c_hoistedExpr_sizes[0] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_b_hoistedExpr_data[eml_b_loop_ub] =
      eml_c_hoistedExpr_data[eml_b_loop_ub];
  }
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_b_hoistedExpr_sizes[1];
  eml_stapesVelR_sizes[0] = 1;
  eml_stapesVelR_sizes[1] = eml_tmp_sizes[1];
  eml_loop_ub = eml_tmp_sizes[1] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    for(eml_y = 0; eml_y <= 0; eml_y = 1) {
      eml_stapesVelR_data[eml_b_loop_ub] =
        eml_frameBufferR[(int32_T)eml_b_hoistedExpr_data[eml_b_loop_ub] - 1];
    }
  }
  if(1.0 > eml_numSamples) {
    eml_loop_ub = 0;
  } else {
    eml_loop_ub = (int32_T)eml_numSamples;
  }
  eml_d_hoistedExpr_sizes[0] = eml_loop_ub;
  eml_loop_ub--;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_yR_data[eml_b_loop_ub] = 1.0 + (real_T)eml_b_loop_ub;
  }
  eml_b_hoistedExpr_sizes[0] = 1;
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_d_hoistedExpr_sizes[0];
  eml_b_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
  eml_loop_ub = eml_d_hoistedExpr_sizes[0] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_b_hoistedExpr_data[eml_b_loop_ub] = eml_yR_data[eml_b_loop_ub];
  }
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_b_hoistedExpr_sizes[1];
  eml_frameBufferL_sizes[0] = 1;
  eml_frameBufferL_sizes[1] = eml_tmp_sizes[1];
  for(eml_loop_ub = 0; eml_loop_ub < 2; eml_loop_ub++) {
    eml_iv19[eml_loop_ub] = (int16_T)eml_frameBufferL_sizes[eml_loop_ub];
  }
  if(1.0 > eml_numSamples) {
    eml_loop_ub = 0;
  } else {
    eml_loop_ub = (int32_T)eml_numSamples;
  }
  eml_e_hoistedExpr_sizes[0] = eml_loop_ub;
  eml_loop_ub--;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_yL_data[eml_b_loop_ub] = 1.0 + (real_T)eml_b_loop_ub;
  }
  eml_b_hoistedExpr_sizes[0] = 1;
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_e_hoistedExpr_sizes[0];
  eml_b_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
  eml_loop_ub = eml_e_hoistedExpr_sizes[0] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_b_hoistedExpr_data[eml_b_loop_ub] = eml_yL_data[eml_b_loop_ub];
  }
  eml_iv20[0] = 1;
  eml_iv20[1] = eml_iv19[1];
  eml_loop_ub = eml_iv20[1] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    for(eml_y = 0; eml_y <= 0; eml_y = 1) {
      eml_frameBufferL[(int32_T)eml_b_hoistedExpr_data[eml_b_loop_ub] - 1] = 0.0;
    }
  }
  if(1.0 > eml_numSamples) {
    eml_loop_ub = 0;
  } else {
    eml_loop_ub = (int32_T)eml_numSamples;
  }
  eml_f_hoistedExpr_sizes[0] = eml_loop_ub;
  eml_loop_ub--;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_d_hoistedExpr_data[eml_b_loop_ub] = 1.0 + (real_T)eml_b_loop_ub;
  }
  eml_b_hoistedExpr_sizes[0] = 1;
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_f_hoistedExpr_sizes[0];
  eml_b_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
  eml_loop_ub = eml_f_hoistedExpr_sizes[0] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_b_hoistedExpr_data[eml_b_loop_ub] =
      eml_d_hoistedExpr_data[eml_b_loop_ub];
  }
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_b_hoistedExpr_sizes[1];
  eml_b_frameBufferL_sizes[0] = 1;
  eml_b_frameBufferL_sizes[1] = eml_tmp_sizes[1];
  for(eml_loop_ub = 0; eml_loop_ub < 2; eml_loop_ub++) {
    eml_iv19[eml_loop_ub] = (int16_T)eml_b_frameBufferL_sizes[eml_loop_ub];
  }
  if(1.0 > eml_numSamples) {
    eml_loop_ub = 0;
  } else {
    eml_loop_ub = (int32_T)eml_numSamples;
  }
  eml_g_hoistedExpr_sizes[0] = eml_loop_ub;
  eml_loop_ub--;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_b_yR_data[eml_b_loop_ub] = 1.0 + (real_T)eml_b_loop_ub;
  }
  eml_loop_ub = eml_g_hoistedExpr_sizes[0] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_b_hoistedExpr_data[eml_b_loop_ub] = eml_b_yR_data[eml_b_loop_ub];
  }
  eml_iv20[0] = 1;
  eml_iv20[1] = eml_iv19[1];
  eml_loop_ub = eml_iv20[1] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    for(eml_y = 0; eml_y <= 0; eml_y = 1) {
      eml_frameBufferR[(int32_T)eml_b_hoistedExpr_data[eml_b_loop_ub] - 1] = 0.0;
    }
  }
  /* % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  ACOUSTIC REFLEX */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  find  rms of smoothed ip signal */
  /*   this will be used to trigger the AR reflex */
  m_power(eml_stapesVelL_data, eml_stapesVelL_sizes, eml_b_hoistedExpr_data,
   eml_b_hoistedExpr_sizes);
  for(eml_loop_ub = 0; eml_loop_ub < 2; eml_loop_ub++) {
    eml_b_filterCoeffs[eml_loop_ub] = eml_filterCoeffs[eml_loop_ub];
    eml_c_filterCoeffs[eml_loop_ub] = eml_filterCoeffs[eml_loop_ub + 2];
    eml_d_filterCoeffs[eml_loop_ub] = eml_filterCoeffs[eml_loop_ub];
    eml_e_filterCoeffs[eml_loop_ub] = eml_filterCoeffs[eml_loop_ub + 2];
  }
  m_filter(eml_b_filterCoeffs, eml_c_filterCoeffs, eml_b_hoistedExpr_data,
   eml_b_hoistedExpr_sizes, eml_filterStatesL[0],
   eml_yL_data, eml_b_frameBufferL_sizes, &eml_b_MOCthreshold_dBOP);
  eml_filterStatesL[0] = eml_b_MOCthreshold_dBOP;
  m_power(eml_stapesVelR_data, eml_stapesVelR_sizes, eml_b_hoistedExpr_data,
   eml_b_hoistedExpr_sizes);
  m_filter(eml_d_filterCoeffs, eml_e_filterCoeffs, eml_b_hoistedExpr_data,
   eml_b_hoistedExpr_sizes, eml_filterStatesR[0],
   eml_yR_data, eml_yR_sizes, &eml_b_MOCthreshold_dBOP);
  eml_filterStatesR[0] = eml_b_MOCthreshold_dBOP;
  /*  restore Pa scale */
  m_refp1_sqrt(eml_yL_data, eml_b_frameBufferL_sizes);
  /* confusing name for parameter - it is a short term RMS. */
  m_refp1_sqrt(eml_yR_data, eml_yR_sizes);
  /* confusing name for parameter - it is a short term RMS. */
  /*  attenuate input (NB cross product used) */
  if(1.0 > eml_numSamples) {
    eml_loop_ub = 0;
  } else {
    eml_loop_ub = (int32_T)eml_numSamples;
  }
  eml_h_hoistedExpr_sizes[0] = eml_loop_ub;
  eml_loop_ub--;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_e_hoistedExpr_data[eml_b_loop_ub] = 1.0 + (real_T)eml_b_loop_ub;
  }
  eml_b_hoistedExpr_sizes[0] = 1;
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_h_hoistedExpr_sizes[0];
  eml_b_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
  eml_loop_ub = eml_h_hoistedExpr_sizes[0] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_b_hoistedExpr_data[eml_b_loop_ub] =
      eml_e_hoistedExpr_data[eml_b_loop_ub];
  }
  eml_stapesVelL_sizes[0] = 1;
  eml_loop_ub = eml_stapesVelL_sizes[1] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_stapesVelL_data[eml_b_loop_ub] /=
      eml_ARampL[(int32_T)eml_b_hoistedExpr_data[eml_b_loop_ub] - 1];
  }
  if(1.0 > eml_numSamples) {
    eml_loop_ub = 0;
  } else {
    eml_loop_ub = (int32_T)eml_numSamples;
  }
  eml_i_hoistedExpr_sizes[0] = eml_loop_ub;
  eml_loop_ub--;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_f_hoistedExpr_data[eml_b_loop_ub] = 1.0 + (real_T)eml_b_loop_ub;
  }
  eml_b_hoistedExpr_sizes[0] = 1;
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_i_hoistedExpr_sizes[0];
  eml_b_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
  eml_loop_ub = eml_i_hoistedExpr_sizes[0] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_b_hoistedExpr_data[eml_b_loop_ub] =
      eml_f_hoistedExpr_data[eml_b_loop_ub];
  }
  eml_stapesVelR_sizes[0] = 1;
  eml_loop_ub = eml_stapesVelR_sizes[1] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_stapesVelR_data[eml_b_loop_ub] /=
      eml_ARampR[(int32_T)eml_b_hoistedExpr_data[eml_b_loop_ub] - 1];
  }
  /* CALC ARamp FOR NEXT FRAME */
  /*  compare levels in the previous segment with AR threshold */
  if(1.0 > eml_numSamples) {
    eml_loop_ub = 0;
  } else {
    eml_loop_ub = (int32_T)eml_numSamples;
  }
  eml_j_hoistedExpr_sizes[0] = eml_loop_ub;
  eml_loop_ub--;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_g_hoistedExpr_data[eml_b_loop_ub] = 1.0 + (real_T)eml_b_loop_ub;
  }
  eml_b_hoistedExpr_sizes[0] = 1;
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_j_hoistedExpr_sizes[0];
  eml_b_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
  eml_loop_ub = eml_j_hoistedExpr_sizes[0] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_b_hoistedExpr_data[eml_b_loop_ub] =
      eml_g_hoistedExpr_data[eml_b_loop_ub];
  }
  eml_loop_ub = eml_b_frameBufferL_sizes[1] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_ARampL[(int32_T)eml_b_hoistedExpr_data[eml_b_loop_ub] - 1] =
      eml_yL_data[eml_b_frameBufferL_sizes[0] * eml_b_loop_ub] /
      eml_ARthresholdPa;
  }
  if(1.0 > eml_numSamples) {
    eml_loop_ub = 0;
  } else {
    eml_loop_ub = (int32_T)eml_numSamples;
  }
  eml_k_hoistedExpr_sizes[0] = eml_loop_ub;
  eml_loop_ub--;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_h_hoistedExpr_data[eml_b_loop_ub] = 1.0 + (real_T)eml_b_loop_ub;
  }
  eml_b_hoistedExpr_sizes[0] = 1;
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_k_hoistedExpr_sizes[0];
  eml_b_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
  eml_loop_ub = eml_k_hoistedExpr_sizes[0] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_b_hoistedExpr_data[eml_b_loop_ub] =
      eml_h_hoistedExpr_data[eml_b_loop_ub];
  }
  eml_loop_ub = eml_yR_sizes[1] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_ARampR[(int32_T)eml_b_hoistedExpr_data[eml_b_loop_ub] - 1] =
      eml_yR_data[eml_yR_sizes[0] * eml_b_loop_ub] /
      eml_ARthresholdPa;
  }
  /*  all sub-treshold values are set to 1 */
  if(1.0 > eml_numSamples) {
    eml_loop_ub = 0;
  } else {
    eml_loop_ub = (int32_T)eml_numSamples;
  }
  eml_l_hoistedExpr_sizes[0] = eml_loop_ub;
  eml_loop_ub--;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_i_hoistedExpr_data[eml_b_loop_ub] = 1.0 + (real_T)eml_b_loop_ub;
  }
  eml_b_hoistedExpr_sizes[0] = 1;
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_l_hoistedExpr_sizes[0];
  eml_b_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
  eml_loop_ub = eml_l_hoistedExpr_sizes[0] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_b_hoistedExpr_data[eml_b_loop_ub] =
      eml_i_hoistedExpr_data[eml_b_loop_ub];
  }
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_b_hoistedExpr_sizes[1];
  eml_b_frameBufferL_sizes[0] = 1;
  eml_b_frameBufferL_sizes[1] = eml_tmp_sizes[1];
  eml_loop_ub = eml_tmp_sizes[1] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    for(eml_y = 0; eml_y <= 0; eml_y = 1) {
      eml_yL_data[eml_b_loop_ub] =
        eml_ARampL[(int32_T)eml_b_hoistedExpr_data[eml_b_loop_ub] - 1];
    }
  }
  eml_yL_sizes[0] = 1;
  eml_yL_sizes[1] = eml_b_frameBufferL_sizes[1];
  eml_loop_ub = eml_b_frameBufferL_sizes[1] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_b_yL_data[eml_b_loop_ub] = (eml_yL_data[eml_b_loop_ub] < 1.0);
  }
  m_eml_li_find(eml_b_yL_data, eml_yL_sizes, eml_tmp_data,
   eml_frameBufferL_sizes);
  eml_b_tmp_sizes[0] = 1;
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_frameBufferL_sizes[1];
  eml_b_tmp_sizes[1] = eml_tmp_sizes[1];
  eml_loop_ub = eml_frameBufferL_sizes[1] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_b_tmp_data[eml_b_loop_ub] = eml_tmp_data[eml_b_loop_ub];
  }
  eml_loop_ub = eml_b_tmp_sizes[1] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    for(eml_y = 0; eml_y <= 0; eml_y = 1) {
      eml_ARampL[eml_b_tmp_data[eml_b_loop_ub] - 1] = 1.0;
    }
  }
  if(1.0 > eml_numSamples) {
    eml_loop_ub = 0;
  } else {
    eml_loop_ub = (int32_T)eml_numSamples;
  }
  eml_m_hoistedExpr_sizes[0] = eml_loop_ub;
  eml_loop_ub--;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_j_hoistedExpr_data[eml_b_loop_ub] = 1.0 + (real_T)eml_b_loop_ub;
  }
  eml_b_hoistedExpr_sizes[0] = 1;
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_m_hoistedExpr_sizes[0];
  eml_b_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
  eml_loop_ub = eml_m_hoistedExpr_sizes[0] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_b_hoistedExpr_data[eml_b_loop_ub] =
      eml_j_hoistedExpr_data[eml_b_loop_ub];
  }
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_b_hoistedExpr_sizes[1];
  eml_b_frameBufferL_sizes[0] = 1;
  eml_b_frameBufferL_sizes[1] = eml_tmp_sizes[1];
  eml_loop_ub = eml_tmp_sizes[1] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    for(eml_y = 0; eml_y <= 0; eml_y = 1) {
      eml_yL_data[eml_b_loop_ub] =
        eml_ARampR[(int32_T)eml_b_hoistedExpr_data[eml_b_loop_ub] - 1];
    }
  }
  eml_b_yL_sizes[0] = 1;
  eml_b_yL_sizes[1] = eml_b_frameBufferL_sizes[1];
  eml_loop_ub = eml_b_frameBufferL_sizes[1] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_b_yL_data[eml_b_loop_ub] = (eml_yL_data[eml_b_loop_ub] < 1.0);
  }
  m_eml_li_find(eml_b_yL_data, eml_b_yL_sizes, eml_tmp_data,
   eml_frameBufferL_sizes);
  eml_b_tmp_sizes[0] = 1;
  eml_tmp_sizes[0] = 1;
  eml_tmp_sizes[1] = eml_frameBufferL_sizes[1];
  eml_b_tmp_sizes[1] = eml_tmp_sizes[1];
  eml_loop_ub = eml_frameBufferL_sizes[1] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    eml_b_tmp_data[eml_b_loop_ub] = eml_tmp_data[eml_b_loop_ub];
  }
  eml_loop_ub = eml_b_tmp_sizes[1] - 1;
  for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
    for(eml_y = 0; eml_y <= 0; eml_y = 1) {
      eml_ARampR[eml_b_tmp_data[eml_b_loop_ub] - 1] = 1.0;
    }
  }
  /* % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  DRNL */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  for(eml_filterCount = 1; (real_T)eml_filterCount <= eml_numChannels;
   eml_filterCount++) {
    eml_yR_sizes[0] = 1;
    eml_yR_sizes[1] = eml_stapesVelL_sizes[1];
    eml_loop_ub = eml_stapesVelL_sizes[1] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_yR_data[eml_b_loop_ub] = eml_stapesVelL_data[eml_b_loop_ub];
    }
    eml_b_yR_sizes[0] = 1;
    eml_b_yR_sizes[1] = eml_stapesVelR_sizes[1];
    eml_loop_ub = eml_stapesVelR_sizes[1] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_c_yR_data[eml_b_loop_ub] = eml_stapesVelR_data[eml_b_loop_ub];
    }
    if(1.0 > eml_numSamples) {
      eml_loop_ub = 0;
    } else {
      eml_loop_ub = (int32_T)eml_numSamples;
    }
    eml_n_hoistedExpr_sizes[0] = eml_loop_ub;
    eml_loop_ub--;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_k_hoistedExpr_data[eml_b_loop_ub] = 1.0 + (real_T)eml_b_loop_ub;
    }
    eml_b_hoistedExpr_sizes[0] = 1;
    eml_tmp_sizes[0] = 1;
    eml_tmp_sizes[1] = eml_n_hoistedExpr_sizes[0];
    eml_b_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
    eml_loop_ub = eml_n_hoistedExpr_sizes[0] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_b_hoistedExpr_data[eml_b_loop_ub] =
        eml_k_hoistedExpr_data[eml_b_loop_ub];
    }
    eml_hoistedExpr_sizes[0] = eml_b_hoistedExpr_sizes[1];
    eml_loop_ub = eml_filterCount - 1;
    eml_b_loop_ub = eml_hoistedExpr_sizes[0] - 1;
    for(eml_y = 0; eml_y <= eml_b_loop_ub; eml_y++) {
      eml_yL_data[eml_y] = eml_MOCcontrol[eml_loop_ub + 11 *
        ((int32_T)eml_b_hoistedExpr_data[eml_y] - 1)];
    }
    /* syntactic shorthand */
    if(eml_useGTF) {
      for(eml_c_loop_ub = 1.0; eml_c_loop_ub <= eml_filterOrder;
       eml_c_loop_ub++) {
        eml_b_MOCthreshold_dBOP = ((real_T)(17 * (eml_filterCount - 1)) + 3.0) +
          (eml_c_loop_ub - 1.0) * 2.0;
        eml_d0 = ((real_T)(17 * (eml_filterCount - 1)) + 4.0) + (eml_c_loop_ub -
          1.0) * 2.0;
        if(eml_b_MOCthreshold_dBOP > eml_d0) {
          eml_b_MOCthreshold_dBOP = 1.0;
          eml_d0 = 0.0;
        }
        eml_o_hoistedExpr_sizes[0] = ((int32_T)eml_d0 -
          (int32_T)eml_b_MOCthreshold_dBOP) + 1;
        eml_loop_ub = (int32_T)eml_d0 - (int32_T)eml_b_MOCthreshold_dBOP;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_l_hoistedExpr_data[eml_b_loop_ub] = eml_b_MOCthreshold_dBOP +
            (real_T)eml_b_loop_ub;
        }
        eml_c_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_o_hoistedExpr_sizes[0];
        eml_c_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_o_hoistedExpr_sizes[0] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_c_tmp_data[eml_b_loop_ub] = eml_l_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_loop_ub = 10 * (eml_filterCount - 1);
        eml_y = 10 * (eml_filterCount - 1);
        for(eml_b_loop_ub = 0; eml_b_loop_ub < 2; eml_b_loop_ub++) {
          eml_b_filterCoeffs[eml_b_loop_ub] = eml_filterCoeffs[(eml_loop_ub + (9
            + eml_b_loop_ub)) - 1];
        }
        for(eml_loop_ub = 0; eml_loop_ub < 3; eml_loop_ub++) {
          eml_f_filterCoeffs[eml_loop_ub] = eml_filterCoeffs[(eml_y + (14 +
            eml_loop_ub)) - 1];
        }
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_c_tmp_sizes[1];
        eml_p_hoistedExpr_sizes[0] = 1;
        eml_p_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          for(eml_y = 0; eml_y <= 0; eml_y = 1) {
            eml_m_hoistedExpr_data[eml_b_loop_ub] =
              eml_filterStatesL[(int32_T)eml_c_tmp_data[eml_b_loop_ub] - 1];
          }
        }
        eml_filterStatesL_sizes[0] = eml_p_hoistedExpr_sizes[1];
        eml_loop_ub = eml_p_hoistedExpr_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_filterStatesL_data[eml_b_loop_ub] =
            eml_m_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_c_yR_sizes[0] = 1;
        eml_c_yR_sizes[1] = eml_yR_sizes[1];
        eml_loop_ub = eml_yR_sizes[0] * eml_yR_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_b_yR_data[eml_b_loop_ub] = eml_yR_data[eml_b_loop_ub];
        }
        m_b_filter(eml_b_filterCoeffs, eml_f_filterCoeffs, eml_b_yR_data,
         eml_c_yR_sizes, eml_filterStatesL_data,
         eml_filterStatesL_sizes, eml_yR_data, eml_yR_sizes, eml_c_filterCoeffs);
        eml_b_MOCthreshold_dBOP = ((real_T)(17 * (eml_filterCount - 1)) + 3.0) +
          (eml_c_loop_ub - 1.0) * 2.0;
        eml_d0 = ((real_T)(17 * (eml_filterCount - 1)) + 4.0) + (eml_c_loop_ub -
          1.0) * 2.0;
        if(eml_b_MOCthreshold_dBOP > eml_d0) {
          eml_b_MOCthreshold_dBOP = 1.0;
          eml_d0 = 0.0;
        }
        eml_q_hoistedExpr_sizes[0] = ((int32_T)eml_d0 -
          (int32_T)eml_b_MOCthreshold_dBOP) + 1;
        eml_loop_ub = (int32_T)eml_d0 - (int32_T)eml_b_MOCthreshold_dBOP;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_n_hoistedExpr_data[eml_b_loop_ub] = eml_b_MOCthreshold_dBOP +
            (real_T)eml_b_loop_ub;
        }
        eml_c_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_q_hoistedExpr_sizes[0];
        eml_c_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_q_hoistedExpr_sizes[0] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_c_tmp_data[eml_b_loop_ub] = eml_n_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_d_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_c_tmp_sizes[1];
        eml_d_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_c_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_d_tmp_data[eml_b_loop_ub] = (int32_T)eml_c_tmp_data[eml_b_loop_ub];
        }
        eml_loop_ub = eml_d_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          for(eml_y = 0; eml_y <= 0; eml_y = 1) {
            eml_filterStatesL[eml_d_tmp_data[eml_b_loop_ub] - 1] =
              eml_c_filterCoeffs[eml_b_loop_ub];
          }
        }
        eml_b_MOCthreshold_dBOP = ((real_T)(17 * (eml_filterCount - 1)) + 3.0) +
          (eml_c_loop_ub - 1.0) * 2.0;
        eml_d0 = ((real_T)(17 * (eml_filterCount - 1)) + 4.0) + (eml_c_loop_ub -
          1.0) * 2.0;
        if(eml_b_MOCthreshold_dBOP > eml_d0) {
          eml_b_MOCthreshold_dBOP = 1.0;
          eml_d0 = 0.0;
        }
        eml_r_hoistedExpr_sizes[0] = ((int32_T)eml_d0 -
          (int32_T)eml_b_MOCthreshold_dBOP) + 1;
        eml_loop_ub = (int32_T)eml_d0 - (int32_T)eml_b_MOCthreshold_dBOP;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_o_hoistedExpr_data[eml_b_loop_ub] = eml_b_MOCthreshold_dBOP +
            (real_T)eml_b_loop_ub;
        }
        eml_c_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_r_hoistedExpr_sizes[0];
        eml_c_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_r_hoistedExpr_sizes[0] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_c_tmp_data[eml_b_loop_ub] = eml_o_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_loop_ub = 10 * (eml_filterCount - 1);
        eml_y = 10 * (eml_filterCount - 1);
        for(eml_b_loop_ub = 0; eml_b_loop_ub < 2; eml_b_loop_ub++) {
          eml_b_filterCoeffs[eml_b_loop_ub] = eml_filterCoeffs[(eml_loop_ub + (9
            + eml_b_loop_ub)) - 1];
        }
        for(eml_loop_ub = 0; eml_loop_ub < 3; eml_loop_ub++) {
          eml_f_filterCoeffs[eml_loop_ub] = eml_filterCoeffs[(eml_y + (14 +
            eml_loop_ub)) - 1];
        }
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_c_tmp_sizes[1];
        eml_s_hoistedExpr_sizes[0] = 1;
        eml_s_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          for(eml_y = 0; eml_y <= 0; eml_y = 1) {
            eml_p_hoistedExpr_data[eml_b_loop_ub] =
              eml_filterStatesR[(int32_T)eml_c_tmp_data[eml_b_loop_ub] - 1];
          }
        }
        eml_filterStatesR_sizes[0] = eml_s_hoistedExpr_sizes[1];
        eml_loop_ub = eml_s_hoistedExpr_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_filterStatesR_data[eml_b_loop_ub] =
            eml_p_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_d_yR_sizes[0] = 1;
        eml_d_yR_sizes[1] = eml_b_yR_sizes[1];
        eml_loop_ub = eml_b_yR_sizes[0] * eml_b_yR_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_b_yR_data[eml_b_loop_ub] = eml_c_yR_data[eml_b_loop_ub];
        }
        m_b_filter(eml_b_filterCoeffs, eml_f_filterCoeffs, eml_b_yR_data,
         eml_d_yR_sizes, eml_filterStatesR_data,
         eml_filterStatesR_sizes, eml_c_yR_data, eml_b_yR_sizes,
         eml_c_filterCoeffs);
        eml_b_MOCthreshold_dBOP = ((real_T)(17 * (eml_filterCount - 1)) + 3.0) +
          (eml_c_loop_ub - 1.0) * 2.0;
        eml_d0 = ((real_T)(17 * (eml_filterCount - 1)) + 4.0) + (eml_c_loop_ub -
          1.0) * 2.0;
        if(eml_b_MOCthreshold_dBOP > eml_d0) {
          eml_b_MOCthreshold_dBOP = 1.0;
          eml_d0 = 0.0;
        }
        eml_t_hoistedExpr_sizes[0] = ((int32_T)eml_d0 -
          (int32_T)eml_b_MOCthreshold_dBOP) + 1;
        eml_loop_ub = (int32_T)eml_d0 - (int32_T)eml_b_MOCthreshold_dBOP;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_q_hoistedExpr_data[eml_b_loop_ub] = eml_b_MOCthreshold_dBOP +
            (real_T)eml_b_loop_ub;
        }
        eml_c_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_t_hoistedExpr_sizes[0];
        eml_c_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_t_hoistedExpr_sizes[0] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_c_tmp_data[eml_b_loop_ub] = eml_q_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_d_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_c_tmp_sizes[1];
        eml_d_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_c_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_d_tmp_data[eml_b_loop_ub] = (int32_T)eml_c_tmp_data[eml_b_loop_ub];
        }
        eml_loop_ub = eml_d_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          for(eml_y = 0; eml_y <= 0; eml_y = 1) {
            eml_filterStatesR[eml_d_tmp_data[eml_b_loop_ub] - 1] =
              eml_c_filterCoeffs[eml_b_loop_ub];
          }
        }
      }
    } else {
      eml_c_loop_ub = eml_filterOrder / 2.0;
      for(eml_nn = 1.0; eml_nn <= eml_c_loop_ub; eml_nn++) {
        eml_b_MOCthreshold_dBOP = ((real_T)(17 * (eml_filterCount - 1)) + 3.0) +
          (eml_nn - 1.0) * 4.0;
        eml_d0 = ((real_T)(17 * (eml_filterCount - 1)) + 6.0) + (eml_nn - 1.0) *
          4.0;
        if(eml_b_MOCthreshold_dBOP > eml_d0) {
          eml_b_MOCthreshold_dBOP = 1.0;
          eml_d0 = 0.0;
        }
        eml_u_hoistedExpr_sizes[0] = ((int32_T)eml_d0 -
          (int32_T)eml_b_MOCthreshold_dBOP) + 1;
        eml_loop_ub = (int32_T)eml_d0 - (int32_T)eml_b_MOCthreshold_dBOP;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_r_hoistedExpr_data[eml_b_loop_ub] = eml_b_MOCthreshold_dBOP +
            (real_T)eml_b_loop_ub;
        }
        eml_e_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_u_hoistedExpr_sizes[0];
        eml_e_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_u_hoistedExpr_sizes[0] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_e_tmp_data[eml_b_loop_ub] = eml_r_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_loop_ub = 10 * (eml_filterCount - 1);
        eml_b_loop_ub = 10 * (eml_filterCount - 1);
        for(eml_y = 0; eml_y < 5; eml_y++) {
          eml_g_filterCoeffs[eml_y] = eml_filterCoeffs[(eml_loop_ub + (9 +
            eml_y)) - 1];
          eml_h_filterCoeffs[eml_y] = eml_filterCoeffs[(eml_b_loop_ub + (14 +
            eml_y)) - 1];
        }
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_e_tmp_sizes[1];
        eml_v_hoistedExpr_sizes[0] = 1;
        eml_v_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          for(eml_y = 0; eml_y <= 0; eml_y = 1) {
            eml_s_hoistedExpr_data[eml_b_loop_ub] =
              eml_filterStatesL[(int32_T)eml_e_tmp_data[eml_b_loop_ub] - 1];
          }
        }
        eml_b_filterStatesL_sizes[0] = eml_v_hoistedExpr_sizes[1];
        eml_loop_ub = eml_v_hoistedExpr_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_b_filterStatesL_data[eml_b_loop_ub] =
            eml_s_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_e_yR_sizes[0] = 1;
        eml_e_yR_sizes[1] = eml_yR_sizes[1];
        eml_loop_ub = eml_yR_sizes[0] * eml_yR_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_b_yR_data[eml_b_loop_ub] = eml_yR_data[eml_b_loop_ub];
        }
        m_c_filter(eml_g_filterCoeffs, eml_h_filterCoeffs, eml_b_yR_data,
         eml_e_yR_sizes, eml_b_filterStatesL_data,
         eml_b_filterStatesL_sizes, eml_yR_data, eml_yR_sizes, eml_dv0);
        eml_b_MOCthreshold_dBOP = ((real_T)(17 * (eml_filterCount - 1)) + 3.0) +
          (eml_nn - 1.0) * 4.0;
        eml_d0 = ((real_T)(17 * (eml_filterCount - 1)) + 6.0) + (eml_nn - 1.0) *
          4.0;
        if(eml_b_MOCthreshold_dBOP > eml_d0) {
          eml_b_MOCthreshold_dBOP = 1.0;
          eml_d0 = 0.0;
        }
        eml_w_hoistedExpr_sizes[0] = ((int32_T)eml_d0 -
          (int32_T)eml_b_MOCthreshold_dBOP) + 1;
        eml_loop_ub = (int32_T)eml_d0 - (int32_T)eml_b_MOCthreshold_dBOP;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_t_hoistedExpr_data[eml_b_loop_ub] = eml_b_MOCthreshold_dBOP +
            (real_T)eml_b_loop_ub;
        }
        eml_e_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_w_hoistedExpr_sizes[0];
        eml_e_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_w_hoistedExpr_sizes[0] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_e_tmp_data[eml_b_loop_ub] = eml_t_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_f_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_e_tmp_sizes[1];
        eml_f_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_e_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_f_tmp_data[eml_b_loop_ub] = (int32_T)eml_e_tmp_data[eml_b_loop_ub];
        }
        eml_loop_ub = eml_f_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          for(eml_y = 0; eml_y <= 0; eml_y = 1) {
            eml_filterStatesL[eml_f_tmp_data[eml_b_loop_ub] - 1] =
              eml_dv0[eml_b_loop_ub];
          }
        }
        eml_b_MOCthreshold_dBOP = ((real_T)(17 * (eml_filterCount - 1)) + 3.0) +
          (eml_nn - 1.0) * 4.0;
        eml_d0 = ((real_T)(17 * (eml_filterCount - 1)) + 6.0) + (eml_nn - 1.0) *
          4.0;
        if(eml_b_MOCthreshold_dBOP > eml_d0) {
          eml_b_MOCthreshold_dBOP = 1.0;
          eml_d0 = 0.0;
        }
        eml_x_hoistedExpr_sizes[0] = ((int32_T)eml_d0 -
          (int32_T)eml_b_MOCthreshold_dBOP) + 1;
        eml_loop_ub = (int32_T)eml_d0 - (int32_T)eml_b_MOCthreshold_dBOP;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_u_hoistedExpr_data[eml_b_loop_ub] = eml_b_MOCthreshold_dBOP +
            (real_T)eml_b_loop_ub;
        }
        eml_e_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_x_hoistedExpr_sizes[0];
        eml_e_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_x_hoistedExpr_sizes[0] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_e_tmp_data[eml_b_loop_ub] = eml_u_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_loop_ub = 10 * (eml_filterCount - 1);
        eml_b_loop_ub = 10 * (eml_filterCount - 1);
        for(eml_y = 0; eml_y < 5; eml_y++) {
          eml_g_filterCoeffs[eml_y] = eml_filterCoeffs[(eml_loop_ub + (9 +
            eml_y)) - 1];
          eml_h_filterCoeffs[eml_y] = eml_filterCoeffs[(eml_b_loop_ub + (14 +
            eml_y)) - 1];
        }
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_e_tmp_sizes[1];
        eml_y_hoistedExpr_sizes[0] = 1;
        eml_y_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          for(eml_y = 0; eml_y <= 0; eml_y = 1) {
            eml_v_hoistedExpr_data[eml_b_loop_ub] =
              eml_filterStatesR[(int32_T)eml_e_tmp_data[eml_b_loop_ub] - 1];
          }
        }
        eml_b_filterStatesR_sizes[0] = eml_y_hoistedExpr_sizes[1];
        eml_loop_ub = eml_y_hoistedExpr_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_b_filterStatesR_data[eml_b_loop_ub] =
            eml_v_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_f_yR_sizes[0] = 1;
        eml_f_yR_sizes[1] = eml_b_yR_sizes[1];
        eml_loop_ub = eml_b_yR_sizes[0] * eml_b_yR_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_b_yR_data[eml_b_loop_ub] = eml_c_yR_data[eml_b_loop_ub];
        }
        m_c_filter(eml_g_filterCoeffs, eml_h_filterCoeffs, eml_b_yR_data,
         eml_f_yR_sizes, eml_b_filterStatesR_data,
         eml_b_filterStatesR_sizes, eml_c_yR_data, eml_b_yR_sizes, eml_dv0);
        eml_b_MOCthreshold_dBOP = ((real_T)(17 * (eml_filterCount - 1)) + 3.0) +
          (eml_nn - 1.0) * 4.0;
        eml_d0 = ((real_T)(17 * (eml_filterCount - 1)) + 6.0) + (eml_nn - 1.0) *
          4.0;
        if(eml_b_MOCthreshold_dBOP > eml_d0) {
          eml_b_MOCthreshold_dBOP = 1.0;
          eml_d0 = 0.0;
        }
        eml_ab_hoistedExpr_sizes[0] = ((int32_T)eml_d0 -
          (int32_T)eml_b_MOCthreshold_dBOP) + 1;
        eml_loop_ub = (int32_T)eml_d0 - (int32_T)eml_b_MOCthreshold_dBOP;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_w_hoistedExpr_data[eml_b_loop_ub] = eml_b_MOCthreshold_dBOP +
            (real_T)eml_b_loop_ub;
        }
        eml_e_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_ab_hoistedExpr_sizes[0];
        eml_e_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_ab_hoistedExpr_sizes[0] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_e_tmp_data[eml_b_loop_ub] = eml_w_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_f_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_e_tmp_sizes[1];
        eml_f_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_e_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_f_tmp_data[eml_b_loop_ub] = (int32_T)eml_e_tmp_data[eml_b_loop_ub];
        }
        eml_loop_ub = eml_f_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          for(eml_y = 0; eml_y <= 0; eml_y = 1) {
            eml_filterStatesR[eml_f_tmp_data[eml_b_loop_ub] - 1] =
              eml_dv0[eml_b_loop_ub];
          }
        }
      }
    }
    m_abs(eml_yR_data, eml_yR_sizes, eml_b_hoistedExpr_data,
     eml_b_hoistedExpr_sizes);
    eml_peakIPL[eml_filterCount - 1] = 20.0 *
    log10(m_max(eml_b_hoistedExpr_data, eml_b_hoistedExpr_sizes) / 0.00002);
    /* peak in in dB SPL */
    m_abs(eml_c_yR_data, eml_b_yR_sizes, eml_b_hoistedExpr_data,
     eml_b_hoistedExpr_sizes);
    eml_peakIPR[eml_filterCount - 1] = 20.0 *
    log10(m_max(eml_b_hoistedExpr_data, eml_b_hoistedExpr_sizes) / 0.00002);
    /* peak in in dB SPL */
    /*      rmsIPL(filterCount) = 20*log10(  sqrt(mean(yL.^2)) /2e-5 ); %rms in in dB SPL */
    /*      rmsIPR(filterCount) = 20*log10(  sqrt(mean(yR.^2)) /2e-5 ); %rms in in dB SPL */
    eml_rmsIPL[eml_filterCount - 1] = 20.0 *
      log10(fabs(eml_yR_data[eml_yR_sizes[1] - 1]) / 0.00002);
    /* rms in in dB SPL for GUI (bit of a hack, but it is smoothed in GUI) */
    eml_rmsIPR[eml_filterCount - 1] = 20.0 *
      log10(fabs(eml_c_yR_data[eml_b_yR_sizes[1] - 1]) / 0.00002);
    /* rms in in dB SPL   */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /*  compression time */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /*  MOC compression first */
    eml_yR_sizes[0] = 1;
    eml_loop_ub = eml_yR_sizes[1];
    eml_loop_ub--;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_yR_data[eml_b_loop_ub] *= eml_yL_data[eml_b_loop_ub];
    }
    eml_b_yR_sizes[0] = 1;
    eml_loop_ub = eml_b_yR_sizes[1];
    eml_loop_ub--;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_c_yR_data[eml_b_loop_ub] *= eml_yL_data[eml_b_loop_ub];
    }
    /*  Instantaneous after */
    m_refp1_DRNL_brokenstick_nl(eml_yR_data, eml_yR_sizes,
     eml_DRNLb[eml_filterCount - 1], eml_DRNLc[eml_filterCount - 1]);
    m_refp1_DRNL_brokenstick_nl(eml_c_yR_data, eml_b_yR_sizes,
     eml_DRNLb[eml_filterCount - 1], eml_DRNLc[eml_filterCount - 1]);
    if(eml_useGTF) {
      for(eml_c_loop_ub = 1.0; eml_c_loop_ub <= eml_filterOrder;
       eml_c_loop_ub++) {
        eml_b_MOCthreshold_dBOP = ((real_T)(17 * (eml_filterCount - 1)) + 11.0)
          + (eml_c_loop_ub - 1.0) * 2.0;
        eml_d0 = ((real_T)(17 * (eml_filterCount - 1)) + 12.0) + (eml_c_loop_ub
          - 1.0) * 2.0;
        if(eml_b_MOCthreshold_dBOP > eml_d0) {
          eml_b_MOCthreshold_dBOP = 1.0;
          eml_d0 = 0.0;
        }
        eml_bb_hoistedExpr_sizes[0] = ((int32_T)eml_d0 -
          (int32_T)eml_b_MOCthreshold_dBOP) + 1;
        eml_loop_ub = (int32_T)eml_d0 - (int32_T)eml_b_MOCthreshold_dBOP;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_x_hoistedExpr_data[eml_b_loop_ub] = eml_b_MOCthreshold_dBOP +
            (real_T)eml_b_loop_ub;
        }
        eml_g_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_bb_hoistedExpr_sizes[0];
        eml_g_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_bb_hoistedExpr_sizes[0] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_g_tmp_data[eml_b_loop_ub] = eml_x_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_loop_ub = 10 * (eml_filterCount - 1);
        eml_y = 10 * (eml_filterCount - 1);
        for(eml_b_loop_ub = 0; eml_b_loop_ub < 2; eml_b_loop_ub++) {
          eml_b_filterCoeffs[eml_b_loop_ub] = eml_filterCoeffs[(eml_loop_ub + (9
            + eml_b_loop_ub)) - 1];
        }
        for(eml_loop_ub = 0; eml_loop_ub < 3; eml_loop_ub++) {
          eml_f_filterCoeffs[eml_loop_ub] = eml_filterCoeffs[(eml_y + (14 +
            eml_loop_ub)) - 1];
        }
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_g_tmp_sizes[1];
        eml_cb_hoistedExpr_sizes[0] = 1;
        eml_cb_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          for(eml_y = 0; eml_y <= 0; eml_y = 1) {
            eml_y_hoistedExpr_data[eml_b_loop_ub] =
              eml_filterStatesL[(int32_T)eml_g_tmp_data[eml_b_loop_ub] - 1];
          }
        }
        eml_c_filterStatesL_sizes[0] = eml_cb_hoistedExpr_sizes[1];
        eml_loop_ub = eml_cb_hoistedExpr_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_c_filterStatesL_data[eml_b_loop_ub] =
            eml_y_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_g_yR_sizes[0] = 1;
        eml_g_yR_sizes[1] = eml_yR_sizes[1];
        eml_loop_ub = eml_yR_sizes[0] * eml_yR_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_b_yR_data[eml_b_loop_ub] = eml_yR_data[eml_b_loop_ub];
        }
        m_b_filter(eml_b_filterCoeffs, eml_f_filterCoeffs, eml_b_yR_data,
         eml_g_yR_sizes, eml_c_filterStatesL_data,
         eml_c_filterStatesL_sizes, eml_yR_data, eml_yR_sizes,
         eml_c_filterCoeffs);
        eml_b_MOCthreshold_dBOP = ((real_T)(17 * (eml_filterCount - 1)) + 11.0)
          + (eml_c_loop_ub - 1.0) * 2.0;
        eml_d0 = ((real_T)(17 * (eml_filterCount - 1)) + 12.0) + (eml_c_loop_ub
          - 1.0) * 2.0;
        if(eml_b_MOCthreshold_dBOP > eml_d0) {
          eml_b_MOCthreshold_dBOP = 1.0;
          eml_d0 = 0.0;
        }
        eml_db_hoistedExpr_sizes[0] = ((int32_T)eml_d0 -
          (int32_T)eml_b_MOCthreshold_dBOP) + 1;
        eml_loop_ub = (int32_T)eml_d0 - (int32_T)eml_b_MOCthreshold_dBOP;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_ab_hoistedExpr_data[eml_b_loop_ub] = eml_b_MOCthreshold_dBOP +
            (real_T)eml_b_loop_ub;
        }
        eml_g_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_db_hoistedExpr_sizes[0];
        eml_g_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_db_hoistedExpr_sizes[0] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_g_tmp_data[eml_b_loop_ub] = eml_ab_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_h_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_g_tmp_sizes[1];
        eml_h_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_g_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_h_tmp_data[eml_b_loop_ub] = (int32_T)eml_g_tmp_data[eml_b_loop_ub];
        }
        eml_loop_ub = eml_h_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          for(eml_y = 0; eml_y <= 0; eml_y = 1) {
            eml_filterStatesL[eml_h_tmp_data[eml_b_loop_ub] - 1] =
              eml_c_filterCoeffs[eml_b_loop_ub];
          }
        }
        eml_b_MOCthreshold_dBOP = ((real_T)(17 * (eml_filterCount - 1)) + 11.0)
          + (eml_c_loop_ub - 1.0) * 2.0;
        eml_d0 = ((real_T)(17 * (eml_filterCount - 1)) + 12.0) + (eml_c_loop_ub
          - 1.0) * 2.0;
        if(eml_b_MOCthreshold_dBOP > eml_d0) {
          eml_b_MOCthreshold_dBOP = 1.0;
          eml_d0 = 0.0;
        }
        eml_eb_hoistedExpr_sizes[0] = ((int32_T)eml_d0 -
          (int32_T)eml_b_MOCthreshold_dBOP) + 1;
        eml_loop_ub = (int32_T)eml_d0 - (int32_T)eml_b_MOCthreshold_dBOP;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_bb_hoistedExpr_data[eml_b_loop_ub] = eml_b_MOCthreshold_dBOP +
            (real_T)eml_b_loop_ub;
        }
        eml_g_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_eb_hoistedExpr_sizes[0];
        eml_g_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_eb_hoistedExpr_sizes[0] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_g_tmp_data[eml_b_loop_ub] = eml_bb_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_loop_ub = 10 * (eml_filterCount - 1);
        eml_y = 10 * (eml_filterCount - 1);
        for(eml_b_loop_ub = 0; eml_b_loop_ub < 2; eml_b_loop_ub++) {
          eml_b_filterCoeffs[eml_b_loop_ub] = eml_filterCoeffs[(eml_loop_ub + (9
            + eml_b_loop_ub)) - 1];
        }
        for(eml_loop_ub = 0; eml_loop_ub < 3; eml_loop_ub++) {
          eml_f_filterCoeffs[eml_loop_ub] = eml_filterCoeffs[(eml_y + (14 +
            eml_loop_ub)) - 1];
        }
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_g_tmp_sizes[1];
        eml_fb_hoistedExpr_sizes[0] = 1;
        eml_fb_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          for(eml_y = 0; eml_y <= 0; eml_y = 1) {
            eml_cb_hoistedExpr_data[eml_b_loop_ub] =
              eml_filterStatesR[(int32_T)eml_g_tmp_data[eml_b_loop_ub] - 1];
          }
        }
        eml_c_filterStatesR_sizes[0] = eml_fb_hoistedExpr_sizes[1];
        eml_loop_ub = eml_fb_hoistedExpr_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_c_filterStatesR_data[eml_b_loop_ub] =
            eml_cb_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_h_yR_sizes[0] = 1;
        eml_h_yR_sizes[1] = eml_b_yR_sizes[1];
        eml_loop_ub = eml_b_yR_sizes[0] * eml_b_yR_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_b_yR_data[eml_b_loop_ub] = eml_c_yR_data[eml_b_loop_ub];
        }
        m_b_filter(eml_b_filterCoeffs, eml_f_filterCoeffs, eml_b_yR_data,
         eml_h_yR_sizes, eml_c_filterStatesR_data,
         eml_c_filterStatesR_sizes, eml_c_yR_data, eml_b_yR_sizes,
         eml_c_filterCoeffs);
        eml_b_MOCthreshold_dBOP = ((real_T)(17 * (eml_filterCount - 1)) + 11.0)
          + (eml_c_loop_ub - 1.0) * 2.0;
        eml_d0 = ((real_T)(17 * (eml_filterCount - 1)) + 12.0) + (eml_c_loop_ub
          - 1.0) * 2.0;
        if(eml_b_MOCthreshold_dBOP > eml_d0) {
          eml_b_MOCthreshold_dBOP = 1.0;
          eml_d0 = 0.0;
        }
        eml_gb_hoistedExpr_sizes[0] = ((int32_T)eml_d0 -
          (int32_T)eml_b_MOCthreshold_dBOP) + 1;
        eml_loop_ub = (int32_T)eml_d0 - (int32_T)eml_b_MOCthreshold_dBOP;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_db_hoistedExpr_data[eml_b_loop_ub] = eml_b_MOCthreshold_dBOP +
            (real_T)eml_b_loop_ub;
        }
        eml_g_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_gb_hoistedExpr_sizes[0];
        eml_g_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_gb_hoistedExpr_sizes[0] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_g_tmp_data[eml_b_loop_ub] = eml_db_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_h_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_g_tmp_sizes[1];
        eml_h_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_g_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_h_tmp_data[eml_b_loop_ub] = (int32_T)eml_g_tmp_data[eml_b_loop_ub];
        }
        eml_loop_ub = eml_h_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          for(eml_y = 0; eml_y <= 0; eml_y = 1) {
            eml_filterStatesR[eml_h_tmp_data[eml_b_loop_ub] - 1] =
              eml_c_filterCoeffs[eml_b_loop_ub];
          }
        }
      }
    } else {
      eml_c_loop_ub = eml_filterOrder / 2.0;
      for(eml_nn = 1.0; eml_nn <= eml_c_loop_ub; eml_nn++) {
        eml_b_MOCthreshold_dBOP = ((real_T)(17 * (eml_filterCount - 1)) + 11.0)
          + (eml_nn - 1.0) * 4.0;
        eml_d0 = ((real_T)(17 * (eml_filterCount - 1)) + 14.0) + (eml_nn - 1.0)
          * 4.0;
        if(eml_b_MOCthreshold_dBOP > eml_d0) {
          eml_y = 1;
          eml_loop_ub = 0;
        } else {
          eml_y = (int32_T)eml_b_MOCthreshold_dBOP;
          eml_loop_ub = (int32_T)eml_d0;
        }
        eml_hb_hoistedExpr_sizes[0] = (eml_loop_ub - eml_y) + 1;
        eml_loop_ub -= eml_y;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_eb_hoistedExpr_data[eml_b_loop_ub] = (real_T)eml_y +
            (real_T)eml_b_loop_ub;
        }
        eml_i_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_hb_hoistedExpr_sizes[0];
        eml_i_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_hb_hoistedExpr_sizes[0] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_i_tmp_data[eml_b_loop_ub] = eml_eb_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_loop_ub = 10 * (eml_filterCount - 1);
        eml_b_loop_ub = 10 * (eml_filterCount - 1);
        for(eml_y = 0; eml_y < 5; eml_y++) {
          eml_g_filterCoeffs[eml_y] = eml_filterCoeffs[(eml_loop_ub + (9 +
            eml_y)) - 1];
          eml_h_filterCoeffs[eml_y] = eml_filterCoeffs[(eml_b_loop_ub + (14 +
            eml_y)) - 1];
        }
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_i_tmp_sizes[1];
        eml_ib_hoistedExpr_sizes[0] = 1;
        eml_ib_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          for(eml_y = 0; eml_y <= 0; eml_y = 1) {
            eml_fb_hoistedExpr_data[eml_b_loop_ub] =
              eml_filterStatesL[(int32_T)eml_i_tmp_data[eml_b_loop_ub] - 1];
          }
        }
        eml_d_filterStatesL_sizes[0] = eml_ib_hoistedExpr_sizes[1];
        eml_loop_ub = eml_ib_hoistedExpr_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_d_filterStatesL_data[eml_b_loop_ub] =
            eml_fb_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_i_yR_sizes[0] = 1;
        eml_i_yR_sizes[1] = eml_yR_sizes[1];
        eml_loop_ub = eml_yR_sizes[0] * eml_yR_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_b_yR_data[eml_b_loop_ub] = eml_yR_data[eml_b_loop_ub];
        }
        m_c_filter(eml_g_filterCoeffs, eml_h_filterCoeffs, eml_b_yR_data,
         eml_i_yR_sizes, eml_d_filterStatesL_data,
         eml_d_filterStatesL_sizes, eml_yR_data, eml_yR_sizes, eml_dv0);
        eml_b_MOCthreshold_dBOP = ((real_T)(17 * (eml_filterCount - 1)) + 11.0)
          + (eml_nn - 1.0) * 4.0;
        eml_d0 = ((real_T)(17 * (eml_filterCount - 1)) + 14.0) + (eml_nn - 1.0)
          * 4.0;
        if(eml_b_MOCthreshold_dBOP > eml_d0) {
          eml_y = 1;
          eml_loop_ub = 0;
        } else {
          eml_y = (int32_T)eml_b_MOCthreshold_dBOP;
          eml_loop_ub = (int32_T)eml_d0;
        }
        eml_jb_hoistedExpr_sizes[0] = (eml_loop_ub - eml_y) + 1;
        eml_loop_ub -= eml_y;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_gb_hoistedExpr_data[eml_b_loop_ub] = (real_T)eml_y +
            (real_T)eml_b_loop_ub;
        }
        eml_i_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_jb_hoistedExpr_sizes[0];
        eml_i_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_jb_hoistedExpr_sizes[0] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_i_tmp_data[eml_b_loop_ub] = eml_gb_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_j_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_i_tmp_sizes[1];
        eml_j_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_i_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_j_tmp_data[eml_b_loop_ub] = (int32_T)eml_i_tmp_data[eml_b_loop_ub];
        }
        eml_loop_ub = eml_j_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          for(eml_y = 0; eml_y <= 0; eml_y = 1) {
            eml_filterStatesL[eml_j_tmp_data[eml_b_loop_ub] - 1] =
              eml_dv0[eml_b_loop_ub];
          }
        }
        eml_b_MOCthreshold_dBOP = ((real_T)(17 * (eml_filterCount - 1)) + 11.0)
          + (eml_nn - 1.0) * 4.0;
        eml_d0 = ((real_T)(17 * (eml_filterCount - 1)) + 14.0) + (eml_nn - 1.0)
          * 4.0;
        if(eml_b_MOCthreshold_dBOP > eml_d0) {
          eml_y = 1;
          eml_loop_ub = 0;
        } else {
          eml_y = (int32_T)eml_b_MOCthreshold_dBOP;
          eml_loop_ub = (int32_T)eml_d0;
        }
        eml_kb_hoistedExpr_sizes[0] = (eml_loop_ub - eml_y) + 1;
        eml_loop_ub -= eml_y;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_hb_hoistedExpr_data[eml_b_loop_ub] = (real_T)eml_y +
            (real_T)eml_b_loop_ub;
        }
        eml_i_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_kb_hoistedExpr_sizes[0];
        eml_i_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_kb_hoistedExpr_sizes[0] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_i_tmp_data[eml_b_loop_ub] = eml_hb_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_loop_ub = 10 * (eml_filterCount - 1);
        eml_b_loop_ub = 10 * (eml_filterCount - 1);
        for(eml_y = 0; eml_y < 5; eml_y++) {
          eml_g_filterCoeffs[eml_y] = eml_filterCoeffs[(eml_loop_ub + (9 +
            eml_y)) - 1];
          eml_h_filterCoeffs[eml_y] = eml_filterCoeffs[(eml_b_loop_ub + (14 +
            eml_y)) - 1];
        }
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_i_tmp_sizes[1];
        eml_lb_hoistedExpr_sizes[0] = 1;
        eml_lb_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          for(eml_y = 0; eml_y <= 0; eml_y = 1) {
            eml_ib_hoistedExpr_data[eml_b_loop_ub] =
              eml_filterStatesR[(int32_T)eml_i_tmp_data[eml_b_loop_ub] - 1];
          }
        }
        eml_d_filterStatesR_sizes[0] = eml_lb_hoistedExpr_sizes[1];
        eml_loop_ub = eml_lb_hoistedExpr_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_d_filterStatesR_data[eml_b_loop_ub] =
            eml_ib_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_j_yR_sizes[0] = 1;
        eml_j_yR_sizes[1] = eml_b_yR_sizes[1];
        eml_loop_ub = eml_b_yR_sizes[0] * eml_b_yR_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_b_yR_data[eml_b_loop_ub] = eml_c_yR_data[eml_b_loop_ub];
        }
        m_c_filter(eml_g_filterCoeffs, eml_h_filterCoeffs, eml_b_yR_data,
         eml_j_yR_sizes, eml_d_filterStatesR_data,
         eml_d_filterStatesR_sizes, eml_c_yR_data, eml_b_yR_sizes, eml_dv0);
        eml_b_MOCthreshold_dBOP = ((real_T)(17 * (eml_filterCount - 1)) + 11.0)
          + (eml_nn - 1.0) * 4.0;
        eml_d0 = ((real_T)(17 * (eml_filterCount - 1)) + 14.0) + (eml_nn - 1.0)
          * 4.0;
        if(eml_b_MOCthreshold_dBOP > eml_d0) {
          eml_y = 1;
          eml_loop_ub = 0;
        } else {
          eml_y = (int32_T)eml_b_MOCthreshold_dBOP;
          eml_loop_ub = (int32_T)eml_d0;
        }
        eml_mb_hoistedExpr_sizes[0] = (eml_loop_ub - eml_y) + 1;
        eml_loop_ub -= eml_y;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_jb_hoistedExpr_data[eml_b_loop_ub] = (real_T)eml_y +
            (real_T)eml_b_loop_ub;
        }
        eml_i_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_mb_hoistedExpr_sizes[0];
        eml_i_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_mb_hoistedExpr_sizes[0] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_i_tmp_data[eml_b_loop_ub] = eml_jb_hoistedExpr_data[eml_b_loop_ub];
        }
        eml_j_tmp_sizes[0] = 1;
        eml_tmp_sizes[0] = 1;
        eml_tmp_sizes[1] = eml_i_tmp_sizes[1];
        eml_j_tmp_sizes[1] = eml_tmp_sizes[1];
        eml_loop_ub = eml_i_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          eml_j_tmp_data[eml_b_loop_ub] = (int32_T)eml_i_tmp_data[eml_b_loop_ub];
        }
        eml_loop_ub = eml_j_tmp_sizes[1] - 1;
        for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
          for(eml_y = 0; eml_y <= 0; eml_y = 1) {
            eml_filterStatesR[eml_j_tmp_data[eml_b_loop_ub] - 1] =
              eml_dv0[eml_b_loop_ub];
          }
        }
      }
    }
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /*  Deal with MOC control signal for next frame */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    m_power(eml_yR_data, eml_yR_sizes, eml_b_hoistedExpr_data,
     eml_b_hoistedExpr_sizes);
    m_power(eml_c_yR_data, eml_b_yR_sizes, eml_yL_data,
     eml_b_frameBufferL_sizes);
    eml_nb_hoistedExpr_sizes[0] = 2;
    eml_nb_hoistedExpr_sizes[1] = eml_b_hoistedExpr_sizes[1];
    eml_loop_ub = eml_b_hoistedExpr_sizes[1] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_kb_hoistedExpr_data[eml_b_loop_ub << 1] =
        eml_b_hoistedExpr_data[eml_b_hoistedExpr_sizes[0] * eml_b_loop_ub];
    }
    eml_loop_ub = eml_b_frameBufferL_sizes[1] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_kb_hoistedExpr_data[1 + (eml_b_loop_ub << 1)] =
        eml_yL_data[eml_b_frameBufferL_sizes[0] * eml_b_loop_ub];
    }
    m_sum(eml_kb_hoistedExpr_data, eml_nb_hoistedExpr_sizes,
     eml_b_hoistedExpr_data, eml_b_hoistedExpr_sizes);
    for(eml_loop_ub = 0; eml_loop_ub < 2; eml_loop_ub++) {
      eml_b_filterCoeffs[eml_loop_ub] = eml_filterCoeffs[eml_loop_ub + 4];
      eml_c_filterCoeffs[eml_loop_ub] = eml_filterCoeffs[eml_loop_ub + 6];
    }
    eml_ob_hoistedExpr_sizes[0] = 1;
    eml_ob_hoistedExpr_sizes[1] = eml_b_hoistedExpr_sizes[1];
    eml_loop_ub = eml_b_hoistedExpr_sizes[0] * eml_b_hoistedExpr_sizes[1] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_b_yR_data[eml_b_loop_ub] = eml_b_hoistedExpr_data[eml_b_loop_ub] / 2.0;
    }
    m_filter(eml_b_filterCoeffs, eml_c_filterCoeffs, eml_b_yR_data,
     eml_ob_hoistedExpr_sizes, eml_filterStatesL[17 * (
      eml_filterCount - 1) + 1], eml_yL_data, eml_b_frameBufferL_sizes,
     &eml_b_MOCthreshold_dBOP);
    eml_filterStatesL[17 * (eml_filterCount - 1) + 1] = eml_b_MOCthreshold_dBOP;
    m_refp1_sqrt(eml_yL_data, eml_b_frameBufferL_sizes);
    /*  restore to meaningful scale (meters - not anymore now in velocity mode) */
    /*      disp(MOCthreshold_dBOP) */
    eml_b_frameBufferL_sizes[0] = 1;
    eml_loop_ub = eml_b_frameBufferL_sizes[1] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_yL_data[eml_b_loop_ub] /= 0.00002;
    }
    m_refp1_log10(eml_yL_data, eml_b_frameBufferL_sizes);
    eml_b_frameBufferL_sizes[0] = 1;
    eml_loop_ub = eml_b_frameBufferL_sizes[1] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_yL_data[eml_b_loop_ub] *= 20.0;
    }
    eml_c_yL_sizes[0] = 1;
    eml_c_yL_sizes[1] = eml_b_frameBufferL_sizes[1];
    eml_b_MOCthreshold_dBOP = eml_MOCthreshold_dBOP[eml_filterCount - 1];
    eml_loop_ub = eml_b_frameBufferL_sizes[1] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_b_yR_data[eml_b_loop_ub] = eml_yL_data[eml_b_loop_ub] -
        eml_b_MOCthreshold_dBOP;
    }
    m_b_max(eml_b_yR_data, eml_c_yL_sizes, eml_b_hoistedExpr_data,
     eml_b_hoistedExpr_sizes);
    eml_b_hoistedExpr_sizes[0] = 1;
    eml_loop_ub = eml_b_hoistedExpr_sizes[1] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_b_hoistedExpr_data[eml_b_loop_ub] *= eml_MOCfactor;
    }
    /* the tiny offset is due to a crummy Embedded Matlab bug (This caused much misery)! */
    eml_b_hoistedExpr_sizes[0] = 1;
    eml_loop_ub = eml_b_hoistedExpr_sizes[1] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_b_hoistedExpr_data[eml_b_loop_ub] =
        -(eml_b_hoistedExpr_data[eml_b_loop_ub] + 1.0E-009);
    }
    if(1.0 > eml_numSamples) {
      eml_loop_ub = 0;
    } else {
      eml_loop_ub = (int32_T)eml_numSamples;
    }
    eml_pb_hoistedExpr_sizes[0] = eml_loop_ub;
    eml_loop_ub--;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_lb_hoistedExpr_data[eml_b_loop_ub] = 1.0 + (real_T)eml_b_loop_ub;
    }
    eml_loop_ub = eml_pb_hoistedExpr_sizes[0] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_yL_data[eml_b_loop_ub] = eml_lb_hoistedExpr_data[eml_b_loop_ub];
    }
    eml_qb_hoistedExpr_sizes[0] = 1;
    eml_qb_hoistedExpr_sizes[1] = eml_b_hoistedExpr_sizes[1];
    eml_loop_ub = eml_b_hoistedExpr_sizes[1] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_d_hoistedExpr_data[eml_b_loop_ub] =
        eml_b_hoistedExpr_data[eml_b_loop_ub] / 20.0;
    }
    m_c_power(eml_d_hoistedExpr_data, eml_qb_hoistedExpr_sizes, eml_b_yR_data,
     eml_tmp_sizes);
    eml_loop_ub = eml_filterCount - 1;
    eml_b_loop_ub = eml_tmp_sizes[1] - 1;
    for(eml_y = 0; eml_y <= eml_b_loop_ub; eml_y++) {
      eml_MOCcontrol[eml_loop_ub + 11 * ((int32_T)eml_yL_data[eml_y] - 1)] =
        eml_b_yR_data[eml_tmp_sizes[0] * eml_y];
    }
    /* Replace the shortened version */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /*  SOme calculations for the GUI monitoring */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    eml_MOCend[eml_filterCount - 1] = eml_MOCcontrol[(eml_filterCount - 1) + 11
      * ((int32_T)eml_numSamples - 1)];
    m_abs(eml_yR_data, eml_yR_sizes, eml_b_hoistedExpr_data,
     eml_b_hoistedExpr_sizes);
    eml_peakOPL[eml_filterCount - 1] = 20.0 *
    log10(m_max(eml_b_hoistedExpr_data, eml_b_hoistedExpr_sizes) / 0.00002);
    m_abs(eml_c_yR_data, eml_b_yR_sizes, eml_b_hoistedExpr_data,
     eml_b_hoistedExpr_sizes);
    eml_peakOPR[eml_filterCount - 1] = 20.0 *
    log10(m_max(eml_b_hoistedExpr_data, eml_b_hoistedExpr_sizes) / 0.00002);
    /*      rmsOPL(filterCount) = 20*log10( sqrt(mean(yL.^2)) /2e-5  );  */
    /*      rmsOPR(filterCount) = 20*log10( sqrt(mean(yR.^2)) /2e-5 );  */
    eml_rmsOPL[eml_filterCount - 1] = 20.0 *
      log10(fabs(eml_yR_data[eml_yR_sizes[1] - 1]) / 0.00002);
    eml_rmsOPR[eml_filterCount - 1] = 20.0 *
      log10(fabs(eml_c_yR_data[eml_b_yR_sizes[1] - 1]) / 0.00002);
    eml_rb_hoistedExpr_sizes[0] = (int32_T)eml_numSamples;
    eml_loop_ub = (int32_T)eml_numSamples - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_mb_hoistedExpr_data[eml_b_loop_ub] = 1.0 + (real_T)eml_b_loop_ub;
    }
    eml_b_hoistedExpr_sizes[0] = 1;
    eml_tmp_sizes[0] = 1;
    eml_tmp_sizes[1] = eml_rb_hoistedExpr_sizes[0];
    eml_b_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
    eml_loop_ub = eml_rb_hoistedExpr_sizes[0] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_b_hoistedExpr_data[eml_b_loop_ub] =
        eml_mb_hoistedExpr_data[eml_b_loop_ub];
    }
    eml_b_MOCthreshold_dBOP = eml_mainGain[eml_filterCount - 1];
    eml_sb_hoistedExpr_sizes[0] = (int32_T)eml_numSamples;
    eml_loop_ub = (int32_T)eml_numSamples - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_nb_hoistedExpr_data[eml_b_loop_ub] = 1.0 + (real_T)eml_b_loop_ub;
    }
    eml_loop_ub = eml_sb_hoistedExpr_sizes[0] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_yL_data[eml_b_loop_ub] = eml_nb_hoistedExpr_data[eml_b_loop_ub];
    }
    eml_tmp_sizes[0] = 1;
    eml_tmp_sizes[1] = eml_b_hoistedExpr_sizes[1];
    eml_c_frameBufferL_sizes[0] = 1;
    eml_c_frameBufferL_sizes[1] = eml_tmp_sizes[1];
    eml_loop_ub = eml_tmp_sizes[1] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      for(eml_y = 0; eml_y <= 0; eml_y = 1) {
        eml_b_yR_data[eml_b_loop_ub] =
          eml_frameBufferL[(int32_T)eml_b_hoistedExpr_data[eml_b_loop_ub] - 1] +
          eml_yR_data[
          eml_yR_sizes[0] * eml_b_loop_ub] * eml_b_MOCthreshold_dBOP;
      }
    }
    eml_loop_ub = eml_c_frameBufferL_sizes[1] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_frameBufferL[(int32_T)eml_yL_data[eml_b_loop_ub] - 1] =
        eml_b_yR_data[eml_b_loop_ub];
    }
    eml_tb_hoistedExpr_sizes[0] = (int32_T)eml_numSamples;
    eml_loop_ub = (int32_T)eml_numSamples - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_ob_hoistedExpr_data[eml_b_loop_ub] = 1.0 + (real_T)eml_b_loop_ub;
    }
    eml_b_hoistedExpr_sizes[0] = 1;
    eml_tmp_sizes[0] = 1;
    eml_tmp_sizes[1] = eml_tb_hoistedExpr_sizes[0];
    eml_b_hoistedExpr_sizes[1] = eml_tmp_sizes[1];
    eml_loop_ub = eml_tb_hoistedExpr_sizes[0] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_b_hoistedExpr_data[eml_b_loop_ub] =
        eml_ob_hoistedExpr_data[eml_b_loop_ub];
    }
    eml_b_MOCthreshold_dBOP = eml_mainGain[eml_filterCount - 1];
    eml_ub_hoistedExpr_sizes[0] = (int32_T)eml_numSamples;
    eml_loop_ub = (int32_T)eml_numSamples - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_pb_hoistedExpr_data[eml_b_loop_ub] = 1.0 + (real_T)eml_b_loop_ub;
    }
    eml_loop_ub = eml_ub_hoistedExpr_sizes[0] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_yL_data[eml_b_loop_ub] = eml_pb_hoistedExpr_data[eml_b_loop_ub];
    }
    eml_tmp_sizes[0] = 1;
    eml_tmp_sizes[1] = eml_b_hoistedExpr_sizes[1];
    eml_frameBufferR_sizes[0] = 1;
    eml_frameBufferR_sizes[1] = eml_tmp_sizes[1];
    eml_loop_ub = eml_tmp_sizes[1] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      for(eml_y = 0; eml_y <= 0; eml_y = 1) {
        eml_b_yR_data[eml_b_loop_ub] =
          eml_frameBufferR[(int32_T)eml_b_hoistedExpr_data[eml_b_loop_ub] - 1] +
          eml_c_yR_data[
          eml_b_yR_sizes[0] * eml_b_loop_ub] * eml_b_MOCthreshold_dBOP;
      }
    }
    eml_loop_ub = eml_frameBufferR_sizes[1] - 1;
    for(eml_b_loop_ub = 0; eml_b_loop_ub <= eml_loop_ub; eml_b_loop_ub++) {
      eml_frameBufferR[(int32_T)eml_yL_data[eml_b_loop_ub] - 1] =
        eml_b_yR_data[eml_b_loop_ub];
    }
  }
  /*  BF channel */
}
void EssexAidProcessVFrameSwitchable_initialize(void)
{
  rt_InitInfAndNaN(8U);
}
void EssexAidProcessVFrameSwitchable_terminate(void)
{
  /* (no terminate code required) */
}
/* End of Embedded MATLAB Coder code generation (EssexAidProcessVFrameSwitchable.c) */
