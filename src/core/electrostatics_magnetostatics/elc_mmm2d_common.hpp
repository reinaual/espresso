
#ifndef ELC_MMM2D_COMMON
#define ELC_MMM2D_COMMON

#ifdef ELECTROSTATICS

struct SCCache {
  double s, c;
};


inline double *block(double *p, int index, int size);
inline double *blwentry(double *p, int index, int e_size);
inline double *abventry(double *p, int index, int e_size);

template <size_t dir> void elc_mmm2d_common_add_force_add_force();

void elc_mmm2d_common_add_P_force();
void elc_mmm2d_common_add_Q_force();

#endif /* ELECTROSTATICS */

#endif /* ELC_MMM2D_COMMON */
