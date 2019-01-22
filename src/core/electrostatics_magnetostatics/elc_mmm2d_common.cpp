
#include "elc_mmm2d_common.hpp"

#include "grid.hpp"
#include "config.hpp"
#include "layered.hpp"
#include <cmath>

#ifdef ELECTROSTATICS

// #define CHECKPOINTS
#if 0
#define LOG_FORCES(x) x
#else
#define LOG_FORCES(x)
#endif

#define POQESP 0
#define POQECP 1
#define POQESM 2
#define POQECM 3

#define PQESSP 0
#define PQESCP 1
#define PQECSP 2
#define PQECCP 3
#define PQESSM 4
#define PQESCM 5
#define PQECSM 6
#define PQECCM 7

// extern std::vector<SCCache> scxcache;
extern double *partblk;
extern double *gblcblk;


/* block indexing - has to fit to the PQ block definitions above.
   size gives the full size of the data block,
   e_size is the size of only the top or bottom half, i.e. half of size.
*/

inline double *block(double *p, int index, int size) {
  return &p[index * size];
}

inline double *blwentry(double *p, int index, int e_size) {
  return &p[2 * index * e_size];
}

inline double *abventry(double *p, int index, int e_size) {
  return &p[(2 * index + 1) * e_size];
}



template <size_t dir> void elc_mmm2d_common_add_force() {
  constexpr const auto size = 4;

  auto ic = 0;
  for (int c = 1; c <= n_layers; c++) {
    auto const np = cells[c].n;
    auto const part = cells[c].part;
    auto const othcblk = block(gblcblk, c - 1, size);

    for (int i = 0; i < np; i++) {
      part[i].f.f[dir] += part[i].p.q * (
                          partblk[size * ic + POQESM] * othcblk[POQECP] -
                          partblk[size * ic + POQECM] * othcblk[POQESP] +
                          partblk[size * ic + POQESP] * othcblk[POQECM] -
                          partblk[size * ic + POQECP] * othcblk[POQESM]);
      part[i].f.f[2] += part[i].p.q * (
                        partblk[size * ic + POQECM] * othcblk[POQECP] +
                        partblk[size * ic + POQESM] * othcblk[POQESP] -
                        partblk[size * ic + POQECP] * othcblk[POQECM] -
                        partblk[size * ic + POQESP] * othcblk[POQESM]);

      LOG_FORCES(fprintf(stderr, "%d: part %d force %10.3g %10.3g %10.3g\n",
                         this_node, part[i].p.identity, part[i].f.f[0],
                         part[i].f.f[1], part[i].f.f[2]));
      ic++;
    }
  }
}


void elc_mmm2d_common_add_P_force() { elc_mmm2d_common_add_force<0>(); }
void elc_mmm2d_common_add_Q_force() { elc_mmm2d_common_add_force<1>(); }


#endif /* ELECTROSTATICS */
