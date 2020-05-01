/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/** \file
 *  Implementation of \ref elc.hpp.
 */
#include "Particle.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "mmm-common.hpp"
#include "pressure.hpp"
#include <cmath>
#include <mpi.h>

#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/elc.hpp"
#include "electrostatics_magnetostatics/p3m.hpp"

#ifdef P3M

/****************************************
 * LOCAL DEFINES
 ****************************************/

/** Largest reasonable cutoff for far formula */
#define MAXIMAL_FAR_CUT 50

/****************************************
 * LOCAL VARIABLES
 ****************************************/

/** \name Inverse box dimensions and derived constants */
/*@{*/
static double ux, uy, uz, height_inverse;
/*@}*/

ELC_struct elc_params = {1e100, 10,    1, 0, true, true, false, 1,
                         1,     false, 0, 0, 0,    0,    0.0};

/****************************************
 * LOCAL ARRAYS
 ****************************************/

/** \name Product decomposition data organization
 *  For the cell blocks it is assumed that the lower blocks part is in the
 *  lower half. This has to have positive sign, so that has to be first.
 */
/*@{*/
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
/*@}*/

/** temporary buffers for product decomposition (chi-sum) */
static std::vector<double> part_fac;
/** image summation array (X-sum) with positive indices from subset L_- and
 * negative indices from subset L_+*/
// static Utils::VectorXd<8> image_sum{};

/** structure for caching sin and cos values */
typedef struct {
  double s, c;
} SCCache;

/** Cached sin/cos values along the x-axis and y-axis */
/*@{*/
static std::vector<SCCache> scxcache;
static std::vector<SCCache> scycache;
/*@}*/

/****************************************
 * LOCAL FUNCTIONS
 ****************************************/

static void distribute(int size);
/** \name p,q <> 0 per frequency code */
/*@{*/
static Utils::VectorXd<8> setup_PQ(int p, int q, double omega,
                                   const ParticleRange &particles);
static void add_PQ_force(const ParticleRange &particles,
                         Utils::VectorXd<8> &part_sum);
static double PQ_energy(double fpq,
                        int n_particles,
                        Utils::VectorXd<8> &part_sum);
/*@}*/
static void add_dipole_force(const ParticleRange &particles);
static double dipole_energy(const ParticleRange &particles);
static double z_energy(const ParticleRange &particles);
static void add_z_force(const ParticleRange &particles);

void ELC_setup_constants() {
  ux = 1 / box_geo.length()[0];
  uy = 1 / box_geo.length()[1];
  uz = 1 / box_geo.length()[2];

  height_inverse = 1 / elc_params.h;
}

/**
 * @brief Calculated cached sin/cos values for one direction.
 *
 * @tparam Index of the dimension to consider (e.g. 0 for x ...).
 *
 * @param particles Particles to calculate values for
 * @param n_freq Number of frequencies to calculate per particle
 * @param u Inverse box length
 * @return Calculated values.
 */
template <size_t dir>
static std::vector<SCCache> sc_cache(const ParticleRange &particles, int n_freq,
                                     double u) {
  auto const n_part = particles.size();
  std::vector<SCCache> ret((n_freq + 1) * n_part);

  for (size_t freq = 0; freq <= n_freq; freq++) {
    double pref = C_2PI * u * freq;

    size_t o = freq * n_part;
    for (auto const &part : particles) {
      auto const arg = pref * part.r.p[dir];
      ret[o++] = {sin(arg), cos(arg)};
    }
  }

  return ret;
}

static void prepare_sc_cache(const ParticleRange &particles, int n_freq_x,
                             double u_x, int n_freq_y, double u_y) {
  scxcache = sc_cache<0>(particles, n_freq_x, u_x);
  scycache = sc_cache<1>(particles, n_freq_y, u_y);
}

/*****************************************************************/
/* data distribution */
/*****************************************************************/

inline void clear_vec(double *pdc, int size) {
  for (int i = 0; i < size; i++)
    pdc[i] = 0;
}

inline void copy_vec(double *pdc_d, double const *pdc_s, int size) {
  for (int i = 0; i < size; i++)
    pdc_d[i] = pdc_s[i];
}

inline void add_vec(double *pdc_d, double const *pdc_s1, double const *pdc_s2,
                    int size) {
  for (int i = 0; i < size; i++)
    pdc_d[i] = pdc_s1[i] + pdc_s2[i];
}

inline void addscale_vec(double *pdc_d, double scale, double const *pdc_s1,
                         double const *pdc_s2, int size) {
  for (int i = 0; i < size; i++)
    pdc_d[i] = scale * pdc_s1[i] + pdc_s2[i];
}

inline void scale_vec(double scale, double *pdc, int size) {
  for (int i = 0; i < size; i++)
    pdc[i] *= scale;
}

inline double *block(double *p, int index, int size) {
  return &p[index * size];
}

void distribute(int size) {
  double send_buf[8];
// copy_vec(send_buf, image_sum, size);
//  MPI_Allreduce(send_buf, image_sum, size, MPI_DOUBLE, MPI_SUM, comm_cart);
}

/*****************************************************************/
/* dipole terms */
/*****************************************************************/

/** Calculate the dipole force.
 *  See @cite yeh99a.
 */

/**
 * @brief Calculate and add the dipole correction force with forced tinfoil
 * boundary condition in the P3M method
 * @param particles Particle to calculate dipole moment and correction for
 */
static void add_dipole_force(const ParticleRange &particles) {
  const double pref = coulomb.prefactor * 4 * Utils::pi() * ux * uy * uz;

  const auto local_particles = particles;
  const double shift = 0.5 * elc_params.h;

  Utils::Vector3d moments{};

  /*
  moments[0] = 0; // sum q_i
  moments[1] = 0; // sum q_i z_i
  moments[2] = 0; // sum q_i (z_i - L/2)
  */

  for (auto const &p : local_particles) {
    const double zpos = p.r.p[2];
    const double q = p.p.q;

    moments[0] += q;
    moments[1] += q * zpos;
    moments[2] += q * (zpos - shift);

    if (elc_params.dielectric_contrast_on) {
      if (zpos < elc_params.space_layer) {
        moments[0] += elc_params.delta_mid_bot * q;
        moments[2] += elc_params.delta_mid_bot * q * (-zpos - shift);
      }
      if (zpos > (elc_params.h - elc_params.space_layer)) {
        moments[0] += elc_params.delta_mid_top * q;
        moments[2] += elc_params.delta_mid_top * q *
                      (2 * elc_params.h - zpos - shift);
      }
    }
  }

  moments[0] *= pref;
  moments[1] *= pref * height_inverse / uz;
  moments[2] *= pref;

  moments = boost::mpi::all_reduce(comm_cart, moments, std::plus<>());

  // Yeh + Berkowitz dipole term @cite yeh99a
  double field_tot = moments[2];

  // Const. potential contribution
  if (elc_params.const_pot) {
    const double field_induced = moments[1];
    const double field_applied = elc_params.pot_diff * height_inverse;
    field_tot -= field_applied + field_induced;
  }

  for (auto &p : local_particles) {
    p.f.f[2] -= field_tot * p.p.q;
  }
}

/**
 * @brief Calculate the dipole correction energy for vacuum boundary conditions
 * @param particles Particle to calculate correction for
 * @return Dipole correction energy
 */
static double dipole_energy(const ParticleRange &particles) {
  const double pref = 2 * Utils::pi() * ux * uy * uz;
  double shift = 0.5 * elc_params.h;

  Utils::VectorXd<7> moments{};
  /*
  moments[0] = 0; // sum q_i               primary box
  moments[1] = 0; // sum q_i               boundary layers
  moments[2] = 0; // sum q_i (z_i - L/2)   primary box
  moments[3] = 0; // sum q_i (z_i - L/2)   boundary layers
  moments[4] = 0; // sum q_i (z_i - L/2)^2 primary box
  moments[5] = 0; // sum q_i (z_i - L/2)^2 boundary layers
  moments[6] = 0; // sum q_i z_i           primary box
   */

  for (auto &p : particles) {
    const double zpos = p.r.p[2];
    const double q = p.p.q;

    moments[0] += q;
    moments[2] += q * (zpos - shift);
    moments[4] += q * (Utils::sqr(zpos - shift));
    moments[6] += q * zpos;

    if (elc_params.dielectric_contrast_on) {
      if (zpos < elc_params.space_layer) {
        moments[1] += elc_params.delta_mid_bot * q;
        moments[3] += elc_params.delta_mid_bot * q * (-zpos - shift);
        moments[5] +=
            elc_params.delta_mid_bot * q * (Utils::sqr(-zpos - shift));
      }
      if (zpos > (elc_params.h - elc_params.space_layer)) {
        moments[1] += elc_params.delta_mid_top * q;
        moments[3] +=
            elc_params.delta_mid_top * q * (2 * elc_params.h - zpos - shift);
        moments[5] += elc_params.delta_mid_top * q *
                      (Utils::sqr(2 * elc_params.h - zpos - shift));
      }
    }
  }

  moments = boost::mpi::all_reduce(comm_cart, moments, std::plus<>());

  double eng = Utils::sqr(moments[2]) - moments[0] * moments[4] -
               Utils::sqr(elc_params.h * moments[0]) / 12;

  if (elc_params.dielectric_contrast_on) {
    if (elc_params.const_pot) {
      // zero potential difference contribution
      eng += pref * height_inverse / uz * Utils::sqr(moments[6]);
      // external potential shift contribution
      eng -= 2 * elc_params.pot_diff * height_inverse * moments[6];
    }

    /* counter the P3M homogeneous background contribution to the
       boundaries. We never need that, since a homogeneous background
       spanning the artificial boundary layers is aphysical. */
    eng += pref * (-(moments[1] * moments[4] + moments[0] * moments[5]) -
                   (1. - 2. / 3.) * moments[0] * moments[1] *
                       Utils::sqr(box_geo.length()[2]));
  }

  return this_node == 0 ? pref * eng : 0;
}

/*****************************************************************/

inline double image_sum_b(double q, double z) {
  double shift = 0.5 * box_geo.length()[2];
  double fac = elc_params.delta_mid_top * elc_params.delta_mid_bot;
  double image_sum =
      (q / (1.0 - fac) * (z - 2.0 * fac * box_geo.length()[2] / (1.0 - fac))) -
      q * shift / (1 - fac);
  return image_sum;
}

inline double image_sum_t(double q, double z) {
  double shift = 0.5 * box_geo.length()[2];
  double fac = elc_params.delta_mid_top * elc_params.delta_mid_bot;
  double image_sum =
      (q / (1.0 - fac) * (z + 2.0 * fac * box_geo.length()[2] / (1.0 - fac))) -
      q * shift / (1 - fac);
  return image_sum;
}

/*****************************************************************/
static double z_energy(const ParticleRange &particles) {
  double pref = coulomb.prefactor * 2 * M_PI * ux * uy;
  int size = 4;

  Utils::Vector4d image_sum{};

  double eng = 0;
  /* for nonneutral systems, this shift gives the background contribution
     (rsp. for this shift, the DM of the background is zero) */
  double shift = 0.5 * box_geo.length()[2];

  if (elc_params.dielectric_contrast_on) {
    if (elc_params.const_pot) {
      for (auto &p : particles) {
        image_sum[0] += p.p.q;
        image_sum[1] += p.p.q * (p.r.p[2] - shift);
        if (p.r.p[2] < elc_params.space_layer) {
          image_sum[2] -= elc_params.delta_mid_bot * p.p.q;
          image_sum[3] -= elc_params.delta_mid_bot * p.p.q * (-p.r.p[2] - shift);
        }
        if (p.r.p[2] > (elc_params.h - elc_params.space_layer)) {
          image_sum[2] += elc_params.delta_mid_top * p.p.q;
          image_sum[3] += elc_params.delta_mid_top * p.p.q *
                        (2 * elc_params.h - p.r.p[2] - shift);
        }
      }
    } else {
      double delta = elc_params.delta_mid_top * elc_params.delta_mid_bot;
      double fac_delta_mid_bot = elc_params.delta_mid_bot / (1 - delta);
      double fac_delta_mid_top = elc_params.delta_mid_top / (1 - delta);
      double fac_delta = delta / (1 - delta);

      for (auto &p : particles) {
        image_sum[0] += p.p.q;
        image_sum[1] += p.p.q * (p.r.p[2] - shift);
        if (elc_params.dielectric_contrast_on) {
          if (p.r.p[2] < elc_params.space_layer) {
            image_sum[2] += fac_delta * (elc_params.delta_mid_bot + 1) * p.p.q;
            image_sum[3] +=
                p.p.q * (image_sum_b(elc_params.delta_mid_bot * delta,
                                     -(2 * elc_params.h + p.r.p[2])) +
                         image_sum_b(delta, -(2 * elc_params.h - p.r.p[2])));
          } else {
            image_sum[2] +=
                fac_delta_mid_bot * (1 + elc_params.delta_mid_top) * p.p.q;
            image_sum[3] +=
                p.p.q * (image_sum_b(elc_params.delta_mid_bot, -p.r.p[2]) +
                         image_sum_b(delta, -(2 * elc_params.h - p.r.p[2])));
          }
          if (p.r.p[2] > (elc_params.h - elc_params.space_layer)) {
            // note the minus sign here which is required due to |z_i-z_j|
            image_sum[2] -= fac_delta * (elc_params.delta_mid_top + 1) * p.p.q;
            image_sum[3] -=
                p.p.q * (image_sum_t(elc_params.delta_mid_top * delta,
                                     4 * elc_params.h - p.r.p[2]) +
                         image_sum_t(delta, 2 * elc_params.h + p.r.p[2]));
          } else {
            // note the minus sign here which is required due to |z_i-z_j|
            image_sum[2] -=
                fac_delta_mid_top * (1 + elc_params.delta_mid_bot) * p.p.q;
            image_sum[3] -=
                p.p.q * (image_sum_t(elc_params.delta_mid_top,
                                     2 * elc_params.h - p.r.p[2]) +
                         image_sum_t(delta, 2 * elc_params.h + p.r.p[2]));
          }
        }
      }
    }
  }
  distribute(size);

  if (this_node == 0)
    eng -= pref * (image_sum[1] * image_sum[2] - image_sum[0] * image_sum[3]);

  return eng;
}

/*****************************************************************/
static void add_z_force(const ParticleRange &particles) {
  double pref = coulomb.prefactor * 2 * M_PI * ux * uy;

  double zforce = 0;

  if (elc_params.dielectric_contrast_on) {
    auto local_particles = particles;
    int size = 1;
    if (elc_params.const_pot) {
      /* just counter the 2 pi |z| contribution stemming from P3M */
      for (auto &p : local_particles) {
        if (p.r.p[2] < elc_params.space_layer)
          zforce -= elc_params.delta_mid_bot * p.p.q;
        if (p.r.p[2] > (elc_params.h - elc_params.space_layer))
          zforce += elc_params.delta_mid_top * p.p.q;
      }
    } else {
      double delta = elc_params.delta_mid_top * elc_params.delta_mid_bot;
      double fac_delta_mid_bot = elc_params.delta_mid_bot / (1 - delta);
      double fac_delta_mid_top = elc_params.delta_mid_top / (1 - delta);
      double fac_delta = delta / (1 - delta);

      for (auto &p : local_particles) {
        if (p.r.p[2] < elc_params.space_layer) {
          zforce += fac_delta * (elc_params.delta_mid_bot + 1) * p.p.q;
        } else {
          zforce +=
              fac_delta_mid_bot * (1 + elc_params.delta_mid_top) * p.p.q;
        }

        if (p.r.p[2] > (elc_params.h - elc_params.space_layer)) {
          // note the minus sign here which is required due to |z_i-z_j|
          zforce -= fac_delta * (elc_params.delta_mid_top + 1) * p.p.q;
        } else {
          // note the minus sign here which is required due to |z_i-z_j|
          zforce -=
              fac_delta_mid_top * (1 + elc_params.delta_mid_bot) * p.p.q;
        }
      }
    }

    zforce *= pref;

    zforce = boost::mpi::all_reduce(comm_cart, zforce, std::plus<>());

    for (auto &p : local_particles) {
      p.f.f[2] += zforce * p.p.q;
    }
  }
}

/*****************************************************************/
/* PQ particle blocks */
/*****************************************************************/

/**
 * @brief Calculate the necessary fourier factors for the summation
 * @param index_x x-direction index for caches cos/sin values
 * @param index_y y-direction index for caches cos/sin values
 * @return Calculated fourier factors
 */
inline Utils::VectorXd<8> fourier_factors(const size_t index_x,
                                          const size_t index_y) {
  Utils::VectorXd<8> temp_sum;
  const auto x_cache = scxcache[index_x], y_cache = scycache[index_y];
  temp_sum[PQESSP] = x_cache.s * y_cache.s;
  temp_sum[PQESCP] = x_cache.s * y_cache.c;
  temp_sum[PQECSP] = x_cache.c * y_cache.s;
  temp_sum[PQECCP] = x_cache.c * y_cache.c;
  temp_sum[PQESSM] = temp_sum[PQESSP];
  temp_sum[PQESCM] = temp_sum[PQESCP];
  temp_sum[PQECSM] = temp_sum[PQECSP];
  temp_sum[PQECCM] = temp_sum[PQECCP];

  return temp_sum;
}

/**
 * @brief Apply the factor_p to the first 4 elements of the array and factor_n to the second half
 * @param vector vector to multiply the factors with
 * @param factor_p factor for the first half
 * @param factor_n factor for the second half
 * @return new vector with applied factors
 */
inline Utils::VectorXd<8> apply_factors(Utils::VectorXd<8> vector,
                                        const double factor_p,
                                        const double factor_n) {
  /* this function is returning a copy of the vector with the first half scaled
   * by the first factor and the second half with the second factor */
  for (unsigned int i = 0; i < 4; ++i) {
    vector[i] *= factor_p;
    vector[i + 4] *= factor_n;
  }
  return vector;
}

/**
 * @brief Calculate the Lprime summation prefactors for the image_sum (X-sum)
 * @param z z position
 * @param omega 2 * pi * f_pq
 * @return Calculated prefactor
 */
inline double sum_prefactor_prime(const double z, const double omega) {
  /* the divisor has the opposite sign compared to the paper but matches the
     previous implementation */
  return -exp(-omega * z) / expm1(2 * omega * box_geo.length()[2]);
}

/**
 * @brief Calulate pq summation parts
 *
 * @param p first fourier index
 * @param q second fourier index
 * @param omega 2 * pi * f_pq
 * @param particles Particle to calculate sum for
 *
 * @return Calculated sums.
 */
static Utils::VectorXd<8> setup_PQ(const int p, const int q, const double omega,
                                   const ParticleRange &particles) {
  Utils::VectorXd<8> part_sum{};

  double const pref = -coulomb.prefactor * 8 * M_PI * ux * uy /
                      (expm1(omega * box_geo.length()[2]));
  double const pref_di = coulomb.prefactor * 8 * M_PI * ux * uy;
  int const size = 8;
  double lclimgebot[8], lclimgetop[8], lclimge[8];
  double fac_delta_mid_bot = 1, fac_delta_mid_top = 1, fac_delta = 1;
  if (elc_params.dielectric_contrast_on) {
    double fac_elc =
        1.0 / (1 - elc_params.delta_mid_top * elc_params.delta_mid_bot *
                       exp(-omega * 2 * elc_params.h));
    fac_delta_mid_bot = elc_params.delta_mid_bot * fac_elc;
    fac_delta_mid_top = elc_params.delta_mid_top * fac_elc;
    fac_delta = fac_delta_mid_bot * elc_params.delta_mid_top;
  }

  clear_vec(lclimge, size);

  size_t ic = 0;
  auto const ox = static_cast<size_t>(p * particles.size());
  auto const oy = static_cast<size_t>(q * particles.size());
  for (auto const &p : particles) {
    const size_t index_x = ox + ic;
    const size_t index_y = oy + ic;

    const double zpos = p.r.p[2];

    // setup vector with fourier factors
    Utils::VectorXd<8> temp_sum = p.p.q * fourier_factors(index_x, index_y);

    // add prefactor scaled part to particle sums (chi-sums in eq (3.3))
    const double chi_factor = exp(omega * zpos);
    const double inv_chi_factor = 1. / chi_factor;
    part_sum += apply_factors(temp_sum, chi_factor, inv_chi_factor);

    // add prefactor scaled part to the image sums (X-sums in eq (3.6))
    const double x_factor = sum_prefactor_prime(zpos, omega);
    const double inv_x_factor =
        sum_prefactor_prime(box_geo.length()[2] - zpos, omega);
    // TODO this might be written more compact with a proxy class...
    const auto out = apply_factors(temp_sum, x_factor, inv_x_factor);
    for (int i = 0; i < 8; i++) {
      part_fac[ic + i] = out[i];
    }

    /*
    if (elc_params.dielectric_contrast_on) {
      if (p.r.p[2] < elc_params.space_layer) { // handle the lower case first
        // change e to take into account the z position of the images

        e = exp(-omega * p.r.p[2]);
        scale = p.p.q * elc_params.delta_mid_bot;

        lclimgebot[PQESSM] = scxcache[index_x].s * scycache[index_y].s / e;
        lclimgebot[PQESCM] = scxcache[index_x].s * scycache[index_y].c / e;
        lclimgebot[PQECSM] = scxcache[index_x].c * scycache[index_y].s / e;
        lclimgebot[PQECCM] = scxcache[index_x].c * scycache[index_y].c / e;

        lclimgebot[PQESSP] = scxcache[index_x].s * scycache[index_y].s * e;
        lclimgebot[PQESCP] = scxcache[index_x].s * scycache[index_y].c * e;
        lclimgebot[PQECSP] = scxcache[index_x].c * scycache[index_y].s * e;
        lclimgebot[PQECCP] = scxcache[index_x].c * scycache[index_y].c * e;

        addscale_vec(image_sum, scale, lclimgebot, image_sum, size);

        e = (exp(omega * (-p.r.p[2] - 2 * elc_params.h)) *
                 elc_params.delta_mid_bot +
             exp(omega * (p.r.p[2] - 2 * elc_params.h))) *
            fac_delta * p.p.q;

      } else {

        e = (exp(omega * (-p.r.p[2])) +
             exp(omega * (p.r.p[2] - 2 * elc_params.h)) *
                 elc_params.delta_mid_top) *
            fac_delta_mid_bot * p.p.q;
      }

      lclimge[PQESSP] += scxcache[index_x].s * scycache[index_y].s * e;
      lclimge[PQESCP] += scxcache[index_x].s * scycache[index_y].c * e;
      lclimge[PQECSP] += scxcache[index_x].c * scycache[index_y].s * e;
      lclimge[PQECCP] += scxcache[index_x].c * scycache[index_y].c * e;

      if (p.r.p[2] > (elc_params.h -
                      elc_params.space_layer)) { // handle the upper case now

        e = exp(omega * (2 * elc_params.h - p.r.p[2]));
        scale = p.p.q * elc_params.delta_mid_top;

        lclimgetop[PQESSM] = scxcache[index_x].s * scycache[index_y].s / e;
        lclimgetop[PQESCM] = scxcache[index_x].s * scycache[index_y].c / e;
        lclimgetop[PQECSM] = scxcache[index_x].c * scycache[index_y].s / e;
        lclimgetop[PQECCM] = scxcache[index_x].c * scycache[index_y].c / e;

        lclimgetop[PQESSP] = scxcache[index_x].s * scycache[index_y].s * e;
        lclimgetop[PQESCP] = scxcache[index_x].s * scycache[index_y].c * e;
        lclimgetop[PQECSP] = scxcache[index_x].c * scycache[index_y].s * e;
        lclimgetop[PQECCP] = scxcache[index_x].c * scycache[index_y].c * e;

        addscale_vec(image_sum, scale, lclimgetop, image_sum, size);

        e = (exp(omega * (p.r.p[2] - 4 * elc_params.h)) *
                 elc_params.delta_mid_top +
             exp(omega * (-p.r.p[2] - 2 * elc_params.h))) *
            fac_delta * p.p.q;

      } else {

        e = (exp(omega * (p.r.p[2] - 2 * elc_params.h)) +
             exp(omega * (-p.r.p[2] - 2 * elc_params.h)) *
                 elc_params.delta_mid_bot) *
            fac_delta_mid_top * p.p.q;
      }

      lclimge[PQESSM] += scxcache[index_x].s * scycache[index_y].s * e;
      lclimge[PQESCM] += scxcache[index_x].s * scycache[index_y].c * e;
      lclimge[PQECSM] += scxcache[index_x].c * scycache[index_y].s * e;
      lclimge[PQECCM] += scxcache[index_x].c * scycache[index_y].c * e;
    }
    */

    ic++;
  }

  /*
  if (elc_params.dielectric_contrast_on) {
    scale_vec(pref_di, lclimge, size);
    add_vec(image_sum, image_sum, lclimge, size);
  } */

  return boost::mpi::all_reduce(comm_cart, part_sum, std::plus<>());
}

static void add_PQ_force(const ParticleRange &particles,
                         const Utils::VectorXd<8> &part_sum) {
  const double pref_z = 2 * Utils::pi() * ux * uy * coulomb.prefactor;

  size_t ic = 0;
  for (auto &p : particles) {
    const size_t offset = 8 * ic;

    p.f.f[2] += pref_z * (part_fac[offset + PQECCM] * part_sum[PQECCP] +
                          part_fac[offset + PQESCM] * part_sum[PQESCP] +
                          part_fac[offset + PQECSM] * part_sum[PQECSP] +
                          part_fac[offset + PQESSM] * part_sum[PQESSP] -
                          part_fac[offset + PQECCP] * part_sum[PQECCM] -
                          part_fac[offset + PQESCP] * part_sum[PQESCM] -
                          part_fac[offset + PQECSP] * part_sum[PQECSM] -
                          part_fac[offset + PQESSP] * part_sum[PQESSM]);
    ic++;
  }
}

/**
 * @brief Calculate the factors for the pq-summation
 * @param fpq factor f_pq
 * @param n_particles number of particles to iterate over
 * @param part_sum summation values
 * @return Calculated energy contribution.
 */
static double PQ_energy(const double fpq, const size_t n_particles,
                        const Utils::VectorXd<8> &part_sum) {
  // return p,q-summation part of the energy according to equation (3.10)
  double eng = 0;

  // calculate the image sums
  Utils::VectorXd<8> image_sum{};
  for (size_t i = 0; i < n_particles; ++i) {
    const size_t offset = i * 8;
    // TODO this can be condensed with a proxy class...
    for (int j = 0; j < 8; j++) {
      image_sum[j] = part_fac[offset + j];
    }
  }

  // first sum in (3.10)
  eng += part_sum[PQECCM] * image_sum[PQECCP] +
         part_sum[PQESCM] * image_sum[PQESCP] +
         part_sum[PQECSM] * image_sum[PQECSP] +
         part_sum[PQESSM] * image_sum[PQESSP];

  // second sum in (3.10)
  eng += image_sum[PQECCM] * part_sum[PQECCP] +
         image_sum[PQESCM] * part_sum[PQESCP] +
         image_sum[PQECSM] * part_sum[PQECSP] +
         image_sum[PQESSM] * part_sum[PQESSP];

  return eng / fpq;
}

/*****************************************************************/
/* main loops */
/*****************************************************************/

void ELC_add_force(const ParticleRange &particles) {
  double omega;

  auto const n_scxcache = size_t(ceil(elc_params.far_cut * box_geo.length()[0]) + 1);
  auto const n_scycache = size_t(ceil(elc_params.far_cut * box_geo.length()[1]) + 1);

  prepare_sc_cache(particles, n_scxcache, ux, n_scycache, uy);
  part_fac.resize(8 * particles.size());

  add_dipole_force(particles);
  if (elc_params.dielectric_contrast_on) {
    add_z_force(particles);
  }

  const int p_range =
      static_cast<size_t>(elc_params.far_cut * box_geo.length()[0] + 1);
  for (size_t p = 0; p < p_range; ++p) {
    for (size_t q = 0; sqrt(Utils::sqr(ux * p) + Utils::sqr(uy * q)) < elc_params.far_cut; ++q) {
      // skip the p^2 + p^2 == 0 term
      if (p == 0 and q == 0) {
        continue;
      }

      const double fpq = sqrt(Utils::sqr(ux * p) + Utils::sqr(uy * q));
      const auto part_sum = setup_PQ(p, q, 2 * Utils::pi() * fpq, particles);

      add_PQ_force(particles, part_sum);
    }
  }
}

double ELC_energy(const ParticleRange &particles) {
  double eng = dipole_energy(particles);
  if (elc_params.dielectric_contrast_on) {
    eng += z_energy(particles);
  }

  auto const n_scxcache = int(ceil(elc_params.far_cut * box_geo.length()[0]) + 1);
  auto const n_scycache = int(ceil(elc_params.far_cut * box_geo.length()[1]) + 1);
  prepare_sc_cache(particles, n_scxcache, ux, n_scycache, uy);
  part_fac.resize(8 * particles.size());
  auto const n_particles = particles.size();

  double pq_energy = 0;
  const int p_range =
      static_cast<size_t>(elc_params.far_cut * box_geo.length()[0] + 1);
  for (size_t p = 0; p < p_range; ++p) {
    for (size_t q = 0; sqrt(Utils::sqr(ux * p) + Utils::sqr(uy * q)) < elc_params.far_cut; ++q) {
      // skip the p^2 + p^2 == 0 term
      if (p == 0 and q == 0) {
        continue;
      }

      const double fpq = sqrt(Utils::sqr(ux * p) + Utils::sqr(uy * q));
      const auto part_sum = setup_PQ(p, q, 2 * Utils::pi() * fpq, particles);

      pq_energy += PQ_energy(fpq, n_particles, part_sum);
    }
  }
  eng -= 0.5 * ux * uy * pq_energy;

  // TODO scale the end results in the coulomb.cpp/coulomb_inline.cpp at one
  // place to remove inconsistencies
  return coulomb.prefactor * eng;
}

int ELC_tune(double error) {
  double err;
  double h = elc_params.h, lz = box_geo.length()[2];
  double min_inv_boxl = std::min(ux, uy);

  if (elc_params.dielectric_contrast_on) {
    // adjust lz according to dielectric layer method
    lz = elc_params.h + elc_params.space_layer;
  }

  if (h < 0)
    return ES_ERROR;

  elc_params.far_cut = min_inv_boxl;

  do {
    const auto prefactor = 2 * Utils::pi() * elc_params.far_cut;

    const auto sum = prefactor + 2 * (ux + uy);
    const auto den = -expm1(-prefactor * lz);
    const auto num1 = exp(prefactor * (h - lz));
    const auto num2 = exp(-prefactor * (h + lz));

    err = 0.5 / den *
          (num1 * (sum + 1 / (lz - h)) / (lz - h) +
           num2 * (sum + 1 / (lz + h)) / (lz + h));

    elc_params.far_cut += min_inv_boxl;
  } while (err > error && elc_params.far_cut < MAXIMAL_FAR_CUT);
  if (elc_params.far_cut >= MAXIMAL_FAR_CUT)
    return ES_ERROR;
  elc_params.far_cut -= min_inv_boxl;
  elc_params.far_cut2 = Utils::sqr(elc_params.far_cut);

  return ES_OK;
}

/****************************************
 * COMMON PARTS
 ****************************************/

int ELC_sanity_checks() {
  if (!box_geo.periodic(0) || !box_geo.periodic(1) || !box_geo.periodic(2)) {
    runtimeErrorMsg() << "ELC requires periodicity 1 1 1";
    return ES_ERROR;
  }
  /* The product of the two dielectric contrasts should be < 1 for ELC to
     work. This is not the case for two parallel boundaries, which can only
     be treated by the constant potential code */
  if (elc_params.dielectric_contrast_on &&
      (fabs(1.0 - elc_params.delta_mid_top * elc_params.delta_mid_bot) <
       ROUND_ERROR_PREC) &&
      !elc_params.const_pot) {
    runtimeErrorMsg() << "ELC with two parallel metallic boundaries requires "
                         "the const_pot option";
    return ES_ERROR;
  }

  // ELC with non-neutral systems and no fully metallic boundaries does not work
  if (elc_params.dielectric_contrast_on && !elc_params.const_pot &&
      p3m.square_sum_q > ROUND_ERROR_PREC) {
    runtimeErrorMsg() << "ELC does not work for non-neutral systems and "
                         "non-metallic dielectric contrast.";
    return ES_ERROR;
  }

  // Disable this line to make ELC work again with non-neutral systems and
  // metallic boundaries
  if (elc_params.dielectric_contrast_on && elc_params.const_pot &&
      p3m.square_sum_q > ROUND_ERROR_PREC) {
    runtimeErrorMsg() << "ELC does not currently support non-neutral "
                         "systems with a dielectric contrast.";
    return ES_ERROR;
  }

  return ES_OK;
}

void ELC_init() {

  ELC_setup_constants();

  if (elc_params.dielectric_contrast_on) {
    // recalculate the space layer size
    // set the space_layer to be 1/3 of the gap size, so that box = layer
    elc_params.space_layer = (1. / 3.) * elc_params.gap_size;
    // but make sure we leave enough space to not have to bother with
    // overlapping
    // realspace P3M
    double maxsl = elc_params.gap_size - p3m.params.r_cut;
    // and make sure the space layer is not bigger than half the actual
    // simulation box,
    // to avoid overlaps
    if (maxsl > .5 * elc_params.h)
      maxsl = .5 * elc_params.h;
    if (elc_params.space_layer > maxsl) {
      if (maxsl <= 0) {
        runtimeErrorMsg() << "P3M real space cutoff too large for ELC w/ "
                             "dielectric contrast";
      } else
        elc_params.space_layer = maxsl;
    }

    // set the space_box
    elc_params.space_box = elc_params.gap_size - 2 * elc_params.space_layer;
    // reset minimal_dist for tuning
    elc_params.minimal_dist =
        std::min(elc_params.space_box, elc_params.space_layer);
  }

  if (elc_params.far_calculated && (elc_params.dielectric_contrast_on)) {
    if (ELC_tune(elc_params.maxPWerror) == ES_ERROR) {
      runtimeErrorMsg() << "ELC auto-retuning failed, gap size too small";
    }
  }
  if (elc_params.dielectric_contrast_on) {
    p3m.params.additional_mesh[0] = 0;
    p3m.params.additional_mesh[1] = 0;
    p3m.params.additional_mesh[2] = elc_params.space_layer;
  } else {
    p3m.params.additional_mesh[0] = 0;
    p3m.params.additional_mesh[1] = 0;
    p3m.params.additional_mesh[2] = 0;
  }
}

int ELC_set_params(double maxPWerror, double gap_size, double far_cut,
                   bool neutralize, double delta_top, double delta_bot,
                   bool const_pot, double pot_diff) {
  elc_params.maxPWerror = maxPWerror;
  elc_params.gap_size = gap_size;
  elc_params.h = box_geo.length()[2] - gap_size;

  if (delta_top != 0.0 || delta_bot != 0.0) {
    elc_params.dielectric_contrast_on = true;

    elc_params.delta_mid_top = delta_top;
    elc_params.delta_mid_bot = delta_bot;

    // neutralize is automatic with dielectric contrast
    elc_params.neutralize = false;
    // initial setup of parameters, may change later when P3M is finally tuned
    // set the space_layer to be 1/3 of the gap size, so that box = layer
    elc_params.space_layer = (1. / 3.) * gap_size;
    // set the space_box
    elc_params.space_box = gap_size - 2 * elc_params.space_layer;
    // reset minimal_dist for tuning
    elc_params.minimal_dist =
        std::min(elc_params.space_box, elc_params.space_layer);

    // Constant potential parameter setup
    if (const_pot) {
      elc_params.const_pot = true;
      elc_params.pot_diff = pot_diff;
    }
  } else {
    // setup without dielectric contrast
    elc_params.dielectric_contrast_on = false;
    elc_params.const_pot = false;
    elc_params.delta_mid_top = 0;
    elc_params.delta_mid_bot = 0;
    elc_params.neutralize = neutralize;
    elc_params.space_layer = 0;
    elc_params.space_box = elc_params.minimal_dist = gap_size;
  }

  ELC_setup_constants();

  Coulomb::elc_sanity_check();

  elc_params.far_cut = far_cut;
  if (far_cut != -1) {
    elc_params.far_cut2 = Utils::sqr(far_cut);
    elc_params.far_calculated = false;
  } else {
    elc_params.far_calculated = true;
    if (ELC_tune(elc_params.maxPWerror) == ES_ERROR) {
      runtimeErrorMsg() << "ELC tuning failed, gap size too small";
    }
  }
  mpi_bcast_coulomb_params();

  return ES_OK;
}

////////////////////////////////////////////////////////////////////////////////////

void ELC_P3M_self_forces(const ParticleRange &particles) {
  Utils::Vector3d pos;
  double q;

  for (auto &p : particles) {
    if (p.r.p[2] < elc_params.space_layer) {
      q = elc_params.delta_mid_bot * p.p.q * p.p.q;

      pos[0] = p.r.p[0];
      pos[1] = p.r.p[1];
      pos[2] = -p.r.p[2];
      auto const d = get_mi_vector(p.r.p, pos, box_geo);

      p3m_add_pair_force(q, d, d.norm(), p.f.f);
    }
    if (p.r.p[2] > (elc_params.h - elc_params.space_layer)) {
      q = elc_params.delta_mid_top * p.p.q * p.p.q;
      pos[0] = p.r.p[0];
      pos[1] = p.r.p[1];
      pos[2] = 2 * elc_params.h - p.r.p[2];
      auto const d = get_mi_vector(p.r.p, pos, box_geo);

      p3m_add_pair_force(q, d, d.norm(), p.f.f);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////

namespace {
void assign_image_charge(const Particle &p) {
  if (p.r.p[2] < elc_params.space_layer) {
    auto const q_eff = elc_params.delta_mid_bot * p.p.q;
    auto const pos = Utils::Vector3d{p.r.p[0], p.r.p[1], -p.r.p[2]};

    p3m_assign_charge(q_eff, pos);
  }

  if (p.r.p[2] > (elc_params.h - elc_params.space_layer)) {
    auto const q_eff = elc_params.delta_mid_top * p.p.q;
    auto const pos =
        Utils::Vector3d{p.r.p[0], p.r.p[1], 2 * elc_params.h - p.r.p[2]};

    p3m_assign_charge(q_eff, pos);
  }
}
} // namespace

void ELC_p3m_charge_assign_both(const ParticleRange &particles) {
  p3m.inter_weights.reset(p3m.params.cao);

  /* prepare local FFT mesh */
  for (int i = 0; i < p3m.local_mesh.size; i++)
    p3m.rs_mesh[i] = 0.0;

  for (auto &p : particles) {
    if (p.p.q != 0.0) {
      p3m_assign_charge(p.p.q, p.r.p, p3m.inter_weights);
      assign_image_charge(p);
    }
  }
}

void ELC_p3m_charge_assign_image(const ParticleRange &particles) {
  /* prepare local FFT mesh */
  for (int i = 0; i < p3m.local_mesh.size; i++)
    p3m.rs_mesh[i] = 0.0;

  for (auto &p : particles) {
    if (p.p.q != 0.0) {
      assign_image_charge(p);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////

void ELC_P3M_dielectric_layers_force_contribution(Particle const &p1,
                                                  Particle const &p2,
                                                  Utils::Vector3d &force1,
                                                  Utils::Vector3d &force2) {
  Utils::Vector3d pos;
  double q;

  if (p1.r.p[2] < elc_params.space_layer) {
    q = elc_params.delta_mid_bot * p1.p.q * p2.p.q;
    pos[0] = p1.r.p[0];
    pos[1] = p1.r.p[1];
    pos[2] = -p1.r.p[2];
    auto const d = get_mi_vector(p2.r.p, pos, box_geo);

    p3m_add_pair_force(q, d, d.norm(), force2);
  }

  if (p1.r.p[2] > (elc_params.h - elc_params.space_layer)) {
    q = elc_params.delta_mid_top * p1.p.q * p2.p.q;
    pos[0] = p1.r.p[0];
    pos[1] = p1.r.p[1];
    pos[2] = 2 * elc_params.h - p1.r.p[2];
    auto const d = get_mi_vector(p2.r.p, pos, box_geo);

    p3m_add_pair_force(q, d, d.norm(), force2);
  }

  if (p2.r.p[2] < elc_params.space_layer) {
    q = elc_params.delta_mid_bot * p1.p.q * p2.p.q;
    pos[0] = p2.r.p[0];
    pos[1] = p2.r.p[1];
    pos[2] = -p2.r.p[2];
    auto const d = get_mi_vector(p1.r.p, pos, box_geo);

    p3m_add_pair_force(q, d, d.norm(), force1);
  }

  if (p2.r.p[2] > (elc_params.h - elc_params.space_layer)) {
    q = elc_params.delta_mid_top * p1.p.q * p2.p.q;
    pos[0] = p2.r.p[0];
    pos[1] = p2.r.p[1];
    pos[2] = 2 * elc_params.h - p2.r.p[2];
    auto const d = get_mi_vector(p1.r.p, pos, box_geo);

    p3m_add_pair_force(q, d, d.norm(), force1);
  }
}

/////////////////////////////////////////////////////////////////////////////////////

double ELC_P3M_dielectric_layers_energy_contribution(Particle const &p1,
                                                     Particle const &p2) {
  Utils::Vector3d pos;
  double q;
  double tp2;
  double eng = 0.0;

  tp2 = p2.r.p[2];

  if (p1.r.p[2] < elc_params.space_layer) {
    q = elc_params.delta_mid_bot * p1.p.q * p2.p.q;
    pos[0] = p1.r.p[0];
    pos[1] = p1.r.p[1];
    pos[2] = -p1.r.p[2];

    eng += p3m_pair_energy(q, get_mi_vector(p2.r.p, pos, box_geo).norm());
  }

  if (p1.r.p[2] > (elc_params.h - elc_params.space_layer)) {
    q = elc_params.delta_mid_top * p1.p.q * p2.p.q;
    pos[0] = p1.r.p[0];
    pos[1] = p1.r.p[1];
    pos[2] = 2 * elc_params.h - p1.r.p[2];

    eng += p3m_pair_energy(q, get_mi_vector(p2.r.p, pos, box_geo).norm());
  }

  if (tp2 < elc_params.space_layer) {
    q = elc_params.delta_mid_bot * p1.p.q * p2.p.q;
    pos[0] = p2.r.p[0];
    pos[1] = p2.r.p[1];
    pos[2] = -tp2;

    eng += p3m_pair_energy(q, get_mi_vector(p1.r.p, pos, box_geo).norm());
  }

  if (tp2 > (elc_params.h - elc_params.space_layer)) {
    q = elc_params.delta_mid_top * p1.p.q * p2.p.q;
    pos[0] = p2.r.p[0];
    pos[1] = p2.r.p[1];
    pos[2] = 2 * elc_params.h - tp2;

    eng += p3m_pair_energy(q, get_mi_vector(p1.r.p, pos, box_geo).norm());
  }

  return (eng);
}

//////////////////////////////////////////////////////////////////////////////////

double ELC_P3M_dielectric_layers_energy_self(ParticleRange const &particles) {
  Utils::Vector3d pos;
  double q;
  double eng = 0.0;

  // Loop cell neighbors
  for (auto const &p : particles) {
    // Loop neighbor cell particles

    if (p.r.p[2] < elc_params.space_layer) {
      q = elc_params.delta_mid_bot * p.p.q * p.p.q;
      pos[0] = p.r.p[0];
      pos[1] = p.r.p[1];
      pos[2] = -p.r.p[2];

      eng += p3m_pair_energy(q, get_mi_vector(p.r.p, pos, box_geo).norm());
    }

    if (p.r.p[2] > (elc_params.h - elc_params.space_layer)) {
      q = elc_params.delta_mid_top * p.p.q * p.p.q;
      pos[0] = p.r.p[0];
      pos[1] = p.r.p[1];
      pos[2] = 2 * elc_params.h - p.r.p[2];

      eng += p3m_pair_energy(q, get_mi_vector(p.r.p, pos, box_geo).norm());
    }
  }
  return (eng);
}

/////////////////////////////////////////////////////////////////////////////////

void ELC_P3M_modify_p3m_sums_both(ParticleRange const &particles) {
  double node_sums[3], tot_sums[3];

  for (int i = 0; i < 3; i++) {
    node_sums[i] = 0.0;
    tot_sums[i] = 0.0;
  }

  for (auto const &p : particles) {
    if (p.p.q != 0.0) {

      node_sums[0] += 1.0;
      node_sums[1] += Utils::sqr(p.p.q);
      node_sums[2] += p.p.q;

      if (p.r.p[2] < elc_params.space_layer) {

        node_sums[0] += 1.0;
        node_sums[1] += Utils::sqr(elc_params.delta_mid_bot * p.p.q);
        node_sums[2] += elc_params.delta_mid_bot * p.p.q;
      }
      if (p.r.p[2] > (elc_params.h - elc_params.space_layer)) {

        node_sums[0] += 1.0;
        node_sums[1] += Utils::sqr(elc_params.delta_mid_top * p.p.q);
        node_sums[2] += elc_params.delta_mid_top * p.p.q;
      }
    }
  }

  MPI_Allreduce(node_sums, tot_sums, 3, MPI_DOUBLE, MPI_SUM, comm_cart);
  p3m.sum_qpart = (int)(tot_sums[0] + 0.1);
  p3m.sum_q2 = tot_sums[1];
  p3m.square_sum_q = Utils::sqr(tot_sums[2]);
}

void ELC_P3M_modify_p3m_sums_image(ParticleRange const &particles) {
  double node_sums[3], tot_sums[3];

  for (int i = 0; i < 3; i++) {
    node_sums[i] = 0.0;
    tot_sums[i] = 0.0;
  }

  for (auto const &p : particles) {
    if (p.p.q != 0.0) {

      if (p.r.p[2] < elc_params.space_layer) {

        node_sums[0] += 1.0;
        node_sums[1] += Utils::sqr(elc_params.delta_mid_bot * p.p.q);
        node_sums[2] += elc_params.delta_mid_bot * p.p.q;
      }
      if (p.r.p[2] > (elc_params.h - elc_params.space_layer)) {

        node_sums[0] += 1.0;
        node_sums[1] += Utils::sqr(elc_params.delta_mid_top * p.p.q);
        node_sums[2] += elc_params.delta_mid_top * p.p.q;
      }
    }
  }

  MPI_Allreduce(node_sums, tot_sums, 3, MPI_DOUBLE, MPI_SUM, comm_cart);

  p3m.sum_qpart = (int)(tot_sums[0] + 0.1);
  p3m.sum_q2 = tot_sums[1];
  p3m.square_sum_q = Utils::sqr(tot_sums[2]);
}

// this function is required in force.cpp for energy evaluation
void ELC_P3M_restore_p3m_sums(ParticleRange const &particles) {
  double node_sums[3], tot_sums[3];

  for (int i = 0; i < 3; i++) {
    node_sums[i] = 0.0;
    tot_sums[i] = 0.0;
  }

  for (auto const &p : particles) {
    if (p.p.q != 0.0) {

      node_sums[0] += 1.0;
      node_sums[1] += Utils::sqr(p.p.q);
      node_sums[2] += p.p.q;
    }
  }

  MPI_Allreduce(node_sums, tot_sums, 3, MPI_DOUBLE, MPI_SUM, comm_cart);

  p3m.sum_qpart = (int)(tot_sums[0] + 0.1);
  p3m.sum_q2 = tot_sums[1];
  p3m.square_sum_q = Utils::sqr(tot_sums[2]);
}

#endif
