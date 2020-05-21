#ifndef ESPRESSO_TEST2D_HPP
#define ESPRESSO_TEST2D_HPP

#ifdef P3M

/** @brief parameters for the test2d method */
typedef struct {
  /** half the gap size to use for mirror charge reflection */
  double gap_half;
} TEST2DParameters;
extern TEST2DParameters test2d_params;

/** @brief set coulomb method and initialize parameters */
void TEST2D_init();

/** @brief calculate necessary parameters */
void TEST2D_on_boxl_change();

/** @brief set coulomb method and reset P3M additional mesh but does not reset
 * mesh or sums */
void TEST2D_disable();

/** @brief sanity check
 * @return error code
 */
int TEST2D_sanity_check();

/**
 * @brief modify the P3M sums to include mirror particles
 * @param particles particles to calculate mirror particles for
 */
void TEST2D_modify_P3M_sums(const ParticleRange &particles);

/**
 * @brief assign charges for real and mirror paricles to p3m mesh
 * @param particles particles to include in mesh calculation
 */
void TEST2D_charge_assign_P3M(const ParticleRange &particles);

/**
 * @brief calculate energy contribution of mirror charges of particle pairs
 * @param p1 first particle
 * @param p2 second particle
 * @return energy contribution of the mirror charges
 */
double TEST2D_pair_energy(Particle const &p1, Particle const &p2);

/**
 * @brief calculate energy contribution of direct mirror charges
 * @param particles particles to calculate energy for
 * @return energy contribution of direct mirror
 */
double TEST2D_self_energy(const ParticleRange &particles);

/**
 * @brief calculate force contribution of direct mirror charges
 * @param particles particles to calculate forces for
 */
void TEST2D_self_forces(const ParticleRange &particles);

/**
 * @brief calculate pairwise force contribution of mirror charges
 * @param pos1 first particle position
 * @param pos2 second particle position
 * @param q1q2 charge factor
 * @return calculated force
 */
Utils::Vector3d TEST2D_force_charge_contribution(const Utils::Vector3d &pos1,
                                                 const Utils::Vector3d &pos2,
                                                 double q1q2);

#endif // P3M
#endif // ESPRESSO_TEST2D_HPP
