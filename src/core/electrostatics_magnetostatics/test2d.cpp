#include "ParticleRange.hpp"
#include "communication.hpp"
#include "common.hpp"
#include "coulomb.hpp"
#include "grid.hpp"
#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/range/numeric.hpp>

#include "electrostatics_magnetostatics/p3m.hpp"
#include "electrostatics_magnetostatics/test2d.hpp"
#include "errorhandling.hpp"

#ifdef P3M

TEST2DParameters test2d_params{0., 0., 0.};

void TEST2D_set_params(double potential_difference) {
  test2d_params.potential_difference = potential_difference;

  TEST2D_init();

  mpi_bcast_coulomb_params();
}

void TEST2D_init() {
  coulomb.method = COULOMB_TEST2D_P3M;

  TEST2D_on_boxl_change();

  p3m.params.additional_mesh[0] = 0;
  p3m.params.additional_mesh[1] = 0;
  p3m.params.additional_mesh[2] = 0.5 * box_geo.length()[2];
}

void TEST2D_on_boxl_change() {
  test2d_params.gap_half = 0.25 * box_geo.length()[2];
  // the prefactor 2 is due to the walls placed at 0 and box_l/2
  test2d_params.potential_difference_per_boxl = 2 * test2d_params.potential_difference / box_geo.length()[2];
}

void TEST2D_disable() {
  coulomb.method = COULOMB_P3M;

  for (int i = 0; i < 3; ++i) {
    p3m.params.additional_mesh[i] = 0;
  }
}

int TEST2D_sanity_check() {
  if (!box_geo.periodic(0) or !box_geo.periodic(1) or !box_geo.periodic(2)) {
    runtimeErrorMsg() << "TEST2D requires periodicity 1 1 1";
    return ES_ERROR;
  }

  return ES_OK;
}

void TEST2D_modify_P3M_sums(const ParticleRange &particles) {
  Utils::Vector2d moments{};

  for (auto const &p : particles) {
    if (p.p.q != 0) {
      moments[0] += 1;
      moments[1] += p.p.q * p.p.q;
    }
  }

  moments = boost::mpi::all_reduce(comm_cart, moments, std::plus<>());

  // mirror charges double the amount of charges in the system
  p3m.sum_qpart = 2 * static_cast<int>(moments[0]);
  // the total charge in the system has to vanish by construction
  p3m.square_sum_q = 0;
  // mirrored charges have the same q^2 value
  p3m.sum_q2 = 2 * moments[1];
}

void TEST2D_assign_image_charge(const Particle &p) {
  if (p.r.p[2] <= test2d_params.gap_half) {
    auto const pos = detail::mirror_position(p.r.p, 0);

    p3m_assign_charge(-p.p.q, pos);
  } else {
    auto const pos = detail::mirror_position(p.r.p, 0.5 * box_geo.length()[2]);

    p3m_assign_charge(-p.p.q, pos);
  }
}

void TEST2D_charge_assign_P3M(const ParticleRange &particles) {
  p3m.inter_weights.reset(p3m.params.cao);

  /* prepare local FFT mesh */
  for (int i = 0; i < p3m.local_mesh.size; i++) {
    p3m.rs_mesh[i] = 0.0;
  }

  for (auto const &p : particles) {
    if (p.p.q != 0.0) {
      p3m_assign_charge(p.p.q, p.r.p, p3m.inter_weights);
      TEST2D_assign_image_charge(p);
    }
  }
}

double TEST2D_image_charge_energy_contribution(const Utils::Vector3d &pos1,
                                         const Utils::Vector3d &pos2,
                                         double q1q2) {

  // bottom image charge
  if (pos1[2] <= test2d_params.gap_half) {
    return p3m_pair_energy(
        -q1q2,
        get_mi_vector(pos2, detail::mirror_position(pos1, 0), box_geo).norm());
  }

  // top image charge
  return p3m_pair_energy(
      -q1q2,
      get_mi_vector(pos2, detail::mirror_position(pos1, 0.5 * box_geo.length()[2]), box_geo)
          .norm());
}

double TEST2D_image_charge_energy(Particle const &p1, Particle const &p2) {
  auto const pos1 = p1.r.p;
  auto const pos2 = p2.r.p;
  auto const q1q2 = p1.p.q * p2.p.q;

  return TEST2D_image_charge_energy_contribution(pos1, pos2, q1q2)
         + TEST2D_image_charge_energy_contribution(pos2, pos1, q1q2);
}

double TEST2D_potential_difference_energy(double z_pos, double charge) {
  return charge * z_pos * test2d_params.potential_difference_per_boxl;
}

double TEST2D_image_charge_energy(const ParticleRange &particles) {
  return boost::accumulate(particles, 0.0, [](double energy, auto const &p) {
        return energy + 0.5 * coulomb.prefactor * TEST2D_image_charge_energy_contribution(p.r.p, p.r.p, p.p.q * p.p.q) + TEST2D_potential_difference_energy(p.r.p[2], p.p.q);
      });
}

Utils::Vector3d TEST2D_image_charge_force_contribution(const Utils::Vector3d &pos1,
                                                 const Utils::Vector3d &pos2,
                                                 double q1q2) {
  Utils::Vector3d force{};
  auto const q = -q1q2;
  if (pos1[2] <= test2d_params.gap_half) {
    // bottom image charge
    auto const d = get_mi_vector(pos2, detail::mirror_position(pos1, 0), box_geo);

    p3m_add_pair_force(q, d, d.norm(), force);
  } else {
    // top image charge
    auto const d = get_mi_vector(
        pos2, detail::mirror_position(pos1, 0.5 * box_geo.length()[2]), box_geo);

    p3m_add_pair_force(q, d, d.norm(), force);
  }

  return force;
}

Utils::Vector3d TEST2D_potential_difference_force(double charge) {
  return {0, 0, -charge * test2d_params.potential_difference_per_boxl};
}

void TEST2D_image_charge_force_contribution(const ParticleRange &particles) {
  for (auto &p : particles) {
    p.f.f += coulomb.prefactor * TEST2D_image_charge_force_contribution(p.r.p, p.r.p, p.p.q * p.p.q)
             + TEST2D_potential_difference_force(p.p.q);
  }
}

#endif // P3M