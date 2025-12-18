#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>

#include "Empiric.hpp"
#include "Estimate.hpp"
#include "MainDist.hpp"

struct ResultRow {
  ld c;
  ld mae_clean;
  ld mae_poison;
};

static ld meanAbsError(const nstu::vector& errs) {
  ld sum = 0;
  for (ld e : errs) sum += std::fabs(e);
  return sum / errs.size();
}

static uint32_t min_u32(uint32_t a, uint32_t b) { return (a < b) ? a : b; }

static void poison(Empiric& em, uint32_t poison_count, MainDist& contam) {
  const uint32_t n = em.size();
  poison_count = min_u32(poison_count, n);
  for (uint32_t i = 0; i < poison_count; ++i) {
    em.changeData(i, contam.Xi());
  }
  em.notify();
}

static void printSmallestWeights(const Empiric& em, const nstu::vector& ws,
                                 int k) {
  const ld inf = 1e300L;
  ld* best_w = new ld[k];
  ld* best_x = new ld[k];
  for (int i = 0; i < k; ++i) {
    best_w[i] = inf;
    best_x[i] = 0;
  }

  for (uint32_t i = 0; i < em.size(); ++i) {
    const ld x = em.at(i);
    const ld w = ws[i];
    if (w >= best_w[k - 1]) continue;
    int pos = k - 1;
    while (pos > 0 && w < best_w[pos - 1]) {
      best_w[pos] = best_w[pos - 1];
      best_x[pos] = best_x[pos - 1];
      --pos;
    }
    best_w[pos] = w;
    best_x[pos] = x;
  }

  for (int i = 0; i < k; ++i) {
    std::cout << best_x[i] << "\t" << best_w[i] << "\n";
  }

  delete[] best_w;
  delete[] best_x;
}

static void printLargestWeights(const Empiric& em, const nstu::vector& ws,
                                int k) {
  const ld neg_inf = -1e300L;
  ld* best_w = new ld[k];
  ld* best_x = new ld[k];
  for (int i = 0; i < k; ++i) {
    best_w[i] = neg_inf;
    best_x[i] = 0;
  }

  for (uint32_t i = 0; i < em.size(); ++i) {
    const ld x = em.at(i);
    const ld w = ws[i];
    if (w <= best_w[k - 1]) continue;
    int pos = k - 1;
    while (pos > 0 && w > best_w[pos - 1]) {
      best_w[pos] = best_w[pos - 1];
      best_x[pos] = best_x[pos - 1];
      --pos;
    }
    best_w[pos] = w;
    best_x[pos] = x;
  }

  for (int i = 0; i < k; ++i) {
    std::cout << best_x[i] << "\t" << best_w[i] << "\n";
  }

  delete[] best_w;
  delete[] best_x;
}

int main() {
  std::cout << std::fixed << std::setprecision(6);
  const unsigned base_seed = 12345;

  const uint32_t n = 500;
  const ld poison_p = 0.2L;
  const uint32_t poison_count =
      static_cast<uint32_t>(poison_p * n + 0.5L);
  const int samples = 20;

  const ld nu_close_to_normal = 75.0L;
  const ld mu_true = 0.0L;

  MainDist ideal(nu_close_to_normal, mu_true, 1.0L);
  MainDist contam(nu_close_to_normal, 6.0L, 3.0L);
  const ld sigma = std::sqrt(ideal.D());

  const ld c_values[] = {0.75L, 1.0L, 1.5L, 1.725L, 2.385L};
  const uint32_t c_count = sizeof(c_values) / sizeof(c_values[0]);

  ResultRow rows[5]{};

  std::cout << "n=" << n << ", samples=" << samples << ", poison_p=" << poison_p
            << ", sigma=" << sigma << "\n";
  std::cout << "ideal: nu=" << ideal.getNu() << " mu=" << ideal.getMu()
            << " lambda=" << ideal.getLambda() << "\n";
  std::cout << "contam: nu=" << contam.getNu() << " mu=" << contam.getMu()
            << " lambda=" << contam.getLambda() << "\n\n";

  nstu::vector mean_clean_errs;
  nstu::vector mean_poison_errs;
  for (int s = 0; s < samples; ++s) {
    srand(base_seed + static_cast<unsigned>(s) * 2);
    Empiric clean(ideal.Xi(n));
    mean_clean_errs.push_back(clean.mean() - mu_true);
    srand(base_seed + static_cast<unsigned>(s) * 2);
    Empiric poisoned(ideal.Xi(n));
    srand(base_seed + static_cast<unsigned>(s) * 2 + 1);
    poison(poisoned, poison_count, contam);
    mean_poison_errs.push_back(poisoned.mean() - mu_true);
  }

  for (uint32_t ci = 0; ci < c_count; ++ci) {
    const ld c = c_values[ci];
    nstu::vector clean_errs;
    nstu::vector poison_errs;

    for (int s = 0; s < samples; ++s) {
      srand(base_seed + static_cast<unsigned>(s) * 2);
      Empiric em(ideal.Xi(n));
      Estimate est(em, sigma, c, 50);

      const ld mu_clean = est.mu();
      clean_errs.push_back(mu_clean - mu_true);

      srand(base_seed + static_cast<unsigned>(s) * 2 + 1);
      poison(em, poison_count, contam);
      const ld mu_poison = est.mu();
      poison_errs.push_back(mu_poison - mu_true);
    }

    rows[ci] = {c, meanAbsError(clean_errs), meanAbsError(poison_errs)};
  }

  ResultRow best = rows[0];
  for (uint32_t i = 0; i < c_count; ++i) {
    if (rows[i].mae_poison < best.mae_poison) best = rows[i];
  }

  std::cout << "MAE(|mu_hat - mu_true|)\n";
  std::cout << "c\tclean\t\tpoison\n";
  std::cout << "mean\t" << meanAbsError(mean_clean_errs) << "\t"
            << meanAbsError(mean_poison_errs) << "\n";
  for (uint32_t i = 0; i < c_count; ++i) {
    std::cout << rows[i].c << "\t" << rows[i].mae_clean << "\t"
              << rows[i].mae_poison << "\n";
  }
  std::cout << "\nbest_c=" << best.c << " (poison MAE=" << best.mae_poison
            << ")\n\n";

  {
    srand(base_seed + 999);
    Empiric em(ideal.Xi(n));
    Estimate est(em, sigma, best.c, 80);
    const ld mean_clean = em.mean();
    const ld mu_clean = est.mu();

    srand(base_seed + 1000);
    poison(em, poison_count, contam);
    const ld mean_poison = em.mean();
    const ld mu_poison = est.mu();

    std::cout << "One sample details (best_c=" << best.c << ")\n";
    std::cout << "clean: mean=" << mean_clean << " estimate=" << mu_clean
              << "\n";
    std::cout << "poison: mean=" << mean_poison << " estimate=" << mu_poison
              << "\n";

    const nstu::vector ws = est.weights();
    std::cout << "\n10 smallest weights (x, w)\n";
    printSmallestWeights(em, ws, 10);

    std::cout << "\n10 largest weights (x, w)\n";
    printLargestWeights(em, ws, 10);
  }

  return 0;
}
