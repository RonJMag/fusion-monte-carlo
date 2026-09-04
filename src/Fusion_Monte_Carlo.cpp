// =====================================================================================
//  Fusion_Monte_Carlo.cpp
//
//  Generalization of Final_Monte_Carlo_Code.cpp (e+e- -> mu+mu- event generator)
//  to thermonuclear fusion:  D + T -> n + alpha   (and D + D -> n + 3He)
//
//  The original three-stage pipeline is preserved and renamed:
//
//     monteCarloEstimate  ->  monteCarloReactivity   (Stage 1: <sigma v>(T) by MC integration)
//     HitOrMiss           ->  sampleReactantPairs    (Stage 2: accept/reject on sigma*v)
//     CreateEvents        ->  createEvents           (Stage 3: 2->2 kinematics with masses + boost)
//     (new)               ->  histogram / FWHM       (Stage 4: Doppler-broadened neutron spectrum)
//     (new)               ->  beam-target mode, D-D channel, slab transport (Stage 5)
//
//  Units: energies in keV, masses in keV/c^2, velocities as beta = v/c, cross sections in mb.
//
//  Build (g++/clang):   g++ -O2 -std=c++17 Fusion_Monte_Carlo.cpp -o fusion_mc
//  Build (MSVC):        cl /EHsc /O2 /std:c++17 Fusion_Monte_Carlo.cpp /Fe:fusion_mc.exe
//
//  Run examples:
//     ./fusion_mc                                   # D-T thermal plasma, T = 10 keV, 100000 samples
//     ./fusion_mc --reaction DD --T 10              # D-D neutron branch
//     ./fusion_mc --mode beam --Ebeam 100 --T 10    # 100 keV deuteron beam into 10 keV tritium
//     ./fusion_mc --scan                            # <sigma v>(T) table vs Bosch-Hale fit
//     ./fusion_mc --transport 50                    # add 50 cm toy lithium-slab neutron transport
// =====================================================================================

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <algorithm>

// ------------------------------------------------------------------ constants (keV, c = 1)
constexpr double PI      = 3.14159265358979323846;
constexpr double C_LIGHT = 2.99792458e8;      // m/s
constexpr double MB_TO_M2 = 1.0e-31;          // 1 millibarn = 1e-31 m^2

constexpr double M_P   = 938272.09;           // keV
constexpr double M_N   = 939565.42;
constexpr double M_D   = 1875612.94;
constexpr double M_T   = 2808921.11;
constexpr double M_HE3 = 2808391.61;
constexpr double M_HE4 = 3727379.38;

// ------------------------------------------------------------------ random numbers
// Replaces rand()/RAND_MAX (15-bit on MSVC) with a 64-bit Mersenne Twister.
struct RNG {
    std::mt19937_64 gen;
    std::uniform_real_distribution<double> uni{0.0, 1.0};
    std::normal_distribution<double> gauss{0.0, 1.0};
    explicit RNG(unsigned long long seed) : gen(seed) {}
    double u()   { return uni(gen);   }   // uniform (0,1)
    double n01() { return gauss(gen); }   // standard normal
};

// ------------------------------------------------------------------ reaction data
// Bosch & Hale, Nuclear Fusion 32 (1992) 611.  sigma(E) in millibarn, E = CM energy in keV.
//   sigma(E) = S(E) / ( E * exp(B_G / sqrt(E)) )
//   S(E)     = (A1 + E(A2 + E(A3 + E(A4 + E A5)))) / (1 + E(B1 + E(B2 + E(B3 + E B4))))
struct Reaction {
    std::string name;
    double m1, m2;           // reactants
    double m3, m4;           // products  (m3 = the particle we histogram, e.g. the neutron)
    std::string p3, p4;
    double BG;               // Gamow constant, keV^1/2
    double A[5], B[4];
    double Emin, Emax;       // validity range of the sigma(E) fit, keV
    // Bosch-Hale <sigma v> fit (for validation only): C1..C7, mrc2 (keV), valid T range
    double C[7];
    double mrc2;
    double Tmin, Tmax;
};

Reaction makeDT() {
    Reaction r;
    r.name = "D-T";
    r.m1 = M_D;  r.m2 = M_T;  r.m3 = M_N;  r.m4 = M_HE4;
    r.p3 = "n";  r.p4 = "alpha";
    r.BG = 34.3827;
    double A[5] = {6.927e4, 7.454e8, 2.050e6, 5.2002e4, 0.0};
    double B[4] = {6.38e1, -9.95e-1, 6.981e-5, 1.728e-4};
    std::copy(A, A + 5, r.A);  std::copy(B, B + 4, r.B);
    r.Emin = 0.5;  r.Emax = 550.0;
    double C[7] = {1.17302e-9, 1.51361e-2, 7.51886e-2, 4.60643e-3, 1.35000e-2, -1.06750e-4, 1.36600e-5};
    std::copy(C, C + 7, r.C);
    r.mrc2 = 1124656.0;  r.Tmin = 0.2;  r.Tmax = 100.0;
    return r;
}

Reaction makeDDn() {                        // D + D -> n + 3He  (50 % branch)
    Reaction r;
    r.name = "D-D(n)";
    r.m1 = M_D;  r.m2 = M_D;  r.m3 = M_N;  r.m4 = M_HE3;
    r.p3 = "n";  r.p4 = "He3";
    r.BG = 31.3970;
    double A[5] = {5.3701e4, 3.3027e2, -1.2706e-1, 2.9327e-5, -2.5151e-9};
    double B[4] = {0.0, 0.0, 0.0, 0.0};
    std::copy(A, A + 5, r.A);  std::copy(B, B + 4, r.B);
    r.Emin = 0.5;  r.Emax = 4900.0;
    double C[7] = {5.43360e-12, 5.85778e-3, 7.68222e-3, 0.0, -2.96400e-6, 0.0, 0.0};
    std::copy(C, C + 7, r.C);
    r.mrc2 = 937814.0;  r.Tmin = 0.2;  r.Tmax = 100.0;
    return r;
}

double sigmaBoschHale(const Reaction& r, double E) {  // millibarn
    if (E <= 0.0 || E > r.Emax) return 0.0;
    double num = r.A[0] + E * (r.A[1] + E * (r.A[2] + E * (r.A[3] + E * r.A[4])));
    double den = 1.0 + E * (r.B[0] + E * (r.B[1] + E * (r.B[2] + E * r.B[3])));
    double S = num / den;
    return S / (E * std::exp(r.BG / std::sqrt(E)));
}

// Reference <sigma v>(T) fit from the same paper, m^3/s.  Used only to validate Stage 1.
double reactivityBoschHaleFit(const Reaction& r, double T) {
    double theta = T / (1.0 - T * (r.C[1] + T * (r.C[3] + T * r.C[5]))
                              / (1.0 + T * (r.C[2] + T * (r.C[4] + T * r.C[6]))));
    double xi = std::cbrt(r.BG * r.BG / (4.0 * theta));
    double sv_cm3 = r.C[0] * theta * std::sqrt(xi / (r.mrc2 * T * T * T)) * std::exp(-3.0 * xi);
    return sv_cm3 * 1.0e-6;
}

inline double reducedMass(const Reaction& r) { return r.m1 * r.m2 / (r.m1 + r.m2); }

// sigma(E) * v(E) in m^3/s, E = relative (CM) kinetic energy in keV
inline double sigmaV(const Reaction& r, double E) {
    double beta = std::sqrt(2.0 * E / reducedMass(r));
    return sigmaBoschHale(r, E) * MB_TO_M2 * beta * C_LIGHT;
}

// =====================================================================================
//  STAGE 1  --  monteCarloReactivity   (was: monteCarloEstimate)
//
//  <sigma v>(T) = Integral_0^inf  sigma(E) v(E) f_M(E; T) dE
//  f_M(E;T) = (2/sqrt(pi)) T^{-3/2} sqrt(E) exp(-E/T)   (Maxwellian in relative energy)
//
//  Two estimators are returned so you can see the effect of importance sampling:
//    (a) uniform E in [0, Ecut]          -- exactly the method of the original code
//    (b) E ~ Exp(T) (importance sampled) -- weight = sigma v (2/sqrt(pi)) sqrt(E/T)
// =====================================================================================
struct MCResult { double estimate, error, wmax; };

MCResult monteCarloEstimate(std::function<double(double)> f, double lowBound, double upBound,
                            long iterations, RNG& rng) {
    double totalSum = 0.0, sumSqr = 0.0, wmax = -1.0;
    for (long i = 0; i < iterations; ++i) {
        double x = lowBound + rng.u() * (upBound - lowBound);
        double fv = f(x) * (upBound - lowBound);
        if (fv > wmax) wmax = fv;
        totalSum += fv;
        sumSqr   += fv * fv;
    }
    double est = totalSum / iterations;
    double var = sumSqr / iterations - est * est;
    return { est, std::sqrt(std::max(var, 0.0) / iterations), wmax };
}

MCResult monteCarloReactivityUniform(const Reaction& r, double T, long N, RNG& rng) {
    double Ecut = std::min(40.0 * T, r.Emax);
    auto integrand = [&](double E) {
        double fM = (2.0 / std::sqrt(PI)) * std::pow(T, -1.5) * std::sqrt(E) * std::exp(-E / T);
        return sigmaV(r, E) * fM;
    };
    return monteCarloEstimate(integrand, 0.0, Ecut, N, rng);
}

MCResult monteCarloReactivityImportance(const Reaction& r, double T, long N, RNG& rng) {
    double totalSum = 0.0, sumSqr = 0.0, wmax = -1.0;
    for (long i = 0; i < N; ++i) {
        double E = -T * std::log(rng.u());                    // E ~ exp(-E/T)/T
        double w = sigmaV(r, E) * (2.0 / std::sqrt(PI)) * std::sqrt(E / T);
        if (w > wmax) wmax = w;
        totalSum += w;  sumSqr += w * w;
    }
    double est = totalSum / N;
    double var = sumSqr / N - est * est;
    return { est, std::sqrt(std::max(var, 0.0) / N), wmax };
}

// =====================================================================================
//  STAGE 2  --  sampleReactantPairs   (was: HitOrMiss)
//
//  Draw ion velocities, form the relative energy, accept the pair with probability
//  sigma(E_rel) v_rel / W_max.  Accepted pairs are distributed as the fusion *rate*,
//  which is exactly what the original HitOrMiss did with (1 + x^2)/W_max.
//
//  thermal mode : both ions Maxwellian at temperature T
//  beam   mode  : ion 1 monoenergetic along +z (neutral-beam or accelerator), ion 2 Maxwellian
// =====================================================================================
struct Pair { double b1[3], b2[3]; double Erel; };

struct Vec3 { double x, y, z; };
inline Vec3 sub(const Vec3& a, const Vec3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

// Deterministic W_max: scan sigma*v on a grid (the original code estimated W_max from the
// same 100 samples it integrated with, which biases the tails; a grid scan avoids that).
double findWmax(const Reaction& r, double Ecut) {
    double w = 0.0;
    for (int i = 1; i <= 20000; ++i) w = std::max(w, sigmaV(r, Ecut * i / 20000.0));
    return w * 1.02;                                           // 2 % safety margin
}

std::vector<Pair> sampleReactantPairs(const Reaction& r, double T, long Nevents, bool beamMode,
                                      double Ebeam, RNG& rng, long& trials) {
    std::vector<Pair> out;  out.reserve(Nevents);
    double mu   = reducedMass(r);
    double Ecut = beamMode ? std::min(Ebeam * 3.0 + 40.0 * T, r.Emax) : std::min(40.0 * T, r.Emax);
    double Wmax = findWmax(r, Ecut);
    double s1 = std::sqrt(T / r.m1), s2 = std::sqrt(T / r.m2);  // thermal beta per component
    double betaBeam = std::sqrt(2.0 * Ebeam / r.m1);
    trials = 0;
    while ((long)out.size() < Nevents) {
        ++trials;
        Pair p;
        if (beamMode) { p.b1[0] = 0; p.b1[1] = 0; p.b1[2] = betaBeam; }
        else          { for (int k = 0; k < 3; ++k) p.b1[k] = s1 * rng.n01(); }
        for (int k = 0; k < 3; ++k) p.b2[k] = s2 * rng.n01();
        double dx = p.b1[0] - p.b2[0], dy = p.b1[1] - p.b2[1], dz = p.b1[2] - p.b2[2];
        p.Erel = 0.5 * mu * (dx * dx + dy * dy + dz * dz);
        if (p.Erel > Ecut) continue;                            // negligible Maxwellian tail
        double prob = sigmaV(r, p.Erel) / Wmax;
        if (rng.u() < prob) out.push_back(p);                   // hit
    }
    return out;
}

// =====================================================================================
//  STAGE 3  --  createEvents   (was: CreateEvents)
//
//  Same construction as the original: back-to-back products at (theta, phi) in the CM frame,
//  but now with (i) massive products and a Q-value, (ii) an angular-distribution hook
//  (isotropic for keV fusion; the original used 1 + cos^2), and (iii) a Lorentz boost from
//  the CM frame to the lab frame using the sampled reactant velocities.
// =====================================================================================
struct FourVec { double E, px, py, pz; };

struct Event {
    double Erel;                 // keV, reactant CM kinetic energy
    double cosThetaCM, phi;      // emission direction of particle 3 in CM
    FourVec p3lab, p4lab;        // lab four-vectors (keV)
    double K3lab, K4lab;         // lab kinetic energies
    double cosTheta3lab;         // lab polar angle of particle 3 w.r.t. +z
    double pT3lab;               // lab transverse momentum of particle 3
};

// Angular distribution hook (CM frame).  D-T and D-D at keV energies are s-wave: isotropic.
// For beam energies above ~100s of keV you would return 1 + a2 P2(cos) + ... here, and the
// rejection step below turns back into the original (1 + x^2)/W_max sampler.
double angularWeight(double /*cosTheta*/, double /*Erel*/) { return 1.0; }
constexpr double ANGULAR_WMAX = 1.0;

FourVec boost(const FourVec& p, const Vec3& beta) {
    double b2 = dot(beta, beta);
    if (b2 < 1e-30) return p;
    double gamma = 1.0 / std::sqrt(1.0 - b2);
    double bp = beta.x * p.px + beta.y * p.py + beta.z * p.pz;
    double f = (gamma - 1.0) * bp / b2 + gamma * p.E;
    return { gamma * (p.E + bp), p.px + f * beta.x, p.py + f * beta.y, p.pz + f * beta.z };
}

std::vector<Event> createEvents(const Reaction& r, const std::vector<Pair>& pairs, RNG& rng) {
    std::vector<Event> ev;  ev.reserve(pairs.size());
    for (const Pair& q : pairs) {
        Event e;  e.Erel = q.Erel;
        // invariant mass of the initial state (reactants are non-relativistic)
        double sqrtS = r.m1 + r.m2 + q.Erel;
        double s = sqrtS * sqrtS;
        double a = s - (r.m3 + r.m4) * (r.m3 + r.m4);
        double b = s - (r.m3 - r.m4) * (r.m3 - r.m4);
        double pStar = std::sqrt(a * b) / (2.0 * sqrtS);       // |p*| of each product in CM
        double E3 = std::sqrt(pStar * pStar + r.m3 * r.m3);
        double E4 = std::sqrt(pStar * pStar + r.m4 * r.m4);
        // angular sampling (hit-or-miss on the angular hook, as in the original)
        double cosTH;
        do { cosTH = -1.0 + 2.0 * rng.u(); } while (rng.u() * ANGULAR_WMAX > angularWeight(cosTH, q.Erel));
        double sinTH = std::sqrt(1.0 - cosTH * cosTH);
        double phi = 2.0 * PI * rng.u();
        e.cosThetaCM = cosTH;  e.phi = phi;
        FourVec p3 = { E3,  pStar * sinTH * std::cos(phi),  pStar * sinTH * std::sin(phi),  pStar * cosTH };
        FourVec p4 = { E4, -p3.px, -p3.py, -p3.pz };
        // CM velocity of the reactant pair (non-relativistic ions)
        Vec3 bcm = { (r.m1 * q.b1[0] + r.m2 * q.b2[0]) / (r.m1 + r.m2),
                     (r.m1 * q.b1[1] + r.m2 * q.b2[1]) / (r.m1 + r.m2),
                     (r.m1 * q.b1[2] + r.m2 * q.b2[2]) / (r.m1 + r.m2) };
        e.p3lab = boost(p3, bcm);  e.p4lab = boost(p4, bcm);
        e.K3lab = e.p3lab.E - r.m3;  e.K4lab = e.p4lab.E - r.m4;
        double pmag = std::sqrt(e.p3lab.px * e.p3lab.px + e.p3lab.py * e.p3lab.py + e.p3lab.pz * e.p3lab.pz);
        e.cosTheta3lab = e.p3lab.pz / pmag;
        e.pT3lab = std::sqrt(e.p3lab.px * e.p3lab.px + e.p3lab.py * e.p3lab.py);
        ev.push_back(e);
    }
    return ev;
}

// =====================================================================================
//  STAGE 4  --  histogram + FWHM of the lab-frame neutron spectrum
// =====================================================================================
struct Hist {
    double lo, hi;  int nb;  std::vector<long> c;
    Hist(double l, double h, int n) : lo(l), hi(h), nb(n), c(n, 0) {}
    void fill(double x) { if (x >= lo && x < hi) ++c[(int)((x - lo) / (hi - lo) * nb)]; }
    double center(int i) const { return lo + (i + 0.5) * (hi - lo) / nb; }
    void write(const std::string& fn, const char* label) const {
        std::ofstream f(fn);  f << label << ",count\n";
        for (int i = 0; i < nb; ++i) f << center(i) << "," << c[i] << "\n";
    }
    double fwhm() const {   // linear interpolation at half maximum
        long mx = *std::max_element(c.begin(), c.end());  int im = (int)(std::max_element(c.begin(), c.end()) - c.begin());
        double half = 0.5 * mx, w = (hi - lo) / nb;
        int i = im; while (i > 0 && c[i] > half) --i;
        double xl = center(i) + (half - c[i]) / (c[i + 1] - c[i]) * w;
        i = im; while (i < nb - 1 && c[i] > half) ++i;
        double xr = center(i - 1) + (c[i - 1] - half) / (c[i - 1] - c[i]) * w;
        return xr - xl;
    }
};

// =====================================================================================
//  STAGE 5c --  toy neutron transport through a lithium slab (exponential free paths).
//  Constant macroscopic cross sections are placeholders for ENDF/B data; the sampling
//  pattern (distance = -ln(u)/Sigma_t, then absorb/scatter by rejection) is the point.
// =====================================================================================
struct TransportSummary { long transmitted = 0, bred = 0, reflected = 0; double meanCollisions = 0; };

TransportSummary transportSlab(const std::vector<Event>& ev, double thickness_cm, RNG& rng, Hist& hEsc) {
    // natural lithium: n = 4.63e22 /cm^3, sigma_total ~ 1.45 b (14 MeV) -> Sigma_t ~ 0.067 /cm
    const double Sigma_t = 0.067;      // 1/cm   (toy, energy independent)
    const double P_breed = 0.20;       // fraction of collisions that are 6Li(n,a)T or 7Li(n,n'a)T  (toy)
    const double A = 7.0;              // target mass number for elastic energy loss
    TransportSummary s;  long collTot = 0;
    for (const Event& e : ev) {
        double x = 0.0, K = e.K3lab, mu = e.cosTheta3lab;  // enter slab along the lab direction (+z normal)
        if (mu <= 0) mu = -mu;                              // slab surrounds the plasma: use |cos|
        bool alive = true;  int ncol = 0;
        while (alive) {
            double d = -std::log(rng.u()) / Sigma_t;
            x += d * mu;
            if (x >= thickness_cm) { ++s.transmitted; hEsc.fill(K / 1000.0); break; }
            if (x < 0.0)           { ++s.reflected;   break; }
            ++ncol;
            if (rng.u() < P_breed) { ++s.bred; alive = false; break; }
            double cosCM = -1.0 + 2.0 * rng.u();            // isotropic elastic scatter in CM
            K *= (A * A + 1.0 + 2.0 * A * cosCM) / ((A + 1.0) * (A + 1.0));
            mu = -1.0 + 2.0 * rng.u();                      // toy: isotropic lab direction after scatter
        }
        collTot += ncol;
    }
    s.meanCollisions = (double)collTot / ev.size();
    return s;
}

// =====================================================================================
//  main
// =====================================================================================
int main(int argc, char** argv) {
    std::string reaction = "DT", mode = "thermal", prefix = "out";
    double T = 10.0, Ebeam = 100.0, slab = 0.0;
    long N = 100000;
    bool scan = false;
    unsigned long long seed = 20260904ULL;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() { return std::string(argv[++i]); };
        if      (a == "--reaction")  reaction = next();
        else if (a == "--T")         T = std::stod(next());
        else if (a == "--N")         N = std::stol(next());
        else if (a == "--mode")      mode = next();
        else if (a == "--Ebeam")     Ebeam = std::stod(next());
        else if (a == "--out")       prefix = next();
        else if (a == "--seed")      seed = std::stoull(next());
        else if (a == "--scan")      scan = true;
        else if (a == "--transport") slab = std::stod(next());
        else { std::cerr << "unknown option " << a << "\n"; return 1; }
    }
    Reaction r = (reaction == "DD") ? makeDDn() : makeDT();
    bool beamMode = (mode == "beam");
    RNG rng(seed);

    printf("==============================================================\n");
    printf(" Fusion Monte Carlo  |  reaction %s  |  mode %s  |  N = %ld\n", r.name.c_str(), mode.c_str(), N);
    printf("==============================================================\n");

    // ---------------- Stage 1: reactivity
    if (scan) {
        printf("\n[Stage 1] <sigma v>(T) scan, %ld samples per point\n", N);
        printf(" %8s %14s %14s %14s %10s\n", "T[keV]", "MC-uniform", "MC-importance", "Bosch-Hale", "ratio(imp)");
        double Ts[] = {1, 2, 5, 8, 10, 15, 20, 30, 50, 70, 100};
        std::ofstream f(prefix + "_reactivity_scan.csv");
        f << "T_keV,sv_uniform_m3s,err_uniform,sv_importance_m3s,err_importance,sv_boschhale_m3s\n";
        for (double t : Ts) {
            MCResult u = monteCarloReactivityUniform(r, t, N, rng);
            MCResult im = monteCarloReactivityImportance(r, t, N, rng);
            double bh = reactivityBoschHaleFit(r, t);
            printf(" %8.1f %14.4e %14.4e %14.4e %10.4f\n", t, u.estimate, im.estimate, bh, im.estimate / bh);
            f << t << "," << u.estimate << "," << u.error << "," << im.estimate << "," << im.error << "," << bh << "\n";
        }
    }
    MCResult svU = monteCarloReactivityUniform(r, T, N, rng);
    MCResult svI = monteCarloReactivityImportance(r, T, N, rng);
    double svBH = reactivityBoschHaleFit(r, T);
    printf("\n[Stage 1] Reactivity at T = %.2f keV (%ld samples)\n", T, N);
    printf("  uniform sampling     : <sigma v> = %.4e +/- %.2e m^3/s  (rel. err %.2f %%)\n", svU.estimate, svU.error, 100 * svU.error / svU.estimate);
    printf("  importance sampling  : <sigma v> = %.4e +/- %.2e m^3/s  (rel. err %.2f %%)\n", svI.estimate, svI.error, 100 * svI.error / svI.estimate);
    printf("  Bosch-Hale fit       : <sigma v> = %.4e m^3/s   -> MC/fit = %.4f\n", svBH, svI.estimate / svBH);

    // ---------------- Stage 2: reactant pairs
    long trials = 0;
    std::vector<Pair> pairs = sampleReactantPairs(r, T, N, beamMode, Ebeam, rng, trials);
    double meanErel = 0; for (auto& p : pairs) meanErel += p.Erel; meanErel /= pairs.size();
    printf("\n[Stage 2] Hit-or-miss on sigma*v: %ld accepted / %ld trials (efficiency %.2f %%)\n", (long)pairs.size(), trials, 100.0 * pairs.size() / trials);
    printf("  mean reacting E_rel = %.2f keV   (Gamow peak; compare 1.5 T = %.2f keV for all pairs)\n", meanErel, 1.5 * T);

    // ---------------- Stage 3: kinematics
    std::vector<Event> ev = createEvents(r, pairs, rng);
    double Q = r.m1 + r.m2 - r.m3 - r.m4;
    double K30;
    {   // exact two-body kinetic energy of particle 3 at E_rel = 0
        double sqrtS = r.m1 + r.m2, s = sqrtS * sqrtS;
        double pS = std::sqrt((s - (r.m3 + r.m4) * (r.m3 + r.m4)) * (s - (r.m3 - r.m4) * (r.m3 - r.m4))) / (2 * sqrtS);
        K30 = std::sqrt(pS * pS + r.m3 * r.m3) - r.m3;
    }
    printf("\n[Stage 3] Q-value = %.1f keV ; %s kinetic energy at rest = %.1f keV ; %s = %.1f keV\n", Q, r.p3.c_str(), K30, r.p4.c_str(), Q - K30);
    printf("  first 3 events (lab frame, keV):\n");
    for (int i = 0; i < 3 && i < (int)ev.size(); ++i) {
        const Event& e = ev[i];
        printf("   Event %d  E_rel=%.2f  cosTH_cm=%+.3f  P_%s=(%.1f, %.1f, %.1f, %.1f)  P_%s=(%.1f, %.1f, %.1f, %.1f)  K_%s=%.1f\n",
               i, e.Erel, e.cosThetaCM, r.p3.c_str(), e.p3lab.E, e.p3lab.px, e.p3lab.py, e.p3lab.pz,
               r.p4.c_str(), e.p4lab.E, e.p4lab.px, e.p4lab.py, e.p4lab.pz, r.p3.c_str(), e.K3lab);
    }

    // ---------------- Stage 4: neutron spectrum
    double mean = 0, m2 = 0;
    for (auto& e : ev) { mean += e.K3lab; }  mean /= ev.size();
    for (auto& e : ev) { m2 += (e.K3lab - mean) * (e.K3lab - mean); }  double sd = std::sqrt(m2 / ev.size());
    double lo = K30 - 8 * std::max(sd, 1.0), hi = K30 + 8 * std::max(sd, 1.0);
    if (beamMode) { lo = K30 - 3000; hi = K30 + 3000; }
    Hist hE(lo, hi, 200), hCos(-1, 1, 40), hPT(0, 1.0, 100), hErel(0, std::min(40.0 * T + (beamMode ? Ebeam * 2 : 0), r.Emax), 100);
    double pStarRef = std::sqrt(K30 * K30 + 2 * K30 * r.m3);
    for (auto& e : ev) { hE.fill(e.K3lab); hCos.fill(e.cosTheta3lab); hPT.fill(e.pT3lab / pStarRef); hErel.fill(e.Erel); }
    double fwhm = hE.fwhm();
    // Brysk (1973): FWHM = sqrt(8 ln2 * 2 m3 K30 T / (m1+m2))
    double brysk = std::sqrt(8.0 * std::log(2.0) * 2.0 * r.m3 * K30 * T / (r.m1 + r.m2));
    printf("\n[Stage 4] Lab-frame %s spectrum (%ld events)\n", r.p3.c_str(), (long)ev.size());
    printf("  mean K_%s = %.1f keV  (shift from rest value: %+.1f keV)\n", r.p3.c_str(), mean, mean - K30);
    printf("  std dev   = %.1f keV ;  FWHM (histogram) = %.1f keV ; 2.355*sd = %.1f keV\n", sd, fwhm, 2.3548 * sd);
    if (!beamMode) printf("  Brysk thermal-broadening prediction: FWHM = %.1f keV  -> ratio MC/Brysk = %.3f\n", brysk, fwhm / brysk);
    else           printf("  (beam mode: spectrum is set by emission angle, not by T; see cos(theta) vs K correlation)\n");

    hE.write(prefix + "_spectrum_hist.csv", "K_keV");
    hCos.write(prefix + "_costheta_lab_hist.csv", "cos_theta_lab");
    hPT.write(prefix + "_pT_over_pstar_hist.csv", "pT_over_pstar");
    hErel.write(prefix + "_Erel_hist.csv", "E_rel_keV");
    {
        std::ofstream f(prefix + "_events.csv");
        f << "E_rel_keV,cosTheta_cm,phi,K3_lab_keV,K4_lab_keV,cosTheta3_lab,pT3_lab_keV,E3,px3,py3,pz3,E4,px4,py4,pz4\n";
        for (auto& e : ev)
            f << e.Erel << "," << e.cosThetaCM << "," << e.phi << "," << e.K3lab << "," << e.K4lab << "," << e.cosTheta3lab << "," << e.pT3lab << ","
              << e.p3lab.E << "," << e.p3lab.px << "," << e.p3lab.py << "," << e.p3lab.pz << ","
              << e.p4lab.E << "," << e.p4lab.px << "," << e.p4lab.py << "," << e.p4lab.pz << "\n";
    }
    printf("  wrote %s_events.csv and histogram CSVs\n", prefix.c_str());

    // ---------------- Stage 5c: toy transport
    if (slab > 0 && r.p3 == "n") {
        Hist hEsc(0, 16, 64);
        TransportSummary ts = transportSlab(ev, slab, rng, hEsc);
        printf("\n[Stage 5] Toy lithium slab, %.0f cm: transmitted %.1f %%, bred (T-producing) %.1f %%, reflected %.1f %%, mean collisions %.2f\n",
               slab, 100.0 * ts.transmitted / ev.size(), 100.0 * ts.bred / ev.size(), 100.0 * ts.reflected / ev.size(), ts.meanCollisions);
        hEsc.write(prefix + "_transport_escape_hist.csv", "K_MeV");
    }
    printf("\nDone.\n");
    return 0;
}
