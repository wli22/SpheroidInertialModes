#ifndef OBLATEMODE_HPP
#define OBLATEMODE_HPP


#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
#include <tuple>
#include <algorithm>
#include "ALegendre.hpp"


template <typename T>
class OblateMode {
    public:

        OblateMode() = default;

        OblateMode(const int m, const int k, const int delta, const T eccentricity);

        OblateMode(const OblateMode<T>& other)
            : m_(other.m_), m_abs_(other.m_abs_), k_(other.k_), delta_(other.delta_), l_(other.l_), eccentricity_(other.eccentricity_), sigmamnk_(other.sigmamnk_) {}
        
        OblateMode(OblateMode<T>&& other) noexcept
            : m_(other.m_), m_abs_(other.m_abs_), k_(other.k_), delta_(other.delta_), l_(other.l_), eccentricity_(other.eccentricity_), sigmamnk_(std::move(other.sigmamnk_)) {}
        
        OblateMode<T>& operator=(const OblateMode<T>& other) {
            if (this != &other) {
                m_            = other.m_;
                m_abs_        = other.m_abs_;
                k_            = other.k_;
                delta_        = other.delta_;
                l_            = other.l_;
                eccentricity_ = other.eccentricity_;
                sigmamnk_     = other.sigmamnk_;
            }
            return *this;
        }

        OblateMode<T>& operator=(OblateMode<T>&& other) noexcept {
            if (this != &other) {
                m_            = other.m_;
                m_abs_        = other.m_abs_;
                k_            = other.k_;
                delta_        = other.delta_;
                l_            = other.l_;
                eccentricity_ = other.eccentricity_;
                sigmamnk_     = std::move(other.sigmamnk_);
            }
            return *this;
        }

        ~OblateMode() = default;

        const std::vector<T>& half_frequencies() const { return sigmamnk_; }
        T half_frequency(int n) const { return sigmamnk_.at(n-1); }
        T eccentricity_value() const { return eccentricity_; }

        std::pair<T, T> pressure_mode(const int n, const T eta, const T tau, const T phi) const;
        std::pair<T, T> velocity_mode_eta(const int n, const T eta, const T tau, const T phi) const;
        std::pair<T, T> velocity_mode_tau(const int n, const T eta, const T tau, const T phi) const;
        std::pair<T, T> velocity_mode_phi(const int n, const T eta, const T tau, const T phi) const;
        std::tuple<T, T, T, T, T, T> velocity_mode(const int n, const T eta, const T tau, const T phi) const;
        std::pair<T, T> divergence_velocity_mode(const int n, const T eta, const T tau, const T phi) const;

        
    private:
        T sigma_equation(const T sigma) const;
        T bisection_search(const T sigmal, const T sigmar) const;
        void sigma_initialization();
        std::pair<T, T> get_X_and_Y(const T eta, const T tau, const T sigma) const;
        std::tuple<T, T, T> get_terms_123(const T X, const T Y) const;

        int m_;
        int m_abs_;
        int k_;
        int delta_;
        int l_;
        T eccentricity_;
        std::vector<T> sigmamnk_;
        static inline ALegendre<double> alegendre_{};
};


template <typename T>
OblateMode<T>::OblateMode(const int m, const int k, const int delta, const T eccentricity) : m_(m), m_abs_(std::abs(m)), k_(k), delta_(delta), l_(m_abs_+2*k_+delta_), eccentricity_(eccentricity) {
    if (delta_ != 0 && delta_ != 1) {
        std::cerr << "Wrong parameter for equatorial symmetricity!\n";
        std::cerr << "delta must be 0 (symmetric) or 1 (antisymmetric)!\n";
        std::cerr << "The input delta = " << delta_ << "\n";
        std::exit(EXIT_FAILURE);
    }

    if (eccentricity_ < static_cast<T>(0) || eccentricity_ >= static_cast<T>(1)) {
        std::cerr << "Wrong eccentricity!\n";
        std::cerr << "The eccentricity shall be in [0, 1).\n";
        std::exit(EXIT_FAILURE);
    }

    if (m_ == 0 && delta_ == 0 && k_ < 2) {
        std::cerr << "Wrong parameter k for ESAS mode!\n";
        std::cerr << "k must be greater than 1!\n";
        std::cerr << "The input k = " << k_ << "\n";
        std::exit(EXIT_FAILURE);
    }

    if (m_ == 0 && delta_ == 1 && k_ < 1) {
        std::cerr << "Wrong parameter k for EAAS mode!\n";
        std::cerr << "k must be greater than 0!\n";
        std::cerr << "The input k = " << k_ << "\n";
        std::exit(EXIT_FAILURE);
    }

    if (m_ != 0 && delta_ == 0 && k_ < 1) {
        std::cerr << "Wrong parameter k for ESnonAS mode!\n";
        std::cerr << "k must be greater than 0!\n";
        std::cerr << "The input k = " << k_ << "\n";
        std::exit(EXIT_FAILURE);
    }

    if (m_ != 0 && delta_ == 1 && k_ < 0) {
        std::cerr << "Wrong parameter k for EAnonAS mode!\n";
        std::cerr << "k must be nonnegative!\n";
        std::cerr << "The input k = " << k_ << "\n";
        std::exit(EXIT_FAILURE);
    }

    if (m_ == 0) {
        if (delta_ == 0) sigmamnk_.resize(k_-1);
        else sigmamnk_.resize(k_);
    } else {
        if (delta_ == 0) sigmamnk_.resize(2*k_);
        else sigmamnk_.resize(2*k_+1);
    }
    sigma_initialization();
}


template <typename T>
T OblateMode<T>::sigma_equation(const T sigma) const {
    const T sigma2 = sigma*sigma;
    const T ecc2 = eccentricity_*eccentricity_;
    const T one_sub_ecc2 = static_cast<T>(1)-ecc2;
    const T scale = std::sqrt(one_sub_ecc2) / std::sqrt(static_cast<T>(1)-ecc2*sigma2);
    const T xi = sigma*scale;
    auto[plm, result]  = alegendre_.dALegendre_d1_compute(l_, m_, xi);
    // dplms /= std::sqrt(T(1)-xi*xi);
    // result = (static_cast<T>(1)-sigma2)*scale*dplms;
    if (m_ != 0) {
        result *= (T(1)-sigma2)*scale/std::sqrt(T(1)-xi*xi);
        result -= static_cast<T>(m_)*one_sub_ecc2*plm;
    }

    // result /= T(l_);

    return result;
}

template <typename T>
T OblateMode<T>::bisection_search(const T sigmal, const T sigmar) const {
    T xl = sigmal, xr = sigmar, xc = static_cast<T>(0);
    T yl = sigma_equation(xl), yr = sigma_equation(xr), yc = static_cast<T>(0);
    constexpr T tol = 1.0E-15;
    int iter = 0;
    const int Niter = 500;
    T ymax = (std::abs(yl) < std::abs(yr)) ? std::abs(yr) : std::abs(yl);
    if (yl * yr > static_cast<T>(0)) {
        std::cerr << "Bad interval [" << xl << ", " << xr << "] for bisection search.\n";
        std::exit(1);
    }
    xc = (xl + xr) / static_cast<T>(2);
    yc = sigma_equation(xc);
    if (ymax < std::abs(yc)) ymax = std::abs(yc);
    while (std::abs(yc) > tol && (iter++) <= Niter) {
        if (yc * yl <= static_cast<T>(0)) {
            xr = xc;
            yr = yc;
        } else {
            xl = xc;
            yl = yc;
        }
        xc = (xl + xr) / static_cast<T>(2);
        yc = sigma_equation(xc);
    }
    // std::cout << "iter = " << iter << "    yc = " << std::scientific << yc << "\n";
    // if (iter > Niter && std::abs(yc) / static_cast<T>(l) > tol) {
    // if (iter > Niter && std::abs(yc) > tol) {
    //     std::cerr << "DIVERGENCE of bisection search!\n";
    //     std::cout << "yc = " << std::scientific << yc << "\n" << "ymax = " << ymax << "\n" << "ratio = " << std::abs(yc) / static_cast<T>(l) << "\n";
    //     std::exit(1);
    // }
    return xc;
}

template <typename T>
void OblateMode<T>::sigma_initialization() {
    T sigmal = static_cast<T>(0), sigmar = static_cast<T>(0);
    T fl = T(0), fr = T(0), flfr = T(0);
    const T step = (k_ <= 150) ? 0.00001 : 0.000005;
    if (m_ == 0) {
        sigmal = 1.0E-8;
        sigmar = sigmal + step;
        for (T &sigma : sigmamnk_) {
            fl = sigma_equation(sigmal);
            fr = sigma_equation(sigmar);
            flfr = fl*fr;
            while (flfr > static_cast<T>(0)) {
                sigmar += step;
                if (sigmar > static_cast<T>(1)) break;
                fr = sigma_equation(sigmar);
                flfr = fl*fr;
            }
            if (sigmar > static_cast<T>(1)) break;
            sigma = bisection_search(sigmal, sigmar);
            sigmal = sigmar;
            sigmar = sigmal + step;
            if (sigmar > static_cast<T>(1)) break;
            // std::cout << "sigma = " << sigma << "\n";
        }
    } else {
        sigmal = 2.0E-8;
        sigmar = sigmal + step;
        std::size_t i = 0;
        while (sigmal < static_cast<T>(1)) {
            fl = sigma_equation(sigmal);
            fr = sigma_equation(sigmar);
            flfr = fl*fr;
            while (flfr >= static_cast<T>(0)) { 
                sigmar += step;
                if (sigmar >= static_cast<T>(1)) break;
                fr = sigma_equation(sigmar);
                flfr = fl*fr;
            }
            if (sigmar >= static_cast<T>(1)) break;
            if (sigmamnk_.size() == i) break;
            sigmamnk_[i++] = bisection_search(sigmal, sigmar);
            sigmal = sigmar;
            sigmar = sigmal + step;
            // std::cout << "sigmamnk = " << sigmamnk_[i-1] << "  idx = " << i-1 << "\n";
        }
        sigmar = -2.0E-8;
        sigmal = sigmar - step;
        while (sigmar > static_cast<T>(-1)) {
            fl = sigma_equation(sigmal);
            fr = sigma_equation(sigmar);
            flfr = fl*fr;
            while (flfr >= static_cast<T>(0)) {
                sigmal -= step;
                if (sigmal <= static_cast<T>(-1)) break;
                fl = sigma_equation(sigmal);
                flfr = fl*fr;
            }
            if (sigmal <= static_cast<T>(-1)) break;
            if (sigmamnk_.size() == i) break;
            sigmamnk_[i++] = bisection_search(sigmal, sigmar);
            // std::cout << "sigmamnk = " << sigmamnk_[i-1] << "  idx = " << i-1 << "\n";
            sigmar = sigmal;
            sigmal = sigmar - step;
        }
        if (sigmamnk_.size() != i) {
            std::cerr << "Failure of solving half-frequencies!\n";
            std::exit(EXIT_FAILURE);
        }
    }
    if (m_ != 0) std::sort(sigmamnk_.begin(), sigmamnk_.end(), [](const T &a, const T &b) { return std::abs(a) < std::abs(b); } );
}

// template <typename T>
// std::pair<T, T> OblateMode<T>::get_X_and_Y(const T eta, const T tau, const T sigma) const {
//     std::pair<T, T> XY;
//     T sigma2 = sigma*sigma;
//     constexpr T sqrt2 = 1.414213562373095048801688724209698079;
//     T ecc2 = eccentricity_*eccentricity_;
//     T eta2 = eta*eta;
//     T tau2 = tau*tau;
//     T Beta2 = (static_cast<T>(1)-sigma2)/sigma2;
//     T AlphaBeta = std::sqrt(static_cast<T>(1)-ecc2*sigma2) / sigma;
//     T AlphaBeta2 = AlphaBeta*AlphaBeta;
//     T Delta = AlphaBeta2+eta2*tau2-Beta2*(eta2+ecc2)*(static_cast<T>(1)-tau2);
//     T ac4 = static_cast<T>(4)*AlphaBeta2*eta2*tau2;
//     T b2_sub_ac4_sqrt = std::sqrt(Delta*Delta-ac4);
//     T tmp = sqrt2*AlphaBeta;
//     XY.first = std::sqrt(Delta+b2_sub_ac4_sqrt) / tmp;
//     XY.second = (tau/sigma > static_cast<T>(0)) ? std::sqrt(Delta-b2_sub_ac4_sqrt) / tmp : -std::sqrt(Delta-b2_sub_ac4_sqrt) / tmp;
//     if (XY.first > static_cast<T>(1)) XY.first = static_cast<T>(1);
//     else if (XY.first < static_cast<T>(-1))XY.first = static_cast<T>(-1);
//     if (XY.second > static_cast<T>(1)) XY.second = static_cast<T>(1);
//     else if (XY.second < static_cast<T>(-1)) XY.second = static_cast<T>(-1);
//     return XY;
// } 

template <typename T>
std::pair<T, T> OblateMode<T>::get_X_and_Y(const T eta, const T tau, const T sigma) const {
    std::pair<T, T> XY;
    long double sigma2 = sigma*sigma;
    constexpr long double sqrt2 = 1.414213562373095048801688724209698079L;
    long double ecc2 = eccentricity_*eccentricity_;
    long double eta2 = eta*eta;
    long double tau2 = tau*tau;
    long double one_sub_tau2 = 1.0L-tau2;
    if (one_sub_tau2 < 0.0L) one_sub_tau2 = 0.0L;
    long double Beta2 = (1.0L-sigma2)/sigma2;
    long double AlphaBeta = std::sqrt(1.0L-ecc2*sigma2) / sigma;
    long double AlphaBeta2 = AlphaBeta*AlphaBeta;
    long double Delta = AlphaBeta2+eta2*tau2-Beta2*(eta2+ecc2)*one_sub_tau2;
    long double ac4 = 4.0L*AlphaBeta2*eta2*tau2;
    long double b2_sub_ac4_sqrt = std::sqrt(Delta*Delta-ac4);
    long double tmp = sqrt2*AlphaBeta;
    XY.first = std::sqrt(Delta+b2_sub_ac4_sqrt) / tmp;
    XY.second = (tau/sigma > static_cast<T>(0)) ? std::sqrt(Delta-b2_sub_ac4_sqrt) / tmp : -std::sqrt(Delta-b2_sub_ac4_sqrt) / tmp;
    if (XY.first > static_cast<T>(1)) XY.first = static_cast<T>(1);
    else if (XY.first < static_cast<T>(-1))XY.first = static_cast<T>(-1);
    if (XY.second > static_cast<T>(1)) XY.second = static_cast<T>(1);
    else if (XY.second < static_cast<T>(-1)) XY.second = static_cast<T>(-1);
    return XY;
} 

template <typename T>
std::tuple<T, T, T> OblateMode<T>::get_terms_123(const T X, const T Y) const {
    // const auto[plmX, dplmX] = associated_legendre_polynomials_d1(m, l, X);
    // const auto[plmY, dplmY] = associated_legendre_polynomials_d1(m, l, Y);
    const auto[plmX, dplmX_one_sub_X2_sqrt, plmX_div_one_sub_X2_sqrt] = alegendre_.ALegendre_compositions(l_, m_, X);
    const auto[plmY, dplmY_one_sub_Y2_sqrt, plmY_div_one_sub_Y2_sqrt] = alegendre_.ALegendre_compositions(l_, m_, Y);
    T X2 = X*X;
    T Y2 = Y*Y;
    T X2_sub_Y2 = X2 - Y2;
    // T one_sub_X2_sqrt = std::sqrt()
    T one_sub_X2 = static_cast<T>(1)-X2;
    T one_sub_Y2 = static_cast<T>(1)-Y2;
    if (one_sub_X2 < static_cast<T>(0)) one_sub_X2 = static_cast<T>(0);
    if (one_sub_Y2 < static_cast<T>(0)) one_sub_Y2 = static_cast<T>(0);
    T one_sub_X2_sqrt = std::sqrt(one_sub_X2);
    T one_sub_Y2_sqrt = std::sqrt(one_sub_Y2);
    // T one_sub_X2_one_sub_Y2_sqrt = std::sqrt(one_sub_X2*one_sub_Y2);
    // T dplmXplmY = dplmX*plmY;
    // T plmXdplmY = plmX*dplmY;
    std::tuple<T, T, T> terms;
    // std::get<0>(terms) = (X*dplmXplmY-Y*plmXdplmY)*one_sub_X2_one_sub_Y2_sqrt/X2_sub_Y2;
    std::get<0>(terms) = (X*one_sub_Y2_sqrt*dplmX_one_sub_X2_sqrt*plmY - Y*one_sub_X2_sqrt*plmX*dplmY_one_sub_Y2_sqrt)/X2_sub_Y2;
    std::get<1>(terms) = (m_ == 0) ? static_cast<T>(0) : static_cast<T>(m_)*plmX_div_one_sub_X2_sqrt*plmY_div_one_sub_Y2_sqrt;
    std::get<2>(terms) = (Y*one_sub_X2_sqrt*dplmX_one_sub_X2_sqrt*plmY - X*one_sub_Y2_sqrt*plmX*dplmY_one_sub_Y2_sqrt)/X2_sub_Y2;
    if (std::isnan(std::get<0>(terms)) || std::isinf(std::get<0>(terms))) std::get<0>(terms) = static_cast<T>(0);
    if (std::isnan(std::get<1>(terms)) || std::isinf(std::get<1>(terms))) std::get<1>(terms) = static_cast<T>(0);
    if (std::isnan(std::get<2>(terms)) || std::isinf(std::get<2>(terms))) std::get<2>(terms) = static_cast<T>(0);
    return terms;
}

template <typename T>
std::pair<T, T> OblateMode<T>::pressure_mode(const int n, const T eta, const T tau, const T phi) const {
    const T sigma = sigmamnk_.at(n-1);
    const auto[X, Y] = get_X_and_Y(eta, tau, sigma);
    const T plmX = alegendre_.ALegendre_compute(l_, m_, X);
    const T plmY = alegendre_.ALegendre_compute(l_, m_, Y);
    const T plmXplmY = plmX*plmY;
    std::pair<T, T> pmnk_val;
    if (m_ == 0) {
        pmnk_val.first = plmXplmY;
        pmnk_val.second = static_cast<T>(0);
    } else {
        const T mphi = static_cast<T>(m_)*phi;
        pmnk_val.first = plmXplmY*std::cos(mphi);
        pmnk_val.second = plmXplmY*std::sin(mphi);
    }
    return pmnk_val;
}

template <typename T>
std::pair<T, T> OblateMode<T>::velocity_mode_eta(const int n, const T eta, const T tau, const T phi) const {
    const T sigma = sigmamnk_.at(n-1);
    const auto[X, Y] = get_X_and_Y(eta, tau, sigma);
    const auto[term1, term2, term3] = get_terms_123(X, Y);
    std::pair<T, T> umnk;
    T sigma2 = sigma*sigma;
    T kernel;
    T ecc2 = eccentricity_*eccentricity_;
    T tau2 = tau*tau;
    T eta2 = eta*eta;
    T dv = std::sqrt(eta2+ecc2*tau2);
    T one_sub_sigma2 = static_cast<T>(1)-sigma2;
    T one_sub_ecc2sigma2 = static_cast<T>(1)-ecc2*sigma2;
    T one_sub_ecc2sigma2_sqrt = std::sqrt(one_sub_ecc2sigma2);
    T scale = eta*std::sqrt(static_cast<T>(1)-tau2) / (dv*std::sqrt(one_sub_ecc2sigma2*one_sub_sigma2));
    if (std::isnan(scale)) scale = static_cast<T>(0);
    kernel = (m_ == 0) ? scale*sigma*term1 : scale*(sigma*term1-term2);
    scale = tau*std::sqrt(eta2+ecc2) / (dv*one_sub_ecc2sigma2_sqrt);
    if (std::isnan(scale)) scale = static_cast<T>(0);
    kernel -= scale*term3;
    kernel /= static_cast<T>(2);
    if (m_ == 0) {
        umnk.first = static_cast<T>(0);
        umnk.second = kernel;
    } else {
        T mphi = static_cast<T>(m_)*phi;
        umnk.first = -kernel*std::sin(mphi);
        umnk.second = kernel*std::cos(mphi);
    }
    return umnk;
}

template <typename T>
std::pair<T, T> OblateMode<T>::velocity_mode_tau(const int n, const T eta, const T tau, const T phi) const {
    const T sigma = sigmamnk_.at(n-1);
    const auto[X, Y] = get_X_and_Y(eta, tau, sigma);
    const auto[term1, term2, term3] = get_terms_123(X, Y);
    std::pair<T, T> umnk;
    T sigma2 = sigma*sigma;
    T kernel;
    T ecc2 = eccentricity_*eccentricity_;
    T tau2 = tau*tau;
    T eta2 = eta*eta;
    T dv = std::sqrt(eta2+ecc2*tau2);
    T one_sub_sigma2 = static_cast<T>(1)-sigma2;
    T one_sub_ecc2sigma2 = static_cast<T>(1)-ecc2*sigma2;
    T one_sub_ecc2sigma2_sqrt = std::sqrt(one_sub_ecc2sigma2);
    T scale = tau*std::sqrt(eta2+ecc2)/(dv*std::sqrt(one_sub_sigma2*one_sub_ecc2sigma2));
    if (std::isnan(scale)) scale = static_cast<T>(0);
    kernel = (m_ == 0) ? scale*sigma*term1 : scale*(sigma*term1-term2);
    scale = eta*std::sqrt(static_cast<T>(1)-tau2)/(dv*one_sub_ecc2sigma2_sqrt);
    if (std::isnan(scale)) scale = static_cast<T>(0);
    kernel += scale*term3;
    kernel /= static_cast<T>(-2);
    if (m_ == 0) {
        umnk.first = static_cast<T>(0);
        umnk.second = kernel;
    } else {
        T mphi = static_cast<T>(m_)*phi;
        umnk.first = -kernel*std::sin(mphi);
        umnk.second = kernel*std::cos(mphi);
    }
    return umnk;
}

template <typename T>
std::pair<T, T> OblateMode<T>::velocity_mode_phi(const int n, const T eta, const T tau, const T phi) const {
    const T sigma = sigmamnk_.at(n-1);
    const auto[X, Y] = get_X_and_Y(eta, tau, sigma);
    const auto[term1, term2, term3] = get_terms_123(X, Y);
    std::pair<T, T> umnk;
    T sigma2 = sigma*sigma;
    T kernel, scale;
    scale = static_cast<T>(-2)*std::sqrt((static_cast<T>(1)-eccentricity_*eccentricity_*sigma2)*(static_cast<T>(1)-sigma2));
    kernel = (m_ == 0) ? term1/scale : (term1-sigma*term2)/scale;
    if (m_ == 0) {
        umnk.first = kernel;
        umnk.second = static_cast<T>(0);
    } else {
        T mphi = static_cast<T>(m_)*phi;
        umnk.first = kernel*std::cos(mphi);
        umnk.second = kernel*std::sin(mphi);
    }
    return umnk;
}

template <typename T>
std::tuple<T, T, T, T, T, T> OblateMode<T>::velocity_mode(const int n, const T eta, const T tau, const T phi) const {
    const T sigma = sigmamnk_.at(n-1);
    const auto[X, Y] = get_X_and_Y(eta, tau, sigma);
    const auto[term1, term2, term3] = get_terms_123(X, Y);
    T sigma2 = sigma*sigma;
    T kernel_uradius, kernel_uang, kernel_uphi;
    std::tuple<T, T, T, T, T, T> umnk_vec;
    T ecc2 = eccentricity_*eccentricity_;
    T tau2 = tau*tau;
    T eta2 = eta*eta;
    T dv = std::sqrt(eta2+ecc2*tau2);
    T one_sub_sigma2 = static_cast<T>(1)-sigma2;
    T one_sub_ecc2sigma2 = static_cast<T>(1)-ecc2*sigma2;
    T one_sub_ecc2sigma2_sqrt = std::sqrt(one_sub_ecc2sigma2);
    T one_sub_ecc2sigma2_mut_one_sub_sigma2_sqrt = std::sqrt(one_sub_ecc2sigma2*one_sub_sigma2);
    T scale1 = eta*std::sqrt(static_cast<T>(1)-tau2) / dv;
    T scale2 = tau*std::sqrt(eta2+ecc2) / dv;
    if (std::isnan(scale1)) scale1 = static_cast<T>(0);
    if (std::isnan(scale2)) scale2 = static_cast<T>(0);
    if (m_ == 0) {
        T sigmaterm1scaled = sigma*term1/one_sub_ecc2sigma2_mut_one_sub_sigma2_sqrt;
        kernel_uradius = scale1*sigmaterm1scaled;
        kernel_uang = scale2*sigmaterm1scaled;
        kernel_uphi = term1/one_sub_ecc2sigma2_mut_one_sub_sigma2_sqrt;
    } else {
        T tmp = (sigma*term1-term2)/one_sub_ecc2sigma2_mut_one_sub_sigma2_sqrt;
        kernel_uradius = tmp*scale1;
        kernel_uang = tmp*scale2;
        kernel_uphi = (term1-sigma*term2)/one_sub_ecc2sigma2_mut_one_sub_sigma2_sqrt;
    }
    kernel_uradius -= term3*scale2/one_sub_ecc2sigma2_sqrt;
    kernel_uang += term3*scale1/one_sub_ecc2sigma2_sqrt;
    kernel_uradius /= static_cast<T>(2);
    kernel_uang /= static_cast<T>(-2);
    kernel_uphi /= static_cast<T>(-2);
    if (m_ == 0) {
        std::get<0>(umnk_vec) = static_cast<T>(0);
        std::get<1>(umnk_vec) = kernel_uradius;
        std::get<2>(umnk_vec) = static_cast<T>(0);
        std::get<3>(umnk_vec) = kernel_uang;
        std::get<4>(umnk_vec) = kernel_uphi;
        std::get<5>(umnk_vec) = static_cast<T>(0);
    } else {
        T mphi = static_cast<T>(m_)*phi;
        T sinmphi = std::sin(mphi);
        T cosmphi = std::cos(mphi);
        std::get<0>(umnk_vec) = -sinmphi*kernel_uradius;
        std::get<1>(umnk_vec) = cosmphi*kernel_uradius;
        std::get<2>(umnk_vec) = -sinmphi*kernel_uang;
        std::get<3>(umnk_vec) = cosmphi*kernel_uang;
        std::get<4>(umnk_vec) = cosmphi*kernel_uphi;
        std::get<5>(umnk_vec) = sinmphi*kernel_uphi;
    }
    return umnk_vec;
}

template <typename T>
std::pair<T, T> OblateMode<T>::divergence_velocity_mode(const int n, const T eta, const T tau, const T phi) const {
    const T sigma = sigmamnk_.at(n-1);
    T sigma2 = sigma*sigma;
    T one_sub_sigma2 = T(1)-sigma2;
    const auto[X, Y] = get_X_and_Y(eta, tau, sigma);
    auto [plmX, dplmX] = alegendre_.dALegendre_d1_compute(l_, m_, X);
    auto [plmY, dplmY] = alegendre_.dALegendre_d1_compute(l_, m_, Y);
    T X2 = X*X;
    T Y2 = Y*Y;
    T one_sub_X2 = T(1)-X2;
    T one_sub_Y2 = T(1)-Y2;
    dplmX /= std::sqrt(one_sub_X2);
    dplmY /= std::sqrt(one_sub_Y2);
    T m2 = T(m_*m_);
    T lscale = T(l_*(l_+1));
    T plmXplmY = plmX*plmY;

    // d(J gad X dot u_{mnk}) / dX
    T div_X = (m2/(one_sub_X2) - lscale)*plmXplmY;
    // d(J gad Y dot u_{mnk}) / dY
    T div_Y = (lscale - m2/(one_sub_Y2))*plmXplmY;
    // d(J gad phi dot u_{mnk}) / dphi
    T div_phi = T(0);

    if (m_ != 0) {
        T XdplmXplmY = X*dplmX*plmY;
        T YplmXdplmY = Y*plmX*dplmY;
        div_X -= T(m_)*(plmXplmY + XdplmXplmY)/sigma;
        div_Y += T(m_)*(plmXplmY + YplmXdplmY)/sigma;
        div_phi = (XdplmXplmY - YplmXdplmY - sigma*T(m_)*(X2-Y2)/(one_sub_X2*one_sub_Y2)*plmXplmY)*T(m_)/sigma;
    }

    T divergence = (div_X+div_Y+div_phi)/(T(2)*one_sub_sigma2)*std::sqrt(T(1)-eccentricity_*eccentricity_*sigma2);
    if (std::isnan(divergence) || std::isinf(divergence)) divergence = T(0);
    std::pair<T, T> divmode;
    if (m_ == 0) {
        divmode.first = T(0);
        divmode.second = divergence;
    } else {
        T mphi = T(m_)*phi;
        T cosmphi = std::cos(mphi);
        T sinmphi = std::sin(mphi);
        divmode.first = -sinmphi*divergence;
        divmode.second = cosmphi*divergence;
    }
    return divmode;
}

#endif 