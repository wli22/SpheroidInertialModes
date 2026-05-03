#ifndef GEOSTROPHICMODE_HPP
#define GEOSTROPHICMODE_HPP

#include <iostream>
#include <cmath>
#include <utility>
#include "ALegendre.hpp"


template <typename T>
class GeostrophicMode {
    public:
        GeostrophicMode(int k, int delta, T eccentricity);  // spherical and oblate : delta = 0;  prolate : delta = 1

        GeostrophicMode(const GeostrophicMode<T>& other) : k_(other.k_), l_(other.l_), delta_(other.delta_), eccentricity_(other.eccentricity_), one_sub_ecc2_(other.one_sub_ecc2_), alpha_(other.alpha_), normalization_(other.normalization_) {}

        GeostrophicMode<T>& operator=(const GeostrophicMode<T>& other) {
            if (this != &other) {
                k_             = other.k_;
                l_             = other.l_;
                delta_         = other.delta_;
                eccentricity_  = other.eccentricity_;
                one_sub_ecc2_  = other.one_sub_ecc2_;
                alpha_         = other.alpha_;
                normalization_ = other.normalization_;
            }
        }

        ~GeostrophicMode() = default;

        T geostrophic_pressure_mode(T s) const;
        T geostrophic_velocity_mode(T s) const;
        T geostrophic_pressure_mode_normalized(T s) const { return geostrophic_pressure_mode(s)/normalization_; }
        T geostrophic_velocity_mode_normalized(T s) const { return geostrophic_velocity_mode(s)/normalization_; }

        T normalization() const { return normalization_; }
        T eccentricity() const { return eccentricity_; }

        void reset_order(int k) {
            if (k < 1) {
                std::cerr << "Wrong order\n";
                std::exit(EXIT_FAILURE);
            }
            k_ = k;
            l_ = k_+k_;
            normalization_computation();
        }

    private:
        void normalization_computation() { 
            T norm = T(l_*(l_+1))/T(l_+l_+1)*M_PIf64;
            if (delta_ == 0) norm *= std::sqrt(one_sub_ecc2_);
            normalization_ = std::sqrt(norm);
        }

        static inline ALegendre<double> alegendre_{};
        int k_;
        int l_;
        int delta_;
        T eccentricity_;
        T one_sub_ecc2_;
        T alpha_;
        T normalization_;
};


template <typename T>
GeostrophicMode<T>::GeostrophicMode(int k, int delta, T eccentricity) : k_(k), l_(k_+k_), delta_(delta), eccentricity_(eccentricity), one_sub_ecc2_(T(1)-eccentricity_*eccentricity_) {
    if (k < 1) {
        std::cerr << "Wrong order k\n"; 
        std::exit(EXIT_FAILURE);
    }

    if (delta_ == 0) alpha_ = T(1);
    else if (delta_ == 1) alpha_ = std::sqrt(one_sub_ecc2_);
    else {
        std::cerr << "Wrong delta\n";
        std::exit(EXIT_FAILURE);
    }

    normalization_computation();
}

template <typename T>
T GeostrophicMode<T>::geostrophic_pressure_mode(T s) const {
    double X = (delta_ == 0) ? static_cast<double>(T(1)-s*s) : static_cast<double>(T(1)-s*s/one_sub_ecc2_);
    if (X < 0.0) X = 0.0;
    X = std::sqrt(X);
    double p2k = alegendre_.ALegendre_compute(l_, 0, X);
    double j0 = alegendre_.ALegendre_compute(l_, 0, 1.0);
    return j0-p2k;
}

template <typename T>
T GeostrophicMode<T>::geostrophic_velocity_mode(T s) const {
    double X = (delta_ == 0) ? static_cast<double>(T(1)-s*s) : static_cast<double>(T(1)-s*s/one_sub_ecc2_);
    if (X < 0.0) X = 0.0;
    X = std::sqrt(X);
    T gmode = T(0);
    if (X != T(0)) {
        auto[plm, dplm] = alegendre_.dALegendre_d1_compute(l_, 0, X);
        gmode = (delta_ == 0) ? dplm/X/T(2) : dplm/X/T(2)/alpha_;
    } else {
        gmode = -T(l_*(l_+1))*alegendre_.ALegendre_compute(l_, 0, 0.0)/T(2)/alpha_;
    }
    return gmode;
}



#endif 