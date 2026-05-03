#ifndef Legendre_HPP
#define Legendre_HPP

#include <cmath>
#include <tuple>
#include <iomanip>

extern "C" {
    void __besseleval_MOD_bessel_eval_init(double* dsize1);
    void __alegendreeval_MOD_alegendre_eval_init(double* disze2);
    void __alegendreeval_MOD_alegendre_eval(double* dnu, double* dmu, double* t, double* vallogp, double* valp);
}

template <typename T>
class ALegendre {
    public:
        ALegendre() { 
            double disze1, disze2; 
            __besseleval_MOD_bessel_eval_init(&disze1); 
            __alegendreeval_MOD_alegendre_eval_init(&disze2); }
        
        ~ALegendre() = default;

        T ALegendre_compute(const int l, const int m, T x);
        std::tuple<T, T> dALegendre_d1_compute(const int l, const int m, T x);
        std::tuple<T, T, T> ALegendre_compositions(const int l, const int m, T x);
        
};

template <typename T>
T ALegendre<T>::ALegendre_compute(const int l, const int m, T x) {
/*  x = cos(theta).
*/
    double x_abs = std::abs(x);
    if (x_abs > T(1)) x = (x < 0) ? T(-1) : T(1);
    x_abs = std::abs(x);

    T plm;

    if (x_abs == 1.0) {
        if (m != 0) plm = T(0);
        else {
            plm = 1.0;
            if (l%2 == 1 && x < T(0)) plm = -1.0;
        }
        return plm;
    }

    double dnu = double(l);
    int m_abs = std::abs(m);
    double dmu = double(m_abs);
    double t;
    bool sign_flag = false;
    
    if (x < T(0)) {
        t = acosl(x_abs);
        if ((l+m_abs)%2 == 1) sign_flag = true;
    } else {
        t = acosl(x);
    }

    double sintheta_sqrt = sqrtl(sinl(t));

    double vallogp, valp;
    __alegendreeval_MOD_alegendre_eval(&dnu, &dmu, &t, &vallogp, &valp);

    // valp *= std::sqrtl(2.0L);
    valp /= sintheta_sqrt;

    if (sign_flag) valp = -valp;
    if (m%2 == -1) valp = -valp;

    T scale = std::sqrt(T(l+l+1)/T(2));

    return T(valp)/scale;
}

template <typename T>
std::tuple<T, T> ALegendre<T>::dALegendre_d1_compute(const int l, const int m, T x) {
    /*  x = cos(theta).
*/  
    double x_abs = std::abs(x);
    if (x_abs > T(1)) x = (x < 0) ? T(-1) : T(1);
    x_abs = std::abs(x);

    int m_abs = (m < 0) ? -m : m;

    double plm, dplm;
    if (x_abs == T(1)) {
        if (m != 0) {
            plm = 0.0;
            if (m_abs == 1) {
                dplm = -std::sqrt(double(l*(l+1)))/2.0;
                if (x < T(0) && (l+m_abs)%2 == 0) {
                    dplm = -dplm;
                }
                if (m%2 == -1) dplm = -dplm;
            } else {
                dplm = 0.0;
            }
        } else {
            plm = 1.0;
            if (l%2 == 1 && x < T(0)) plm = -1.0;
            dplm = 0.0;
        }
    } else {
        double dnu = double(l);
        double dmu = double(m_abs);
        double t;

        if (x < T(0)) t = acosl(-x);
        else t = acosl(x);

        double sintheta_sqrt = sqrtl(sinl(t));

        double vallogp;
        __alegendreeval_MOD_alegendre_eval(&dnu, &dmu, &t, &vallogp, &plm);
        plm /= sintheta_sqrt;

        if (m_abs == 0) {
            double mu = 1.0;
            __alegendreeval_MOD_alegendre_eval(&dnu, &mu, &t, &vallogp, &dplm);
            dplm *= sqrtl(l*(l+1));
        } else {
            double dmu_sub_one = double(m_abs-1);
            double dmu_add_one = double(m_abs+1);
            double valp_sub_one, valp_add_one;
            __alegendreeval_MOD_alegendre_eval(&dnu, &dmu_sub_one, &t, &vallogp, &valp_sub_one);
            __alegendreeval_MOD_alegendre_eval(&dnu, &dmu_add_one, &t, &vallogp, &valp_add_one);
            dplm = (std::sqrt(double((l-m_abs)*(l+m_abs+1)))*valp_add_one - std::sqrt(double((l+m_abs)*(l-m_abs+1)))*valp_sub_one)/T(2);
        }
        dplm /= sintheta_sqrt;

        double scale = std::sqrt(double(l+l+1)/2.0);
        plm /= scale;
        dplm /= scale;

        if (x < T(0)) {
            if ((l+m_abs)%2 == 1) plm = -plm;
            else dplm = -dplm;
        }

        if (m%2 == -1) {
            plm = -plm;
            dplm = -dplm;
        }
    }

    // return {P_l^m, sqrt{1-x^2} dP_l^m/dx}
    return {T(plm), T(dplm)};
}

template <typename T>
std::tuple<T, T, T> ALegendre<T>::ALegendre_compositions(const int l, const int m, T x) {
    double x_abs = std::abs(x);
    if (x_abs > T(1)) x = (x < 0) ? -1.0 : 1.0;
    x_abs = std::abs(x);

    int m_abs = (m < 0) ? -m : m;

    double plm, dplm, plm_div;
    if (x_abs == 1.0) {
        if (m != 0) {
            plm = 0.0;

            if (m_abs == 1) {
                plm_div = std::sqrt(double(l*(l+1)))/2.0;
                dplm = -plm_div;
                if (x < T(0)) {
                    if ((l+m_abs)%2 == 1) plm_div = -plm_div;
                    else dplm = -dplm;
                }
                if (m%2 == -1) { plm_div = -plm_div; dplm = -dplm; }
            } else {
                plm_div = 0.0;
                dplm = 0.0;
            }
        } else {
            plm = 1.0;
            if (l%2 == 1 && x < T(0)) plm = -1.0;
            dplm = 0.0;
            plm_div = 0.0;
        }
    } else {
        double dnu = double(l);
        double dmu = double(m_abs);
        double t;

        if (x < T(0)) t = acosl(-x);
        else t = acosl(x);

        double sintheta_sqrt = sqrtl(sinl(t));
        double scale1 = sintheta_sqrt*std::sqrt(double(l+l+1)/2.0);

        double vallogp;
        __alegendreeval_MOD_alegendre_eval(&dnu, &dmu, &t, &vallogp, &plm);
        plm /= scale1;

        if (m_abs == 0) {
            double valp, mu = 1.0;
            __alegendreeval_MOD_alegendre_eval(&dnu, &mu, &t, &vallogp, &valp);
            dplm = valp/scale1*sqrtl(l*(l+1));
            plm_div = 0.0;
        } else {
            double ll = double(l+1);
            double dmu_sub_one = double(m_abs-1);
            double dmu_add_one = double(m_abs+1);
            double valp_sub_one, valp_add_one;
            __alegendreeval_MOD_alegendre_eval(&dnu, &dmu_sub_one, &t, &vallogp, &valp_sub_one);
            __alegendreeval_MOD_alegendre_eval(&dnu, &dmu_add_one, &t, &vallogp, &valp_add_one);
            dplm = (std::sqrt(double((l-m_abs)*(l+m_abs+1)))*valp_add_one - std::sqrt(double((l+m_abs)*(l-m_abs+1)))*valp_sub_one)/T(2);
            dplm /= scale1;

            __alegendreeval_MOD_alegendre_eval(&ll, &dmu_sub_one, &t, &vallogp, &valp_sub_one);
            __alegendreeval_MOD_alegendre_eval(&ll, &dmu_add_one, &t, &vallogp, &valp_add_one);

        
            plm_div = (std::sqrt(double((l+m_abs+2)*(l+m_abs+1)))*valp_add_one + std::sqrt(double((l-m_abs+1)*(l-m_abs+2)))*valp_sub_one)/(m_abs+m_abs);
            // plm_div *= std::sqrt(double(l+l+1)/double(l+l+3))/double(m_abs+m_abs);
            plm_div /= sintheta_sqrt*std::sqrt(double(l+l+3)/2.0);
        }
        // dplm /= sintheta_sqrt;
        // plm_div /= sintheta_sqrt;

        if (x < T(0)) {
            if ((l+m_abs)%2 == 1) { plm = -plm; plm_div = -plm_div; }
            else dplm = -dplm;
        }

        if (m%2 == -1) {
            plm = -plm;
            dplm = -dplm;
            plm_div = -plm_div;
        }
    }

    return {T(plm), T(dplm), T(plm_div)};
}

#endif