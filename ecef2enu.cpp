// cgcs2000_enu_convert.cpp
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>

// ---------- CGCS2000 椭球参数 ----------
namespace CGCS2000 {
    constexpr double a = 6378137.0;                 // 长半轴 (m)
    constexpr double f = 1.0 / 298.257222101;       // 扁率
    constexpr double e2 = 2 * f - f * f;            // 第一偏心率平方
    constexpr double deg2rad = M_PI / 180.0;
    constexpr double rad2deg = 180.0 / M_PI;
}

// ---------- 辅助函数 ----------
double normalizeLongitude(double lon_deg) {
    lon_deg = fmod(lon_deg, 360.0);
    if (lon_deg >= 180.0) lon_deg -= 360.0;
    if (lon_deg < -180.0) lon_deg += 360.0;
    return lon_deg;
}

// LLA -> ECEF (复用)
void llaToEcef(double lat_deg, double lon_deg, double alt_m,
               double &X, double &Y, double &Z) {
    using namespace CGCS2000;
    double lat = lat_deg * deg2rad;
    double lon = normalizeLongitude(lon_deg) * deg2rad;
    double sinLat = sin(lat), cosLat = cos(lat);
    double sinLon = sin(lon), cosLon = cos(lon);

    double N = a / sqrt(1.0 - e2 * sinLat * sinLat);
    double nh = N + alt_m;
    X = nh * cosLat * cosLon;
    Y = nh * cosLat * sinLon;
    Z = (N * (1.0 - e2) + alt_m) * sinLat;
}

// ECEF -> LLA (复用，用于验证，但ENU转换本身不需要)
void ecefToLla(double X, double Y, double Z,
               double &lat_deg, double &lon_deg, double &alt_m) {
    using namespace CGCS2000;
    double p = sqrt(X*X + Y*Y);
    if (p < 1e-12) {
        lon_deg = 0.0;
        lat_deg = (Z >= 0) ? 90.0 : -90.0;
        alt_m = fabs(Z) - a * (1.0 - e2);
        return;
    }
    double lon = atan2(Y, X);
    lon_deg = normalizeLongitude(lon * rad2deg);
    double lat = atan2(Z, p * (1.0 - e2));
    double prevLat;
    int iter = 0;
    const int MAX_ITER = 30;
    const double CONV = 1e-12;
    do {
        prevLat = lat;
        double sinLat = sin(lat), cosLat = cos(lat);
        double N = a / sqrt(1.0 - e2 * sinLat * sinLat);
        double h = p / cosLat - N;
        lat = atan2(Z, p * (1.0 - e2 * N / (N + h)));
        iter++;
    } while (fabs(lat - prevLat) > CONV && iter < MAX_ITER);
    if (iter >= MAX_ITER) throw std::runtime_error("ECEF->LLA 迭代未收敛");

    double sinLat = sin(lat), cosLat = cos(lat);
    double N = a / sqrt(1.0 - e2 * sinLat * sinLat);
    if (fabs(cosLat) < 1e-12)
        alt_m = Z / sinLat - N * (1.0 - e2);
    else
        alt_m = p / cosLat - N;
    lat_deg = lat * rad2deg;
}

// ---------- ECEF <-> ENU 核心转换 ----------
/**
 * ECEF -> ENU
 * @param X,Y,Z       待转换点的ECEF坐标 (m)
 * @param ref_lat_deg 参考点纬度 (度)
 * @param ref_lon_deg 参考点经度 (度)
 * @param ref_alt_m   参考点椭球高 (m)
 * @param E,N,U       输出的ENU坐标 (m)
 */
void ecefToEnu(double X, double Y, double Z,
               double ref_lat_deg, double ref_lon_deg, double ref_alt_m,
               double &E, double &N, double &U) {
    // 1. 参考点ECEF
    double X0, Y0, Z0;
    llaToEcef(ref_lat_deg, ref_lon_deg, ref_alt_m, X0, Y0, Z0);

    // 2. 相对向量
    double dx = X - X0;
    double dy = Y - Y0;
    double dz = Z - Z0;

    // 3. 旋转矩阵元素
    using namespace CGCS2000;
    double lat = ref_lat_deg * deg2rad;
    double lon = ref_lon_deg * deg2rad;  // 不需要归一化，但建议归一化
    double sinLat = sin(lat), cosLat = cos(lat);
    double sinLon = sin(lon), cosLon = cos(lon);

    // R 矩阵
    // E = -sinLon*dx + cosLon*dy
    // N = -sinLat*cosLon*dx - sinLat*sinLon*dy + cosLat*dz
    // U =  cosLat*cosLon*dx + cosLat*sinLon*dy + sinLat*dz
    E = -sinLon * dx + cosLon * dy;
    N = -sinLat * cosLon * dx - sinLat * sinLon * dy + cosLat * dz;
    U =  cosLat * cosLon * dx + cosLat * sinLon * dy + sinLat * dz;
}

/**
 * ENU -> ECEF
 * @param E,N,U        ENU坐标 (m)
 * @param ref_lat_deg  参考点纬度 (度)
 * @param ref_lon_deg  参考点经度 (度)
 * @param ref_alt_m    参考点椭球高 (m)
 * @param X,Y,Z        输出的ECEF坐标 (m)
 */
void enuToEcef(double E, double N, double U,
               double ref_lat_deg, double ref_lon_deg, double ref_alt_m,
               double &X, double &Y, double &Z) {
    // 1. 参考点ECEF
    double X0, Y0, Z0;
    llaToEcef(ref_lat_deg, ref_lon_deg, ref_alt_m, X0, Y0, Z0);

    // 2. 旋转矩阵 (与ecefToEnu相同)
    using namespace CGCS2000;
    double lat = ref_lat_deg * deg2rad;
    double lon = ref_lon_deg * deg2rad;
    double sinLat = sin(lat), cosLat = cos(lat);
    double sinLon = sin(lon), cosLon = cos(lon);

    // 逆旋转: dx = R^T * [E;N;U]
    double dx = -sinLon * E - sinLat * cosLon * N + cosLat * cosLon * U;
    double dy =  cosLon * E - sinLat * sinLon * N + cosLat * sinLon * U;
    double dz =                        cosLat * N + sinLat * U;

    X = X0 + dx;
    Y = Y0 + dy;
    Z = Z0 + dz;
}

// ---------- 测试示例 ----------
int main() {
    // 参考点：北京天安门 (39.9042°N, 116.4074°E, 椭球高50m)
    double ref_lat = 39.9042, ref_lon = 116.4074, ref_alt = 50.0;

    // 测试点：在参考点正东方向约100m，北方向约200m，高10m (ENU)
    double E_test = 100.0, N_test = 200.0, U_test = 10.0;

    // ENU -> ECEF
    double X, Y, Z;
    enuToEcef(E_test, N_test, U_test, ref_lat, ref_lon, ref_alt, X, Y, Z);

    // ECEF -> ENU 恢复
    double E_rec, N_rec, U_rec;
    ecefToEnu(X, Y, Z, ref_lat, ref_lon, ref_alt, E_rec, N_rec, U_rec);

    std::cout << std::setprecision(12) << std::fixed;
    std::cout << "参考点: lat=" << ref_lat << "°, lon=" << ref_lon << "°, alt=" << ref_alt << "m\n";
    std::cout << "原始ENU: E=" << E_test << " m, N=" << N_test << " m, U=" << U_test << " m\n";
    std::cout << "\nENU -> ECEF:\n";
    std::cout << "X = " << X << " m\n";
    std::cout << "Y = " << Y << " m\n";
    std::cout << "Z = " << Z << " m\n";
    std::cout << "\nECEF -> ENU (恢复):\n";
    std::cout << "E = " << E_rec << " m\n";
    std::cout << "N = " << N_rec << " m\n";
    std::cout << "U = " << U_rec << " m\n";

    std::cout << "\n误差:\n";
    std::cout << "ΔE = " << (E_rec - E_test) * 1000 << " mm\n";
    std::cout << "ΔN = " << (N_rec - N_test) * 1000 << " mm\n";
    std::cout << "ΔU = " << (U_rec - U_test) * 1000 << " mm\n";

    return 0;
}