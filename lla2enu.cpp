#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>

// ==================== 椭球参数 ====================
namespace CGCS2000 {
    constexpr double a = 6378137.0;                 // 长半轴 (m)
    constexpr double f = 1.0 / 298.257222101;       // 扁率
    constexpr double e2 = 2 * f - f * f;            // 第一偏心率平方
    constexpr double deg2rad = M_PI / 180.0;
    constexpr double rad2deg = 180.0 / M_PI;
}

// ==================== 辅助函数 ====================
double normalizeLongitude(double lon_deg) {
    lon_deg = std::fmod(lon_deg, 360.0);
    if (lon_deg >= 180.0) lon_deg -= 360.0;
    if (lon_deg < -180.0) lon_deg += 360.0;
    return lon_deg;
}

// ==================== LLA <-> ECEF ====================
void llaToEcef(double lat_deg, double lon_deg, double alt_m,
               double &X, double &Y, double &Z) {
    using namespace CGCS2000;
    double lat = lat_deg * deg2rad;
    double lon = normalizeLongitude(lon_deg) * deg2rad;
    double sinLat = std::sin(lat), cosLat = std::cos(lat);
    double sinLon = std::sin(lon), cosLon = std::cos(lon);

    double N = a / std::sqrt(1.0 - e2 * sinLat * sinLat);
    double nh = N + alt_m;
    X = nh * cosLat * cosLon;
    Y = nh * cosLat * sinLon;
    Z = (N * (1.0 - e2) + alt_m) * sinLat;
}

void ecefToLla(double X, double Y, double Z,
               double &lat_deg, double &lon_deg, double &alt_m) {
    using namespace CGCS2000;
    double p = std::sqrt(X*X + Y*Y);
    if (p < 1e-12) { // 极轴附近
        lon_deg = 0.0;
        lat_deg = (Z >= 0) ? 90.0 : -90.0;
        alt_m = std::fabs(Z) - a * (1.0 - e2);
        return;
    }

    double lon = std::atan2(Y, X);
    lon_deg = normalizeLongitude(lon * rad2deg);

    double lat = std::atan2(Z, p * (1.0 - e2));
    double prevLat;
    int iter = 0;
    const int MAX_ITER = 30;
    const double CONV = 1e-12;
    do {
        prevLat = lat;
        double sinLat = std::sin(lat), cosLat = std::cos(lat);
        double N = a / std::sqrt(1.0 - e2 * sinLat * sinLat);
        double h = p / cosLat - N;
        lat = std::atan2(Z, p * (1.0 - e2 * N / (N + h)));
        iter++;
    } while (std::fabs(lat - prevLat) > CONV && iter < MAX_ITER);

    if (iter >= MAX_ITER)
        throw std::runtime_error("ECEF->LLA 迭代未收敛");

    double sinLat = std::sin(lat), cosLat = std::cos(lat);
    double N = a / std::sqrt(1.0 - e2 * sinLat * sinLat);
    if (std::fabs(cosLat) < 1e-12)
        alt_m = Z / sinLat - N * (1.0 - e2);
    else
        alt_m = p / cosLat - N;

    lat_deg = lat * rad2deg;
}

// ==================== LLA <-> ENU 直接接口 ====================
/**
 * 将目标点的大地坐标 (lat, lon, alt) 转换为相对于参考点 (ref_lat, ref_lon, ref_alt) 的 ENU 坐标
 */
void llaToEnu(double lat_deg, double lon_deg, double alt_m,
              double ref_lat_deg, double ref_lon_deg, double ref_alt_m,
              double &E, double &N, double &U) {
    // 1. 两点转ECEF
    double X, Y, Z, X0, Y0, Z0;
    llaToEcef(lat_deg, lon_deg, alt_m, X, Y, Z);
    llaToEcef(ref_lat_deg, ref_lon_deg, ref_alt_m, X0, Y0, Z0);

    // 2. 相对向量
    double dx = X - X0;
    double dy = Y - Y0;
    double dz = Z - Z0;

    // 3. 旋转矩阵（基于参考点）
    using namespace CGCS2000;
    double lat = ref_lat_deg * deg2rad;
    double lon = ref_lon_deg * deg2rad;
    double sinLat = std::sin(lat), cosLat = std::cos(lat);
    double sinLon = std::sin(lon), cosLon = std::cos(lon);

    E = -sinLon * dx + cosLon * dy;
    N = -sinLat * cosLon * dx - sinLat * sinLon * dy + cosLat * dz;
    U =  cosLat * cosLon * dx + cosLat * sinLon * dy + sinLat * dz;
}

/**
 * 将 ENU 坐标转换为目标点的大地坐标 (lat, lon, alt)
 */
void enuToLla(double E, double N, double U,
              double ref_lat_deg, double ref_lon_deg, double ref_alt_m,
              double &lat_deg, double &lon_deg, double &alt_m) {
    // 1. 参考点转ECEF
    double X0, Y0, Z0;
    llaToEcef(ref_lat_deg, ref_lon_deg, ref_alt_m, X0, Y0, Z0);

    // 2. 逆旋转（R^T）
    using namespace CGCS2000;
    double lat = ref_lat_deg * deg2rad;
    double lon = ref_lon_deg * deg2rad;
    double sinLat = std::sin(lat), cosLat = std::cos(lat);
    double sinLon = std::sin(lon), cosLon = std::cos(lon);

    double dx = -sinLon * E - sinLat * cosLon * N + cosLat * cosLon * U;
    double dy =  cosLon * E - sinLat * sinLon * N + cosLat * sinLon * U;
    double dz =                        cosLat * N + sinLat * U;

    // 3. 目标点ECEF
    double X = X0 + dx;
    double Y = Y0 + dy;
    double Z = Z0 + dz;

    // 4. ECEF -> LLA
    ecefToLla(X, Y, Z, lat_deg, lon_deg, alt_m);
}

// ==================== 测试 ====================
int main() {
    // 参考点：北京天安门 (39.9042°N, 116.4074°E, 椭球高50m)
    double ref_lat = 39.9042, ref_lon = 116.4074, ref_alt = 50.0;

    // 目标点：参考点正东100m，正北200m，高10m（ENU坐标）
    double E_test = 100.0, N_test = 200.0, U_test = 10.0;

    // ENU -> LLA
    double lat, lon, alt;
    enuToLla(E_test, N_test, U_test, ref_lat, ref_lon, ref_alt, lat, lon, alt);

    // LLA -> ENU 恢复
    double E_rec, N_rec, U_rec;
    llaToEnu(lat, lon, alt, ref_lat, ref_lon, ref_alt, E_rec, N_rec, U_rec);

    std::cout << std::setprecision(12) << std::fixed;
    std::cout << "参考点: lat=" << ref_lat << "°, lon=" << ref_lon << "°, alt=" << ref_alt << " m\n";
    std::cout << "原始ENU: E=" << E_test << " m, N=" << N_test << " m, U=" << U_test << " m\n";
    std::cout << "\nENU -> LLA:\n";
    std::cout << "Lat = " << lat << " °\n";
    std::cout << "Lon = " << lon << " °\n";
    std::cout << "Alt = " << alt << " m\n";
    std::cout << "\nLLA -> ENU (恢复):\n";
    std::cout << "E = " << E_rec << " m\n";
    std::cout << "N = " << N_rec << " m\n";
    std::cout << "U = " << U_rec << " m\n";

    std::cout << "\n误差 (恢复后与原值比较):\n";
    std::cout << "ΔE = " << (E_rec - E_test) * 1000 << " mm\n";
    std::cout << "ΔN = " << (N_rec - N_test) * 1000 << " mm\n";
    std::cout << "ΔU = " << (U_rec - U_test) * 1000 << " mm\n";

    return 0;
}