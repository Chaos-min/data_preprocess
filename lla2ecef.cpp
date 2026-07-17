// cgcs2000_lla_ecef.cpp
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>

// CGCS2000 椭球参数（与WGS84基本一致，但扁率略有不同）
namespace CGCS2000 {
    constexpr double a = 6378137.0;                 // 长半轴 (m)
    constexpr double f = 1.0 / 298.257222101;       // 扁率
    constexpr double e2 = 2 * f - f * f;            // 第一偏心率平方
    constexpr double deg2rad = M_PI / 180.0;
    constexpr double rad2deg = 180.0 / M_PI;
}

// 将角度归一化到 [-180, 180)
double normalizeLongitude(double lon_deg) {
    lon_deg = fmod(lon_deg, 360.0);
    if (lon_deg >= 180.0) lon_deg -= 360.0;
    if (lon_deg < -180.0) lon_deg += 360.0;
    return lon_deg;
}

// LLA -> ECEF
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

// ECEF -> LLA (Bowring迭代)
void ecefToLla(double X, double Y, double Z,
               double &lat_deg, double &lon_deg, double &alt_m) {
    using namespace CGCS2000;
    double p = sqrt(X*X + Y*Y);
    if (p < 1e-12) { // 在极轴附近
        lon_deg = 0.0;
        double signZ = (Z >= 0) ? 1.0 : -1.0;
        lat_deg = signZ * 90.0;
        alt_m = fabs(Z) - a * (1.0 - e2);
        return;
    }

    double lon = atan2(Y, X);
    lon_deg = lon * rad2deg;
    lon_deg = normalizeLongitude(lon_deg);

    // 初始纬度估计
    double lat = atan2(Z, p * (1.0 - e2));
    double prevLat;
    int iter = 0;
    const int MAX_ITER = 30;
    const double CONV = 1e-12; // 弧度收敛阈值

    do {
        prevLat = lat;
        double sinLat = sin(lat), cosLat = cos(lat);
        double N = a / sqrt(1.0 - e2 * sinLat * sinLat);
        double h = p / cosLat - N; // 临时高度，用于更新纬度
        lat = atan2(Z, p * (1.0 - e2 * N / (N + h)));
        iter++;
    } while (fabs(lat - prevLat) > CONV && iter < MAX_ITER);

    if (iter >= MAX_ITER) {
        throw std::runtime_error("ECEF->LLA: 迭代未收敛");
    }

    // 最终计算高度
    double sinLat = sin(lat), cosLat = cos(lat);
    double N = a / sqrt(1.0 - e2 * sinLat * sinLat);
    // 若纬度接近 ±90°，用 Z 计算高度更稳定
    if (fabs(cosLat) < 1e-12) {
        alt_m = Z / sinLat - N * (1.0 - e2);
    } else {
        alt_m = p / cosLat - N;
    }

    lat_deg = lat * rad2deg;
    // 经度已经归一化，高度单位米
}

// ---------- 测试示例 ----------
int main() {
    // 测试点：北京天安门（约39.9042°N, 116.4074°E, 海拔50m）
    double lat0 = 39.9042, lon0 = 116.4074, alt0 = 50.0;
    double X, Y, Z;
    llaToEcef(lat0, lon0, alt0, X, Y, Z);
    std::cout << std::setprecision(12) << std::fixed;
    std::cout << "LLA -> ECEF:\n";
    std::cout << "X = " << X << " m\n";
    std::cout << "Y = " << Y << " m\n";
    std::cout << "Z = " << Z << " m\n";

    double lat1, lon1, alt1;
    ecefToLla(X, Y, Z, lat1, lon1, alt1);
    std::cout << "\nECEF -> LLA (恢复):\n";
    std::cout << "Lat = " << lat1 << " deg\n";
    std::cout << "Lon = " << lon1 << " deg\n";
    std::cout << "Alt = " << alt1 << " m\n";

    // 检查误差
    std::cout << "\n误差:\n";
    std::cout << "Δlat = " << (lat1 - lat0) * 1e6 << " μdeg\n";
    std::cout << "Δlon = " << (lon1 - lon0) * 1e6 << " μdeg\n";
    std::cout << "Δalt = " << (alt1 - alt0) * 1000 << " mm\n";

    return 0;
}