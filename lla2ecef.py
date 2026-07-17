import math
import sys

class CGCS2000:
    """CGCS2000椭球参数（与WGS84非常接近）"""
    a = 6378137.0          # 长半轴 (m)
    f = 1.0 / 298.257222101
    e2 = 2*f - f*f         # 第一偏心率平方
    deg2rad = math.pi / 180.0
    rad2deg = 180.0 / math.pi

    @staticmethod
    def normalize_longitude(lon_deg):
        """将经度归一化到 [-180, 180) 度"""
        lon_deg = lon_deg % 360.0
        if lon_deg >= 180.0:
            lon_deg -= 360.0
        if lon_deg < -180.0:
            lon_deg += 360.0
        return lon_deg

def lla_to_ecef(lat_deg, lon_deg, alt_m):
    """
    将大地坐标 (纬度, 经度, 椭球高) 转换为 ECEF 直角坐标 (X, Y, Z)
    单位：度, 度, 米 -> 米
    """
    a = CGCS2000.a
    e2 = CGCS2000.e2
    lat = math.radians(lat_deg)
    lon = math.radians(CGCS2000.normalize_longitude(lon_deg))
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    sin_lon = math.sin(lon)
    cos_lon = math.cos(lon)

    N = a / math.sqrt(1.0 - e2 * sin_lat * sin_lat)
    nh = N + alt_m
    X = nh * cos_lat * cos_lon
    Y = nh * cos_lat * sin_lon
    Z = (N * (1.0 - e2) + alt_m) * sin_lat
    return X, Y, Z

def ecef_to_lla(X, Y, Z):
    """
    将 ECEF 直角坐标 (X, Y, Z) 转换为大地坐标 (纬度, 经度, 椭球高)
    使用 Bowring 迭代法，收敛快、稳定
    返回：(lat_deg, lon_deg, alt_m)
    """
    a = CGCS2000.a
    e2 = CGCS2000.e2
    p = math.hypot(X, Y)

    # 极轴附近处理
    if p < 1e-12:
        lon_deg = 0.0
        lat_deg = 90.0 if Z >= 0 else -90.0
        alt_m = abs(Z) - a * (1.0 - e2)
        return lat_deg, lon_deg, alt_m

    lon = math.atan2(Y, X)
    lon_deg = math.degrees(lon)
    lon_deg = CGCS2000.normalize_longitude(lon_deg)

    # 初始纬度估计
    lat = math.atan2(Z, p * (1.0 - e2))
    max_iter = 30
    conv = 1e-12  # 弧度
    for _ in range(max_iter):
        prev_lat = lat
        sin_lat = math.sin(lat)
        cos_lat = math.cos(lat)
        N = a / math.sqrt(1.0 - e2 * sin_lat * sin_lat)
        # 用当前高度更新纬度
        h = p / cos_lat - N
        lat = math.atan2(Z, p * (1.0 - e2 * N / (N + h)))
        if abs(lat - prev_lat) < conv:
            break
    else:
        raise RuntimeError("ECEF->LLA 迭代未收敛")

    # 最终计算高度
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    N = a / math.sqrt(1.0 - e2 * sin_lat * sin_lat)
    if abs(cos_lat) < 1e-12:  # 极区
        alt_m = Z / sin_lat - N * (1.0 - e2)
    else:
        alt_m = p / cos_lat - N

    lat_deg = math.degrees(lat)
    return lat_deg, lon_deg, alt_m

# ---------- 测试示例 ----------
if __name__ == "__main__":
    # 测试点：北京天安门
    lat0, lon0, alt0 = 39.9042, 116.4074, 50.0
    X, Y, Z = lla_to_ecef(lat0, lon0, alt0)
    print("LLA -> ECEF:")
    print(f"X = {X:.6f} m")
    print(f"Y = {Y:.6f} m")
    print(f"Z = {Z:.6f} m")

    lat1, lon1, alt1 = ecef_to_lla(X, Y, Z)
    print("\nECEF -> LLA (恢复):")
    print(f"Lat = {lat1:.12f} deg")
    print(f"Lon = {lon1:.12f} deg")
    print(f"Alt = {alt1:.6f} m")

    print("\n误差:")
    print(f"Δlat = {(lat1 - lat0) * 1e6:.3f} μdeg")
    print(f"Δlon = {(lon1 - lon0) * 1e6:.3f} μdeg")
    print(f"Δalt = {(alt1 - alt0) * 1000:.3f} mm")