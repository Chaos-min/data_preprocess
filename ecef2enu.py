import math

# ---------- CGCS2000 椭球参数 ----------
class CGCS2000:
    a = 6378137.0
    f = 1.0 / 298.257222101
    e2 = 2*f - f*f
    deg2rad = math.pi / 180.0
    rad2deg = 180.0 / math.pi

    @staticmethod
    def normalize_longitude(lon_deg):
        lon_deg = lon_deg % 360.0
        if lon_deg >= 180.0:
            lon_deg -= 360.0
        if lon_deg < -180.0:
            lon_deg += 360.0
        return lon_deg

# ---------- LLA <-> ECEF (复用前文) ----------
def lla_to_ecef(lat_deg, lon_deg, alt_m):
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
    a = CGCS2000.a
    e2 = CGCS2000.e2
    p = math.hypot(X, Y)
    if p < 1e-12:
        lon_deg = 0.0
        lat_deg = 90.0 if Z >= 0 else -90.0
        alt_m = abs(Z) - a * (1.0 - e2)
        return lat_deg, lon_deg, alt_m
    lon = math.atan2(Y, X)
    lon_deg = math.degrees(lon)
    lon_deg = CGCS2000.normalize_longitude(lon_deg)
    lat = math.atan2(Z, p * (1.0 - e2))
    for _ in range(30):
        prev_lat = lat
        sin_lat = math.sin(lat)
        cos_lat = math.cos(lat)
        N = a / math.sqrt(1.0 - e2 * sin_lat * sin_lat)
        h = p / cos_lat - N
        lat = math.atan2(Z, p * (1.0 - e2 * N / (N + h)))
        if abs(lat - prev_lat) < 1e-12:
            break
    else:
        raise RuntimeError("ECEF->LLA 迭代未收敛")
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    N = a / math.sqrt(1.0 - e2 * sin_lat * sin_lat)
    if abs(cos_lat) < 1e-12:
        alt_m = Z / sin_lat - N * (1.0 - e2)
    else:
        alt_m = p / cos_lat - N
    lat_deg = math.degrees(lat)
    return lat_deg, lon_deg, alt_m

# ---------- ECEF <-> ENU ----------
def ecef_to_enu(X, Y, Z, ref_lat_deg, ref_lon_deg, ref_alt_m):
    """
    将ECEF坐标转换为以(ref_lat_deg, ref_lon_deg, ref_alt_m)为原点的ENU坐标
    返回 (E, N, U)
    """
    # 参考点ECEF
    X0, Y0, Z0 = lla_to_ecef(ref_lat_deg, ref_lon_deg, ref_alt_m)
    dx = X - X0
    dy = Y - Y0
    dz = Z - Z0

    lat = math.radians(ref_lat_deg)
    lon = math.radians(ref_lon_deg)
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    sin_lon = math.sin(lon)
    cos_lon = math.cos(lon)

    E = -sin_lon * dx + cos_lon * dy
    N = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz
    U =  cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz
    return E, N, U

def enu_to_ecef(E, N, U, ref_lat_deg, ref_lon_deg, ref_alt_m):
    """
    将ENU坐标转换回ECEF
    返回 (X, Y, Z)
    """
    X0, Y0, Z0 = lla_to_ecef(ref_lat_deg, ref_lon_deg, ref_alt_m)

    lat = math.radians(ref_lat_deg)
    lon = math.radians(ref_lon_deg)
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    sin_lon = math.sin(lon)
    cos_lon = math.cos(lon)

    # 逆旋转
    dx = -sin_lon * E - sin_lat * cos_lon * N + cos_lat * cos_lon * U
    dy =  cos_lon * E - sin_lat * sin_lon * N + cos_lat * sin_lon * U
    dz =                        cos_lat * N + sin_lat * U

    X = X0 + dx
    Y = Y0 + dy
    Z = Z0 + dz
    return X, Y, Z

# ---------- 测试 ----------
if __name__ == "__main__":
    # 参考点
    ref_lat, ref_lon, ref_alt = 39.9042, 116.4074, 50.0
    # 测试ENU
    E_test, N_test, U_test = 100.0, 200.0, 10.0

    X, Y, Z = enu_to_ecef(E_test, N_test, U_test, ref_lat, ref_lon, ref_alt)
    E_rec, N_rec, U_rec = ecef_to_enu(X, Y, Z, ref_lat, ref_lon, ref_alt)

    print(f"参考点: lat={ref_lat}°, lon={ref_lon}°, alt={ref_alt}m")
    print(f"原始ENU: E={E_test}m, N={N_test}m, U={U_test}m")
    print("\nENU -> ECEF:")
    print(f"X = {X:.6f} m")
    print(f"Y = {Y:.6f} m")
    print(f"Z = {Z:.6f} m")
    print("\nECEF -> ENU (恢复):")
    print(f"E = {E_rec:.6f} m")
    print(f"N = {N_rec:.6f} m")
    print(f"U = {U_rec:.6f} m")
    print("\n误差:")
    print(f"ΔE = {(E_rec - E_test)*1000:.3f} mm")
    print(f"ΔN = {(N_rec - N_test)*1000:.3f} mm")
    print(f"ΔU = {(U_rec - U_test)*1000:.3f} mm")