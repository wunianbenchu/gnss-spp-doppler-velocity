#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* ========== 常量定义 ========== */
#define MAX_LINE_LEN       4096
#define MAX_OBS_PER_SYS    128
#define MAX_SAT_PER_EPOCH  64
#define MAX_EPOCHS         4096
#define MAX_KEPLER_EPHEM   8192
#define MAX_GPS_SATS       32

#define GPS_F1    1575.42e6
#define GPS_F2    1227.60e6
#define C_LIGHT   299792458.0
#define GM        3.986004418e14
#define OMEGA_E   7.2921151467e-5
#define WGS84_A   6378137.0
#define WGS84_E2  0.0066943799901413165

#define SYS_GAL 0
#define SYS_BDS 1
#define SYS_GLO 2
#define SYS_GPS 3

/* GPS L1/L2 码优先级（实验一） */
static const char* GPS_L1_PRI = "CSLXPWYM";
static const char* GPS_L2_PRI = "CDSLXPWYM";

/* ========== 数据结构 ========== */

/* 单颗卫星双频观测（实验一、九） */
typedef struct {
    double C;      /* 伪距(m) */
    double L;      /* 载波相位(周) */
    double D;      /* 多普勒(Hz) */
    double S;      /* 信噪比(dBHz) */
    char   code;   /* 跟踪模式码 */
} FreqObs;

/* 单颗卫星完整观测 */
typedef struct {
    int    prn;
    FreqObs freq[2];  /* freq[0]=L1, freq[1]=L2 */
} SatObs;

/* 一个历元的观测数据 */
typedef struct {
    int    year, month, day, hour, minute;
    double second;
    SatObs sats[MAX_GPS_SATS];
    int    nsat;
} ObsEpoch;

/* 开普勒广播星历（实验二、五） */
typedef struct {
    char   sys;
    int    prn;
    int    year, month, day, hour, minute;
    double second;
    double a0, a1, a2;                         /* 3个钟差参数 */
    double sqrt_a, e, i0, Omega0, omega, M0;   /* 6个开普勒参数 */
    double i_dot, Omega_dot, delta_n;          /* 3个摄动参数 */
    double Cuc, Cus, Crc, Crs, Cic, Cis;       /* 6个谐波摄动参数 */
    double tgd;                                 /* TGD(s) */
    int    toe_week;                           /* TOE周数(模1024) */
    double toe;                                /* TOE周内秒 */
} KeplerEphem;

/* ========== 全局变量 ========== */
static KeplerEphem kepler_list[MAX_KEPLER_EPHEM];
static int         kepler_cnt = 0;
static ObsEpoch    epoch_list[MAX_EPOCHS];
static int         epoch_cnt = 0;
static double      true_x, true_y, true_z;    /* 真值坐标（从O文件头读取） */
static char        obs_types[4][MAX_OBS_PER_SYS][4] = { 0 };
static int         sys_obs_count[4] = { 0 };

/* ========== 工具函数 ========== */

/* D->E 替换（实验二） */
void replace_D_to_E(char* str) {
    while (*str) {
        if (*str == 'D' || *str == 'd') *str = 'E';
        str++;
    }
}

/* 处理 RINEX 负数前缺空格（实验二、五） */
void preprocess_rinex_line(char* line) {
    char buf[MAX_LINE_LEN * 2] = { 0 };
    int j = 0, len = (int)strlen(line);
    for (int i = 0; i < len; i++) {
        if (line[i] == '-' && i > 0 && line[i - 1] != ' ' &&
            line[i - 1] != 'E' && line[i - 1] != 'e')
            buf[j++] = ' ';
        buf[j++] = line[i];
    }
    buf[j] = '\0';
    strcpy(line, buf);
}

/* 系统索引（实验一） */
int get_sys_idx(char c) {
    switch (c) {
    case 'E': return SYS_GAL;
    case 'C': return SYS_BDS;
    case 'R': return SYS_GLO;
    case 'G': return SYS_GPS;
    default:  return -1;
    }
}

/* 提取 RINEX 头文件标签（实验一） */
void get_rinex_label(const char* line, char* label) {
    memset(label, 0, 21);
    size_t len = strlen(line);
    if (len >= 60) {
        strncpy(label, line + 60, 20);
        label[20] = '\0';
    }
    char* start = label;
    while (*start == ' ') start++;
    char* end = start + strlen(start) - 1;
    while (end > start && (*end == ' ' || *end == '\r' || *end == '\n')) end--;
    *(end + 1) = '\0';
    memmove(label, start, strlen(start) + 1);
}

/* 获取码优先级索引（实验一） */
int get_priority(char code, const char* pri_str) {
    char* pos = strchr(pri_str, code);
    return pos ? (int)(pos - pri_str) : -1;
}

/* 判断闰年（实验四） */
static int isLeapYear(int y) {
    return (y % 4 == 0 && y % 100 != 0) || (y % 400 == 0);
}

/* 某月天数（实验四） */
static int daysInMonth(int year, int month) {
    int d[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
    if (month == 2 && isLeapYear(year)) return 29;
    return d[month - 1];
}

/* 计算简化儒略日 MJD（实验四、五） */
double calc_mjd(int year, int month, int day,
                int hour, int minute, double second) {
    int y = year, m = month;
    if (m <= 2) { y--; m += 12; }
    long long a = y / 100, b = 2 - a + a / 4;
    double jd = floor(365.25 * (y + 4716)) +
                floor(30.6001 * (m + 1)) + day + b - 1524.5;
    return jd - 2400000.5 + hour / 24.0 +
           minute / 1440.0 + second / 86400.0;
}

/* GPS秒总数（从1980-01-06 00:00:00起，实验五） */
double calc_gps_seconds(int year, int month, int day,
                         int hour, int minute, double second) {
    double mjd_current = calc_mjd(year, month, day, hour, minute, second);
    double mjd_gps     = calc_mjd(1980, 1, 6, 0, 0, 0.0);
    return (mjd_current - mjd_gps) * 86400.0;
}

/* ========== 步骤1: 读取观测文件（实验一 + 实验九） ========== */
int read_obs_file(const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        printf("无法打开观测文件: %s\n", filename);
        return 0;
    }
    char line[MAX_LINE_LEN] = { 0 };
    char label[21] = { 0 };
    int  header_end = 0;

    /* 解析头文件 */
    while (fgets(line, MAX_LINE_LEN, fp)) {
        line[strcspn(line, "\r\n")] = '\0';
        get_rinex_label(line, label);

        if (strcmp(label, "END OF HEADER") == 0) {
            header_end = 1;
            break;
        }

        /* 观测类型（实验一） */
        if (strcmp(label, "SYS / # / OBS TYPES") == 0) {
            int cur_sys = -1;
            if (line[0] != ' ') {
                char sys_c = line[0];
                cur_sys = get_sys_idx(sys_c);
                if (cur_sys == -1) continue;
                sscanf(line + 4, "%d", &sys_obs_count[cur_sys]);
                memset(obs_types[cur_sys], 0, sizeof(obs_types[cur_sys]));

                int idx = 0;
                char* p = line + 6;
                while (idx < sys_obs_count[cur_sys] && idx < MAX_OBS_PER_SYS) {
                    while (*p == ' ') p++;
                    if (!*p || p - line >= 60) break;
                    strncpy(obs_types[cur_sys][idx], p, 3);
                    obs_types[cur_sys][idx][3] = '\0';
                    idx++;
                    p += 3;
                }
            }
            else {
                /* 续行：找已填充的 cur_sys */
                for (int s = 0; s < 4; s++) {
                    if (sys_obs_count[s] > 0) {
                        int cnt = 0;
                        while (cnt < sys_obs_count[s] && obs_types[s][cnt][0] != '\0') cnt++;
                        if (cnt < sys_obs_count[s]) { cur_sys = s; break; }
                    }
                }
                if (cur_sys == -1) continue;
                int idx = 0;
                while (idx < sys_obs_count[cur_sys] && obs_types[cur_sys][idx][0] != '\0') idx++;
                char* p = line + 6;
                while (idx < sys_obs_count[cur_sys] && idx < MAX_OBS_PER_SYS) {
                    while (*p == ' ') p++;
                    if (!*p || p - line >= 60) break;
                    strncpy(obs_types[cur_sys][idx], p, 3);
                    obs_types[cur_sys][idx][3] = '\0';
                    idx++;
                    p += 3;
                }
            }
        }

        /* 测站近似坐标（作为真值） */
        if (strcmp(label, "APPROX POSITION XYZ") == 0) {
            double xv, yv, zv;
            char x_str[15] = { 0 }, y_str[15] = { 0 }, z_str[15] = { 0 };
            strncpy(x_str, line, 14);      x_str[14] = '\0';
            strncpy(y_str, line + 14, 14); y_str[14] = '\0';
            strncpy(z_str, line + 28, 14); z_str[14] = '\0';
            sscanf(x_str, "%lf", &xv);
            sscanf(y_str, "%lf", &yv);
            sscanf(z_str, "%lf", &zv);
            true_x = xv; true_y = yv; true_z = zv;
            printf("测站近似坐标(真值): X=%.4f Y=%.4f Z=%.4f\n",
                   true_x, true_y, true_z);
        }

        memset(line, 0, MAX_LINE_LEN);
    }

    if (!header_end) { fclose(fp); return 0; }

    /* 解析观测数据：只取 GPS，取整点（minute==0, second≈0），只取 02:00~06:00 */
    int year, month, day, hour = 0, minute = 0, epoch_flag, sat_num;
    double second = 0.0;
    static SatObs sat_list[MAX_SAT_PER_EPOCH] = { 0 };
    int sat_cnt = 0;
    size_t line_len;

    while (fgets(line, MAX_LINE_LEN, fp)) {
        line[strcspn(line, "\r\n")] = '\0';
        line_len = strlen(line);

        if (line[0] == '>') {
            /* 保存上一个历元 */
            if (sat_cnt > 0 &&
                hour >= 2 && hour <= 6 &&
                minute == 0 && fabs(second) < 1e-3) {
                ObsEpoch* ep = &epoch_list[epoch_cnt];
                ep->year = year; ep->month = month; ep->day = day;
                ep->hour = hour; ep->minute = minute; ep->second = second;
                ep->nsat = 0;
                for (int i = 0; i < sat_cnt && ep->nsat < MAX_GPS_SATS; i++) {
                    if (sat_list[i].prn > 0 &&
                        fabs(sat_list[i].freq[0].C) >= 1000.0 &&
                        fabs(sat_list[i].freq[1].C) >= 1000.0) {
                        ep->sats[ep->nsat++] = sat_list[i];
                    }
                }
                if (ep->nsat > 0) epoch_cnt++;
            }
            memset(sat_list, 0, sizeof(sat_list));
            sat_cnt = 0;

            sscanf(line + 1, "%d %d %d %d %d %lf %d %d",
                   &year, &month, &day, &hour, &minute, &second,
                   &epoch_flag, &sat_num);
            if (year < 100) year += 1900;

            /* 非观测历元跳过 */
            if (epoch_flag != 0) {
                for (int i = 0; i < sat_num; i++)
                    fgets(line, MAX_LINE_LEN, fp);
                continue;
            }
            /* 跳过不在时间范围内的历元 */
            if (hour < 2 || hour > 6) continue;
        }
        else {
            if (sat_cnt >= MAX_SAT_PER_EPOCH || line_len < 3) continue;

            char sys = line[0];
            int prn;
            if (sscanf(line + 1, "%2d", &prn) != 1) continue;

            /* 只处理 GPS */
            if (sys != 'G') continue;

            int sys_idx = get_sys_idx(sys);
            if (sys_idx == -1 || sys_obs_count[sys_idx] == 0) continue;

            int obs_cnt = sys_obs_count[sys_idx];
            SatObs* sat = &sat_list[sat_cnt];
            sat->prn = prn;

            /* 提取所有观测值（实验一：字段位置 = 3 + i*16，值占14列） */
            double obs_vals[MAX_OBS_PER_SYS] = { 0 };
            for (int i = 0; i < obs_cnt && i < MAX_OBS_PER_SYS; i++) {
                int val_start = 3 + i * 16;
                if (val_start + 14 > (int)line_len) continue;
                char val_buf[15] = { 0 };
                strncpy(val_buf, line + val_start, 14);
                val_buf[14] = '\0';
                if (sscanf(val_buf, "%lf", &obs_vals[i]) != 1)
                    obs_vals[i] = 0.0;
            }

            /* GPS双频码优先级选择（实验一） */
            const char* pri_str[2] = { GPS_L1_PRI, GPS_L2_PRI };
            int  best_c_idx[2] = { -1, -1 };
            int  best_pri[2]   = { -1, -1 };
            char best_code[2]  = { 0, 0 };

            for (int f = 0; f < 2; f++) {
                char target_band = (f == 0) ? '1' : '2';
                for (int i = 0; i < obs_cnt && i < MAX_OBS_PER_SYS; i++) {
                    char* type = obs_types[sys_idx][i];
                    if (type[0] != 'C' || type[1] != target_band ||
                        fabs(obs_vals[i]) < 1000.0)
                        continue;
                    int p = get_priority(type[2], pri_str[f]);
                    if (p > best_pri[f]) {
                        best_pri[f]  = p;
                        best_code[f] = type[2];
                        best_c_idx[f] = i;
                    }
                }
            }

            /* 填充双频观测 */
            for (int f = 0; f < 2; f++) {
                if (best_c_idx[f] == -1) continue;
                sat->freq[f].code = best_code[f];
                sat->freq[f].C    = obs_vals[best_c_idx[f]];
                /* 尝试匹配同名载波、多普勒、信噪比 */
                char band = (f == 0) ? '1' : '2';
                char code = best_code[f];
                for (int i = 0; i < obs_cnt && i < MAX_OBS_PER_SYS; i++) {
                    char* type = obs_types[sys_idx][i];
                    if (type[1] != band || type[2] != code) continue;
                    if (type[0] == 'L') sat->freq[f].L = obs_vals[i];
                    if (type[0] == 'D') sat->freq[f].D = obs_vals[i];
                    if (type[0] == 'S') sat->freq[f].S = obs_vals[i];
                }
            }
            sat_cnt++;
        }
        memset(line, 0, MAX_LINE_LEN);
    }

    /* 最后一个历元 */
    if (sat_cnt > 0 &&
        hour >= 2 && hour <= 6 &&
        minute == 0 && fabs(second) < 1e-3) {
        ObsEpoch* ep = &epoch_list[epoch_cnt];
        ep->year = year; ep->month = month; ep->day = day;
        ep->hour = hour; ep->minute = minute; ep->second = second;
        ep->nsat = 0;
        for (int i = 0; i < sat_cnt && ep->nsat < MAX_GPS_SATS; i++) {
            if (sat_list[i].prn > 0 &&
                fabs(sat_list[i].freq[0].C) >= 1000.0 &&
                fabs(sat_list[i].freq[1].C) >= 1000.0) {
                ep->sats[ep->nsat++] = sat_list[i];
            }
        }
        if (ep->nsat > 0) epoch_cnt++;
    }

    fclose(fp);
    printf("观测文件解析完成: %d 个历元\n", epoch_cnt);
    return 1;
}

/* ========== 步骤2: 读取广播星历文件（实验二） ========== */
int read_nav_file(const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        printf("无法打开导航文件: %s\n", filename);
        return 0;
    }
    char line[MAX_LINE_LEN];
    int  header_end = 0;

    /* 跳过文件头 */
    while (fgets(line, MAX_LINE_LEN, fp)) {
        line[strcspn(line, "\r\n")] = '\0';
        if (strstr(line, "END OF HEADER")) {
            header_end = 1;
            break;
        }
    }

    if (!header_end) { fclose(fp); return 0; }

    /* 读取星历数据（实验二） */
    while (fgets(line, MAX_LINE_LEN, fp)) {
        line[strcspn(line, "\r\n")] = '\0';
        if (strlen(line) < 3) continue;

        char sys = line[0];
        int  prn;
        if (sscanf(line + 1, "%2d", &prn) != 1) continue;

        /* 只处理 GPS */
        if (sys != 'G') {
            for (int i = 0; i < 7; i++)
                fgets(line, MAX_LINE_LEN, fp);
            continue;
        }

        if (kepler_cnt >= MAX_KEPLER_EPHEM) break;

        KeplerEphem tmp;
        memset(&tmp, 0, sizeof(tmp));
        tmp.sys = sys;
        tmp.prn = prn;
        preprocess_rinex_line(line);

        if (sscanf(line, "%*c%2d %d %d %d %d %d %lf %lf %lf %lf",
            &tmp.prn, &tmp.year, &tmp.month, &tmp.day,
            &tmp.hour, &tmp.minute, &tmp.second,
            &tmp.a0, &tmp.a1, &tmp.a2) != 10) {
            for (int i = 0; i < 7; i++)
                fgets(line, MAX_LINE_LEN, fp);
            continue;
        }
        if (tmp.year < 100) tmp.year += 2000;

        /* 读取后续7行28个参数（实验二） */
        double buf[28] = { 0 };
        for (int i = 0; i < 7; i++) {
            if (!fgets(line, MAX_LINE_LEN, fp)) break;
            line[strcspn(line, "\r\n")] = '\0';
            replace_D_to_E(line);
            sscanf(line, "%lf %lf %lf %lf",
                   &buf[i * 4], &buf[i * 4 + 1],
                   &buf[i * 4 + 2], &buf[i * 4 + 3]);
        }

        /* 参数映射（实验二：RINEX 3.14 标准映射） */
        /* buf[0] = IODE */
        tmp.Crs       = buf[1];
        tmp.delta_n   = buf[2];
        tmp.M0        = buf[3];
        tmp.Cuc       = buf[4];
        tmp.e         = buf[5];
        tmp.Cus       = buf[6];
        tmp.sqrt_a    = buf[7];
        tmp.toe       = buf[8];
        tmp.Cic       = buf[9];
        tmp.Omega0    = buf[10];
        tmp.Cis       = buf[11];
        tmp.i0        = buf[12];
        tmp.Crc       = buf[13];
        tmp.omega     = buf[14];
        tmp.Omega_dot = buf[15];
        tmp.i_dot     = buf[16];
        tmp.toe_week  = (int)buf[18];
        /* buf[20] = SV accuracy, buf[21] = SV health */
        tmp.tgd       = buf[22];  /* GPS TGD */

        kepler_list[kepler_cnt++] = tmp;
    }

    fclose(fp);
    printf("广播星历解析完成: %d 条GPS星历\n", kepler_cnt);
    return 1;
}

/* ========== 步骤3: 卫星位置与钟差计算（实验五） ========== */
/**
 * 计算GPS卫星在时刻 t（GPS总秒数，从1980-01-06起）
 * 的ECEF位置、钟差（含相对论效应改正）
 * 严格模仿实验五的 compute_satellite_position
 */
void compute_sat_pos_clock(KeplerEphem* eph, double t,
                           double* x, double* y, double* z,
                           double* clk_bias, double* clk_rel) {
    /* 轨道长半轴 */
    double a = eph->sqrt_a * eph->sqrt_a;

    /* 平均角速度（含摄动改正） */
    double n0 = sqrt(GM / (a * a * a));
    double n  = n0 + eph->delta_n;

    /* GPS 周数翻转解析 */
    int full_week = (int)(t / 604800.0);
    int toe_week_resolved = eph->toe_week;
    if (toe_week_resolved < full_week - 512)
        toe_week_resolved += 1024;
    else if (toe_week_resolved > full_week + 512)
        toe_week_resolved -= 1024;
    double toe_full = toe_week_resolved * 604800.0 + eph->toe;
    double tk = t - toe_full;

    /* 偏近点角迭代（牛顿法，实验五，10次） */
    double M = eph->M0 + n * tk;
    double E = M;
    for (int i = 0; i < 10; i++) {
        E = M + eph->e * sin(E);
    }

    double sin_E = sin(E);
    double cos_E = cos(E);

    /* 真近点角 */
    double v = atan2(sqrt(1.0 - eph->e * eph->e) * sin_E,
                     cos_E - eph->e);

    /* 升交角距 */
    double phi = v + eph->omega;

    /* 摄动修正项（实验五） */
    double sin_2phi = sin(2.0 * phi);
    double cos_2phi = cos(2.0 * phi);
    double du = eph->Cuc * cos_2phi + eph->Cus * sin_2phi;
    double dr = eph->Crc * cos_2phi + eph->Crs * sin_2phi;
    double di = eph->Cic * cos_2phi + eph->Cis * sin_2phi;

    /* 改正后的升交角距、向径、倾角 */
    double u = phi + du;
    double r = a * (1.0 - eph->e * cos_E) + dr;
    double i = eph->i0 + eph->i_dot * tk + di;

    /* 轨道面坐标 */
    double xk = r * cos(u);
    double yk = r * sin(u);

    /* 升交点经度（含地球自转修正，实验五） */
    double Omega = eph->Omega0 +
                   (eph->Omega_dot - OMEGA_E) * tk -
                   OMEGA_E * eph->toe;

    /* ECEF坐标 */
    double sin_Omega = sin(Omega);
    double cos_Omega = cos(Omega);
    double sin_i = sin(i);
    double cos_i = cos(i);

    *x = xk * cos_Omega - yk * cos_i * sin_Omega;
    *y = xk * sin_Omega + yk * cos_i * cos_Omega;
    *z = yk * sin_i;

    /* 钟差：二阶多项式（实验五） */
    *clk_bias = eph->a0 + eph->a1 * tk + eph->a2 * tk * tk;

    /* 相对论效应改正: dt_rel = -2 * sqrt(GM*a) * e * sin(E) / c^2
       （实验五严格公式） */
    *clk_rel = -2.0 * sqrt(GM * a) * eph->e * sin_E /
               (C_LIGHT * C_LIGHT);
}

/* 为给定PRN和GPS秒找最近星历索引（实验五、九） */
int find_ephemeris(int prn, double t) {
    int    best_idx  = -1;
    double best_diff = 1e30;
    int    full_week = (int)(t / 604800.0);

    for (int i = 0; i < kepler_cnt; i++) {
        KeplerEphem* e = &kepler_list[i];
        if (e->sys != 'G' || e->prn != prn) continue;

        /* GPS 周数翻转解析 */
        int toe_week_resolved = e->toe_week;
        if (toe_week_resolved < full_week - 512)
            toe_week_resolved += 1024;
        else if (toe_week_resolved > full_week + 512)
            toe_week_resolved -= 1024;

        double toe_full = toe_week_resolved * 604800.0 + e->toe;
        double t_diff = fabs(t - toe_full);
        if (t_diff < 86400.0 && t_diff < best_diff) {
            best_diff = t_diff;
            best_idx  = i;
        }
    }
    return best_idx;
}

/* ========== 步骤4: 卫星速度计算（实验十，中心差分法） ========== */
/**
 * 使用中心差分计算卫星在ECEF中的速度：
 * v(t) ≈ [r(t+0.001) - r(t-0.001)] / 0.002
 */
void compute_sat_vel(KeplerEphem* eph, double t,
                     double* vx, double* vy, double* vz) {
    double xp, yp, zp, xm, ym, zm;
    double clk_dummy, clk_rel_dummy;
    double dt = 0.001;  /* 1ms 步长 */

    compute_sat_pos_clock(eph, t + dt, &xp, &yp, &zp, &clk_dummy, &clk_rel_dummy);
    compute_sat_pos_clock(eph, t - dt, &xm, &ym, &zm, &clk_dummy, &clk_rel_dummy);

    *vx = (xp - xm) / (2.0 * dt);
    *vy = (yp - ym) / (2.0 * dt);
    *vz = (zp - zm) / (2.0 * dt);
}

/* 卫星钟速（实验十）
   dt_s = a1 + 2*a2*(t - toc)  (s/s) */
double compute_sat_clock_drift(KeplerEphem* eph, double t) {
    double toc = calc_gps_seconds(eph->year, eph->month, eph->day,
        eph->hour, eph->minute, eph->second);
    return eph->a1 + 2.0 * eph->a2 * (t - toc);
}

/* ========== 矩阵运算（实验四） ========== */

/* 4x4 矩阵求逆（伴随矩阵法，实验四） */
int mat_inv_4x4(const double* M, double* Inv) {
    double det =
        M[0]  * (M[5]  * (M[10] * M[15] - M[11] * M[14]) -
                 M[6]  * (M[9]  * M[15] - M[11] * M[13]) +
                 M[7]  * (M[9]  * M[14] - M[10] * M[13])) -
        M[1]  * (M[4]  * (M[10] * M[15] - M[11] * M[14]) -
                 M[6]  * (M[8]  * M[15] - M[11] * M[12]) +
                 M[7]  * (M[8]  * M[14] - M[10] * M[12])) +
        M[2]  * (M[4]  * (M[9]  * M[15] - M[11] * M[13]) -
                 M[5]  * (M[8]  * M[15] - M[11] * M[12]) +
                 M[7]  * (M[8]  * M[13] - M[9]  * M[12])) -
        M[3]  * (M[4]  * (M[9]  * M[14] - M[10] * M[13]) -
                 M[5]  * (M[8]  * M[14] - M[10] * M[12]) +
                 M[6]  * (M[8]  * M[13] - M[9]  * M[12]));

    if (fabs(det) < 1e-15) return 0;

    double adj[16];
    adj[0]  =  (M[5]  * (M[10] * M[15] - M[11] * M[14]) -
                M[6]  * (M[9]  * M[15] - M[11] * M[13]) +
                M[7]  * (M[9]  * M[14] - M[10] * M[13]));
    adj[1]  = -(M[1]  * (M[10] * M[15] - M[11] * M[14]) -
                M[2]  * (M[9]  * M[15] - M[11] * M[13]) +
                M[3]  * (M[9]  * M[14] - M[10] * M[13]));
    adj[2]  =  (M[1]  * (M[6]  * M[15] - M[7]  * M[14]) -
                M[2]  * (M[5]  * M[15] - M[7]  * M[13]) +
                M[3]  * (M[5]  * M[14] - M[6]  * M[13]));
    adj[3]  = -(M[1]  * (M[6]  * M[11] - M[7]  * M[10]) -
                M[2]  * (M[5]  * M[11] - M[7]  * M[9])  +
                M[3]  * (M[5]  * M[10] - M[6]  * M[9]));
    adj[4]  = -(M[4]  * (M[10] * M[15] - M[11] * M[14]) -
                M[6]  * (M[8]  * M[15] - M[11] * M[12]) +
                M[7]  * (M[8]  * M[14] - M[10] * M[12]));
    adj[5]  =  (M[0]  * (M[10] * M[15] - M[11] * M[14]) -
                M[2]  * (M[8]  * M[15] - M[11] * M[12]) +
                M[3]  * (M[8]  * M[14] - M[10] * M[12]));
    adj[6]  = -(M[0]  * (M[6]  * M[15] - M[7]  * M[14]) -
                M[2]  * (M[4]  * M[15] - M[7]  * M[12]) +
                M[3]  * (M[4]  * M[14] - M[6]  * M[12]));
    adj[7]  =  (M[0]  * (M[6]  * M[11] - M[7]  * M[10]) -
                M[2]  * (M[4]  * M[11] - M[7]  * M[8])  +
                M[3]  * (M[4]  * M[10] - M[6]  * M[8]));
    adj[8]  =  (M[4]  * (M[9]  * M[15] - M[11] * M[13]) -
                M[5]  * (M[8]  * M[15] - M[11] * M[12]) +
                M[7]  * (M[8]  * M[13] - M[9]  * M[12]));
    adj[9]  = -(M[0]  * (M[9]  * M[15] - M[11] * M[13]) -
                M[1]  * (M[8]  * M[15] - M[11] * M[12]) +
                M[3]  * (M[8]  * M[13] - M[9]  * M[12]));
    adj[10] =  (M[0]  * (M[5]  * M[15] - M[7]  * M[13]) -
                M[1]  * (M[4]  * M[15] - M[7]  * M[12]) +
                M[3]  * (M[4]  * M[13] - M[5]  * M[12]));
    adj[11] = -(M[0]  * (M[5]  * M[11] - M[7]  * M[9])  -
                M[1]  * (M[4]  * M[11] - M[7]  * M[8])  +
                M[3]  * (M[4]  * M[9]  - M[5]  * M[8]));
    adj[12] = -(M[4]  * (M[9]  * M[14] - M[10] * M[13]) -
                M[5]  * (M[8]  * M[14] - M[10] * M[12]) +
                M[6]  * (M[8]  * M[13] - M[9]  * M[12]));
    adj[13] =  (M[0]  * (M[9]  * M[14] - M[10] * M[13]) -
                M[1]  * (M[8]  * M[14] - M[10] * M[12]) +
                M[2]  * (M[8]  * M[13] - M[9]  * M[12]));
    adj[14] = -(M[0]  * (M[5]  * M[14] - M[6]  * M[13]) -
                M[1]  * (M[4]  * M[14] - M[6]  * M[12]) +
                M[2]  * (M[4]  * M[13] - M[5]  * M[12]));
    adj[15] =  (M[0]  * (M[5]  * M[10] - M[6]  * M[9])  -
                M[1]  * (M[4]  * M[10] - M[6]  * M[8])  +
                M[2]  * (M[4]  * M[9]  - M[5]  * M[8]));

    for (int i = 0; i < 16; i++) Inv[i] = adj[i] / det;
    return 1;
}

/* ========== 坐标转换（实验三） ========== */

/* XYZ -> BLH (WGS84)（实验三的迭代法） */
void xyz2blh(double X, double Y, double Z,
             double* B, double* L, double* H) {
    double p = sqrt(X * X + Y * Y);
    *L = atan2(Y, X);
    if (*L < 0) *L += 2.0 * M_PI;

    double a  = WGS84_A;
    double e2 = WGS84_E2;

    /* 初值 */
    double tanB_val = Z / p * (1.0 + e2 * a /
                        sqrt(p * p + Z * Z * (1.0 - e2)));
    *B = atan(tanB_val);

    /* 迭代求解（实验三：5次迭代） */
    double N;
    for (int iter = 0; iter < 5; iter++) {
        N  = a / sqrt(1.0 - e2 * sin(*B) * sin(*B));
        *B = atan((Z + N * e2 * sin(*B)) / p);
    }
    N  = a / sqrt(1.0 - e2 * sin(*B) * sin(*B));
    *H = p / cos(*B) - N;
}

/* XYZ差值转到ENU（实验三的旋转矩阵法） */
void dxyz2enu(double dx, double dy, double dz,
              double B, double L,
              double* dE, double* dN, double* dU) {
    double sB = sin(B), cB = cos(B);
    double sL = sin(L), cL = cos(L);

    /* 旋转矩阵 R: [E,N,U]^T = R * [dX,dY,dZ]^T
       （实验三）
       R = [-sL,      cL,      0   ]
           [-sB*cL,  -sB*sL,   cB  ]
           [ cB*cL,   cB*sL,   sB  ] */
    *dE = -sL * dx + cL * dy;
    *dN = -sB * cL * dx - sB * sL * dy + cB * dz;
    *dU =  cB * cL * dx + cB * sL * dy + sB * dz;
}

/* 构建 ENU 旋转矩阵 R(3x3)（实验三、十）
   [E,N,U]^T = R * [dX,dY,dZ]^T
   R = [-sL,       cL,       0   ]
       [-sB*cL,   -sB*sL,    cB  ]
       [ cB*cL,    cB*sL,    sB  ] */
void build_enu_rot(double B, double L, double R[9]) {
    double sB = sin(B), cB = cos(B);
    double sL = sin(L), cL = cos(L);

    R[0] = -sL;      R[1] = cL;       R[2] = 0.0;
    R[3] = -sB * cL; R[4] = -sB * sL;  R[5] = cB;
    R[6] =  cB * cL; R[7] =  cB * sL;  R[8] = sB;
}

/* 3x3 矩阵乘法 C = A * B（实验十） */
void mat_mul_3x3(const double A[9], const double B[9], double C[9]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            C[i * 3 + j] = 0.0;
            for (int k = 0; k < 3; k++)
                C[i * 3 + j] += A[i * 3 + k] * B[k * 3 + j];
        }
    }
}

/* 3x3 矩阵转置（实验十） */
void mat_trans_3x3(const double A[9], double AT[9]) {
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            AT[i * 3 + j] = A[j * 3 + i];
}

/* ========== 主函数 ========== */
int main() {
    FILE* fp_out = fopen("output.txt", "w");
    if (!fp_out) {
        printf("无法创建 output.txt\n");
        return 1;
    }

    /* ---- 步骤1: 读取观测文件 ---- */
    printf("===== 读取观测文件 (RINEX 3.14) =====\n");
    if (!read_obs_file("jfng1590.24o")) {
        printf("观测文件读取失败！\n");
        fclose(fp_out);
        return 1;
    }

    /* ---- 步骤2: 读取广播星历文件 ---- */
    printf("===== 读取导航文件 (RINEX 3.14) =====\n");
    if (!read_nav_file("brdc1590.24p")) {
        printf("导航文件读取失败！\n");
        fclose(fp_out);
        return 1;
    }

    if (epoch_cnt == 0) {
        printf("无有效历元数据！\n");
        fclose(fp_out);
        return 1;
    }

    /* 测站大地坐标（用于高度角和ENU转换） */
    double B_sta, L_sta, H_sta;
    xyz2blh(true_x, true_y, true_z, &B_sta, &L_sta, &H_sta);
    printf("测站大地坐标: B=%.6f° L=%.6f° H=%.3f\n",
           B_sta * 180.0 / M_PI, L_sta * 180.0 / M_PI, H_sta);

    /* 频率相关常数 */
    double f1_sq = GPS_F1 * GPS_F1;
    double f2_sq = GPS_F2 * GPS_F2;
    double f_denom = f1_sq - f2_sq;
    double IF_coef1 = f1_sq / f_denom;   /* f1²/(f1²-f2²) */
    double IF_coef2 = f2_sq / f_denom;   /* f2²/(f1²-f2²) */

    /* 上历元速度初值（用于加速收敛） */
    double prev_vr_x = 0.0, prev_vr_y = 0.0, prev_vr_z = 0.0;
    double prev_cdt_r = 0.0;

    /* 输出文件表头 */
    fprintf(fp_out, "%-22s %-5s %-4s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s\n",
            "时间", "系统", "Sat", "PDOP", "dE(m)", "dN(m)", "dU(m)",
            "Vx(m/s)", "Vy(m/s)", "Vz(m/s)", "sigmaE", "sigmaN", "sigmaU");
    fprintf(fp_out,
            "----------------------------------------------------------------------------------------------------------------------------\n");

    printf("\n===== 开始历元循环处理 =====\n");

    /* ---- 步骤3: 历元循环 ---- */
    for (int iep = 0; iep < epoch_cnt; iep++) {
        ObsEpoch* ep = &epoch_list[iep];

        /* 当前历元GPS接收时刻（总秒数） */
        double t_recv_gps = calc_gps_seconds(ep->year, ep->month, ep->day,
                                              ep->hour, ep->minute, ep->second);

        printf("\n[历元 %d/%d] %04d-%02d-%02d %02d:%02d:%05.2f GPS时\n",
               iep + 1, epoch_cnt, ep->year, ep->month, ep->day,
               ep->hour, ep->minute, ep->second);

        /* ---- 预计算每颗可见卫星的轨道和钟差 ---- */
        double sat_x[MAX_GPS_SATS], sat_y[MAX_GPS_SATS], sat_z[MAX_GPS_SATS];
        double sat_vx[MAX_GPS_SATS], sat_vy[MAX_GPS_SATS], sat_vz[MAX_GPS_SATS];
        double sat_clk[MAX_GPS_SATS], sat_rel[MAX_GPS_SATS];
        double sat_clk_drift[MAX_GPS_SATS];   /* 卫星钟速 (s/s) */
        double P1[MAX_GPS_SATS], P2[MAX_GPS_SATS];
        double D1[MAX_GPS_SATS], D2[MAX_GPS_SATS];
        int    eph_idx[MAX_GPS_SATS];
        int    valid_sat  = 0;
        int    valid_prn[MAX_GPS_SATS];

        for (int i = 0; i < ep->nsat; i++) {
            int prn = ep->sats[i].prn;
            eph_idx[valid_sat] = find_ephemeris(prn, t_recv_gps);

            if (eph_idx[valid_sat] == -1) {
                printf("  PRN G%02d: 未找到有效星历，跳过\n", prn);
                continue;
            }

            KeplerEphem* eph = &kepler_list[eph_idx[valid_sat]];

            /* ---- 迭代求卫星位置和钟差（含钟差修正） ---- */
            /* 初始信号发射时间: t_send = t_recv - 0.067 (≈20000km/c) */
            double t_send = t_recv_gps - 0.067;
            double sx, sy, sz, cb, cr;

            for (int iter = 0; iter < 3; iter++) {
                compute_sat_pos_clock(eph, t_send, &sx, &sy, &sz, &cb, &cr);

                /* 计算站星距（用真值坐标近似） */
                double dx = sx - true_x;
                double dy = sy - true_y;
                double dz = sz - true_z;
                double rho = sqrt(dx * dx + dy * dy + dz * dz);

                /* t_send = t_recv - rho/c - dt_sat
                   dt_sat = cb（多项式部分）+ cr（相对论效应部分） */
                double dt_sat = cb + cr;
                t_send = t_recv_gps - rho / C_LIGHT - dt_sat;
            }

            /* 最后一次用收敛后的发射时刻计算位置、钟差 */
            compute_sat_pos_clock(eph, t_send, &sx, &sy, &sz, &cb, &cr);

            /* 卫星速度（中心差分法，实验十） */
            double svx, svy, svz;
            compute_sat_vel(eph, t_send, &svx, &svy, &svz);

            /* 卫星钟速（实验十） */
            double dt_s = compute_sat_clock_drift(eph, t_send);

            sat_x[valid_sat]   = sx;
            sat_y[valid_sat]   = sy;
            sat_z[valid_sat]   = sz;
            sat_vx[valid_sat]  = svx;
            sat_vy[valid_sat]  = svy;
            sat_vz[valid_sat]  = svz;
            sat_clk[valid_sat] = cb;
            sat_rel[valid_sat] = cr;
            sat_clk_drift[valid_sat] = dt_s;

            /* C1->P1 码偏差改正 */
            double p1_raw = ep->sats[i].freq[0].C;
            double p2_raw = ep->sats[i].freq[1].C;
            char   code1  = ep->sats[i].freq[0].code;
            char   code2  = ep->sats[i].freq[1].code;

            /* 若L1为C/A码(C1C)，用TGD改正为P1
               P1 = C1 + TGD * c * (f1^2 / (f1^2 - f2^2)) */
            if (code1 == 'C') {
                double tgd_scalar = f1_sq / (f1_sq - f2_sq);
                p1_raw += eph->tgd * C_LIGHT * tgd_scalar;
            }

            /* L2若为C2码，类似改正 */
            if (code2 == 'C') {
                double tgd_scalar2 = f1_sq / (f1_sq - f2_sq);
                p2_raw += eph->tgd * C_LIGHT * tgd_scalar2;
            }

            P1[valid_sat]       = p1_raw;
            P2[valid_sat]       = p2_raw;
            D1[valid_sat]       = ep->sats[i].freq[0].D;
            D2[valid_sat]       = ep->sats[i].freq[1].D;
            valid_prn[valid_sat] = prn;
            valid_sat++;
        }

        if (valid_sat < 4) {
            printf("  可见卫星不足4颗（%d颗），跳过\n", valid_sat);
            fprintf(fp_out, "%04d-%02d-%02d %02d:%02d:%05.2f %-5s %-4s"
                    " %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s\n",
                    ep->year, ep->month, ep->day,
                    ep->hour, ep->minute, ep->second,
                    "GPS", "0", "N/A", "N/A", "N/A", "N/A",
                    "N/A", "N/A", "N/A", "N/A", "N/A", "N/A");
            continue;
        }

        /* ============================================================ */
        /*  Part A: 伪距单点定位 (SPP)（实验九）                         */
        /* ============================================================ */

        /* 测站坐标初值（使用真值作为初始近似） */
        double X = true_x, Y = true_y, Z = true_z;
        double dt_rcv = 0.0;  /* 接收机钟差 c*δt (米) */

        int spp_converged = 0;
        double pdop = 0.0;
        double dE = 0.0, dN = 0.0, dU = 0.0;
        int n_eq_spp = 0;

        for (int iter = 0; iter < 10; iter++) {
            int    n_eq = 0;
            double A_vals[MAX_GPS_SATS * 4] = { 0 };
            double b_vals[MAX_GPS_SATS]     = { 0 };
            double w_vals[MAX_GPS_SATS]     = { 0 };

            /* 当前迭代点的大地坐标（用于高度角计算） */
            double B_cur, L_cur, H_cur;
            xyz2blh(X, Y, Z, &B_cur, &L_cur, &H_cur);
            double sB_cur = sin(B_cur), cB_cur = cos(B_cur);
            double sL_cur = sin(L_cur), cL_cur = cos(L_cur);

            for (int i = 0; i < valid_sat; i++) {
                /* 无电离层组合 */
                double PC = (f1_sq * P1[i] - f2_sq * P2[i]) /
                            (f1_sq - f2_sq);

                /* 地球自转改正（Sagnac效应）
                   tau = PC / c, α = ωe * tau */
                double tau   = PC / C_LIGHT;
                double alpha = OMEGA_E * tau;
                double cos_a = cos(alpha), sin_a = sin(alpha);
                double xs_corr = sat_x[i] * cos_a + sat_y[i] * sin_a;
                double ys_corr = -sat_x[i] * sin_a + sat_y[i] * cos_a;
                double zs_corr = sat_z[i];

                /* 星地几何距离 */
                double dx  = xs_corr - X;
                double dy  = ys_corr - Y;
                double dz  = zs_corr - Z;
                double rho = sqrt(dx * dx + dy * dy + dz * dz);
                double ux  = dx / rho;
                double uy  = dy / rho;
                double uz  = dz / rho;

                /* 高度角（用当前估计位置计算） */
                double e_enu = -sL_cur * dx + cL_cur * dy;
                double n_enu = -sB_cur * cL_cur * dx
                               - sB_cur * sL_cur * dy + cB_cur * dz;
                double u_enu = cB_cur * cL_cur * dx
                               + cB_cur * sL_cur * dy + sB_cur * dz;
                double rho_check = sqrt(e_enu * e_enu +
                                        n_enu * n_enu +
                                        u_enu * u_enu);
                double E_rad = 0.01;
                if (rho_check > 1e-12)
                    E_rad = asin(u_enu / rho_check);
                double E_deg = E_rad * 180.0 / M_PI;

                /* 截止高度角 10° */
                if (E_deg < 10.0) continue;

                /* 对流层改正
                   天顶干延迟 ZHD=2.3m
                   映射函数 M = 1.001/sqrt(0.002001+sin²E) */
                double ZHD    = 2.3;
                double M_wet  = 1.001 / sqrt(0.002001 +
                                   sin(E_rad) * sin(E_rad));
                double d_trop = ZHD * M_wet;

                /* 构建误差方程
                   观测方程: PC = ρ + c*dt_rcv - c*dt_sat - c*dt_rel + d_trop
                   残差: b = PC - (ρ + c*dt_rcv + d_trop - c*dt_sat - c*dt_rel) */
                double rho_calc = rho + dt_rcv + d_trop
                                - sat_clk[i] * C_LIGHT
                                - sat_rel[i] * C_LIGHT;
                double residual = PC - rho_calc;

                /* 设计矩阵
                   A_row = [-ux, -uy, -uz, 1]  (第4列对应c·δt) */
                A_vals[n_eq * 4 + 0] = -ux;
                A_vals[n_eq * 4 + 1] = -uy;
                A_vals[n_eq * 4 + 2] = -uz;
                A_vals[n_eq * 4 + 3] = 1.0;
                b_vals[n_eq] = residual;

                /* 高度角定权 σ² = (0.3/sinE)², w = 1/σ² */
                double sigma = 0.3 / sin(E_rad);
                w_vals[n_eq] = 1.0 / (sigma * sigma);
                n_eq++;
            }

            if (n_eq < 4) {
                printf("  SPP迭代%d: 可用卫星不足4颗(%d)，退出\n",
                       iter + 1, n_eq);
                break;
            }

            /* 加权最小二乘 N = A^T * W * A, U = A^T * W * b */
            double N[16] = { 0 };
            double U[4]  = { 0 };
            for (int i = 0; i < n_eq; i++) {
                double w  = w_vals[i];
                double a0 = A_vals[i * 4 + 0];
                double a1 = A_vals[i * 4 + 1];
                double a2 = A_vals[i * 4 + 2];
                double a3 = A_vals[i * 4 + 3];
                double bi = b_vals[i];

                N[0]  += w * a0 * a0;
                N[1]  += w * a0 * a1;
                N[2]  += w * a0 * a2;
                N[3]  += w * a0 * a3;
                N[4]  += w * a1 * a0;
                N[5]  += w * a1 * a1;
                N[6]  += w * a1 * a2;
                N[7]  += w * a1 * a3;
                N[8]  += w * a2 * a0;
                N[9]  += w * a2 * a1;
                N[10] += w * a2 * a2;
                N[11] += w * a2 * a3;
                N[12] += w * a3 * a0;
                N[13] += w * a3 * a1;
                N[14] += w * a3 * a2;
                N[15] += w * a3 * a3;

                U[0] += w * a0 * bi;
                U[1] += w * a1 * bi;
                U[2] += w * a2 * bi;
                U[3] += w * a3 * bi;
            }

            /* 求逆（伴随矩阵法） */
            double N_inv[16];
            if (!mat_inv_4x4(N, N_inv)) {
                printf("  SPP迭代%d: 矩阵不可逆，退出\n", iter + 1);
                break;
            }

            /* 求解 x = N^-1 * U */
            double dx_est[4] = { 0 };
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    dx_est[i] += N_inv[i * 4 + j] * U[j];

            /* 更新参数 */
            X      += dx_est[0];
            Y      += dx_est[1];
            Z      += dx_est[2];
            dt_rcv += dx_est[3];

            /* 检查收敛：无穷范数 < 1e-4 */
            double max_dx = 0.0;
            for (int i = 0; i < 4; i++) {
                double adx = fabs(dx_est[i]);
                if (adx > max_dx) max_dx = adx;
            }
            printf("  SPP迭代%d: dX=%.3f dY=%.3f dZ=%.3f cdt=%.3f "
                   "最大改正=%.6f\n",
                   iter + 1, dx_est[0], dx_est[1], dx_est[2],
                   dx_est[3], max_dx);

            if (max_dx < 1e-4) {
                spp_converged = 1;
                n_eq_spp = n_eq;

                /* 计算PDOP PDOP = sqrt(Qxx + Qyy + Qzz), Q = N^-1 */
                pdop = sqrt(N_inv[0] + N_inv[5] + N_inv[10]);

                /* 计算ENU误差 */
                double B_true, L_true, H_true;
                xyz2blh(true_x, true_y, true_z,
                        &B_true, &L_true, &H_true);
                double dx_err = X - true_x;
                double dy_err = Y - true_y;
                double dz_err = Z - true_z;
                dxyz2enu(dx_err, dy_err, dz_err,
                         B_true, L_true, &dE, &dN, &dU);

                printf("  => SPP定位结果: X=%.4f Y=%.4f Z=%.4f "
                       "PDOP=%.3f sat=%d\n",
                       X, Y, Z, pdop, n_eq);
                printf("  => SPP误差ENU: dE=%.3f dN=%.3f dU=%.3f m\n",
                       dE, dN, dU);
                break;
            }
        }

        if (!spp_converged) {
            printf("  SPP迭代未收敛（达到最大次数）\n");
        }

        /* ============================================================ */
        /*  Part B: 多普勒测速 (Velocity)（实验十）                      */
        /* ============================================================ */

        /* 使用SPP定位结果作为测站坐标（更精确），若SPP未收敛则用真值 */
        double sta_x_vel = spp_converged ? X : true_x;
        double sta_y_vel = spp_converged ? Y : true_y;
        double sta_z_vel = spp_converged ? Z : true_z;

        /* 速度初值 */
        double vr_x = prev_vr_x, vr_y = prev_vr_y, vr_z = prev_vr_z;
        double cdt_r = prev_cdt_r;  /* c * dt_r (m/s) */

        int vel_converged = 0;
        double fin_vr_x = vr_x, fin_vr_y = vr_y, fin_vr_z = vr_z;
        double fin_cdt_r = cdt_r;
        double Q_vel[9] = { 0 };
        double sigma0 = 0.0;
        double sigma_E = 0.0, sigma_N = 0.0, sigma_U = 0.0;

        /* 先计算地球自转改正后的卫星位置和视线向量 */
        double sat_xp[MAX_GPS_SATS], sat_yp[MAX_GPS_SATS], sat_zp[MAX_GPS_SATS];
        double sat_vxp[MAX_GPS_SATS], sat_vyp[MAX_GPS_SATS], sat_vzp[MAX_GPS_SATS];
        double e_x[MAX_GPS_SATS], e_y[MAX_GPS_SATS], e_z[MAX_GPS_SATS];
        double elev[MAX_GPS_SATS];
        double dopp_IF[MAX_GPS_SATS];
        int    n_vel = 0;
        int    vel_idx[MAX_GPS_SATS];

        for (int i = 0; i < valid_sat; i++) {
            /* 信号传播时间（用SPP定位结果） */
            double dx1 = sat_x[i] - sta_x_vel;
            double dy1 = sat_y[i] - sta_y_vel;
            double dz1 = sat_z[i] - sta_z_vel;
            double rho = sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
            double tau = rho / C_LIGHT;

            /* 地球自转角度 */
            double theta = OMEGA_E * tau;
            double cos_t = cos(theta), sin_t = sin(theta);

            /* 旋转卫星坐标 */
            double xs_p = sat_x[i] * cos_t + sat_y[i] * sin_t;
            double ys_p = -sat_x[i] * sin_t + sat_y[i] * cos_t;
            double zs_p = sat_z[i];

            /* 旋转卫星速度 */
            double vxs_p = sat_vx[i] * cos_t + sat_vy[i] * sin_t;
            double vys_p = -sat_vx[i] * sin_t + sat_vy[i] * cos_t;
            double vzs_p = sat_vz[i];

            /* 精化几何距离 */
            double dxp = xs_p - sta_x_vel;
            double dyp = ys_p - sta_y_vel;
            double dzp = zs_p - sta_z_vel;
            double rho_p = sqrt(dxp * dxp + dyp * dyp + dzp * dzp);

            /* 视线单位向量 (station -> satellite) */
            double ex = dxp / rho_p;
            double ey = dyp / rho_p;
            double ez = dzp / rho_p;

            /* 转换到站心 ENU 计算高度角 */
            double sB = sin(B_sta), cB = cos(B_sta);
            double sL = sin(L_sta), cL = cos(L_sta);

            double e_E = -sL * dxp + cL * dyp;
            double e_N = -sB * cL * dxp - sB * sL * dyp + cB * dzp;
            double e_U = cB * cL * dxp + cB * sL * dyp + sB * dzp;
            double rho_check = sqrt(e_E * e_E + e_N * e_N + e_U * e_U);
            double E_rad = asin(e_U / (rho_check > 1e-12 ? rho_check : 1.0));

            /* 截止高度角 10° */
            if (E_rad < 10.0 * M_PI / 180.0) continue;

            /* 无电离层组合多普勒
               dρ_IF = (f1²·dρ_L1 - f2²·dρ_L2) / (f1² - f2²)
               dρ = -c/f * D_Hz */
            double drho_L1 = -C_LIGHT / GPS_F1 * D1[i];
            double drho_L2 = -C_LIGHT / GPS_F2 * D2[i];
            double drho_IF = IF_coef1 * drho_L1 - IF_coef2 * drho_L2;

            /* 存储 */
            sat_xp[n_vel] = xs_p;
            sat_yp[n_vel] = ys_p;
            sat_zp[n_vel] = zs_p;
            sat_vxp[n_vel] = vxs_p;
            sat_vyp[n_vel] = vys_p;
            sat_vzp[n_vel] = vzs_p;
            e_x[n_vel] = ex;
            e_y[n_vel] = ey;
            e_z[n_vel] = ez;
            elev[n_vel] = E_rad;
            dopp_IF[n_vel] = drho_IF;
            vel_idx[n_vel] = i;
            n_vel++;
        }

        if (n_vel >= 4) {
            for (int iter = 0; iter < 15; iter++) {
                double N_v[16] = { 0 };
                double U_v[4]  = { 0 };

                for (int i = 0; i < n_vel; i++) {
                    int orig_i = vel_idx[i];

                    /* 几何距离变化率（不含接收机钟速的几何项） */
                    double drho_geom = (sat_vxp[i] - vr_x) * e_x[i]
                                     + (sat_vyp[i] - vr_y) * e_y[i]
                                     + (sat_vzp[i] - vr_z) * e_z[i];

                    /* 理论值: dρ_theo = (v_s' - v_r)·e - c·dt_s + c·dt_r */
                    double drho_theo = drho_geom
                                     - C_LIGHT * sat_clk_drift[orig_i]
                                     + cdt_r;

                    /* 残差: l = dρ_IF - dρ_theo */
                    double l_val = dopp_IF[i] - drho_theo;

                    /* 权: p_i = sin²(E) */
                    double p = sin(elev[i]) * sin(elev[i]);

                    /* 设计矩阵行: B_i = [-e_x, -e_y, -e_z, 1.0] */
                    double b0 = -e_x[i];
                    double b1 = -e_y[i];
                    double b2 = -e_z[i];
                    double b3 = 1.0;

                    /* N = B^T P B */
                    N_v[0]  += p * b0 * b0;
                    N_v[1]  += p * b0 * b1;
                    N_v[2]  += p * b0 * b2;
                    N_v[3]  += p * b0 * b3;
                    N_v[4]  += p * b1 * b0;
                    N_v[5]  += p * b1 * b1;
                    N_v[6]  += p * b1 * b2;
                    N_v[7]  += p * b1 * b3;
                    N_v[8]  += p * b2 * b0;
                    N_v[9]  += p * b2 * b1;
                    N_v[10] += p * b2 * b2;
                    N_v[11] += p * b2 * b3;
                    N_v[12] += p * b3 * b0;
                    N_v[13] += p * b3 * b1;
                    N_v[14] += p * b3 * b2;
                    N_v[15] += p * b3 * b3;

                    /* U = B^T P l */
                    U_v[0] += p * b0 * l_val;
                    U_v[1] += p * b1 * l_val;
                    U_v[2] += p * b2 * l_val;
                    U_v[3] += p * b3 * l_val;
                }

                /* 求逆 */
                double N_v_inv[16];
                if (!mat_inv_4x4(N_v, N_v_inv)) {
                    break;
                }

                /* 求解 ΔX = N⁻¹ U */
                double dx_v[4] = { 0 };
                for (int r = 0; r < 4; r++)
                    for (int c = 0; c < 4; c++)
                        dx_v[r] += N_v_inv[r * 4 + c] * U_v[c];

                /* 更新速度参数 */
                vr_x  += dx_v[0];
                vr_y  += dx_v[1];
                vr_z  += dx_v[2];
                cdt_r += dx_v[3];

                /* 收敛判断: ||Δv|| < 1e-4 */
                double norm_dv = sqrt(dx_v[0] * dx_v[0] + dx_v[1] * dx_v[1] + dx_v[2] * dx_v[2]);
                if (norm_dv < 1e-4) {
                    vel_converged = 1;

                    fin_vr_x = vr_x;
                    fin_vr_y = vr_y;
                    fin_vr_z = vr_z;
                    fin_cdt_r = cdt_r;

                    /* 精度评定 */
                    double vtpv = 0.0;
                    for (int i = 0; i < n_vel; i++) {
                        int orig_i = vel_idx[i];

                        double drho_geom = (sat_vxp[i] - vr_x) * e_x[i]
                                         + (sat_vyp[i] - vr_y) * e_y[i]
                                         + (sat_vzp[i] - vr_z) * e_z[i];
                        double drho_theo = drho_geom
                                         - C_LIGHT * sat_clk_drift[orig_i]
                                         + cdt_r;
                        double l_val = dopp_IF[i] - drho_theo;

                        double b0 = -e_x[i], b1 = -e_y[i], b2 = -e_z[i], b3 = 1.0;
                        double v_res = b0 * dx_v[0] + b1 * dx_v[1] + b2 * dx_v[2] + b3 * dx_v[3] - l_val;
                        double p = sin(elev[i]) * sin(elev[i]);
                        vtpv += p * v_res * v_res;
                    }
                    sigma0 = sqrt(vtpv / (n_vel - 4));

                    /* 协因数阵 Q_vel = Q(0..2, 0..2) */
                    for (int r = 0; r < 3; r++)
                        for (int c = 0; c < 3; c++)
                            Q_vel[r * 3 + c] = N_v_inv[r * 4 + c];

                    /* 构建 ENU 旋转矩阵，计算 ENU 速度中误差 */
                    double R_enu[9];
                    build_enu_rot(B_sta, L_sta, R_enu);
                    double RT[9];
                    mat_trans_3x3(R_enu, RT);
                    double temp[9], Q_ENU[9];
                    mat_mul_3x3(RT, Q_vel, temp);
                    mat_mul_3x3(temp, R_enu, Q_ENU);

                    sigma_E = sigma0 * sqrt(fabs(Q_ENU[0]));
                    sigma_N = sigma0 * sqrt(fabs(Q_ENU[4]));
                    sigma_U = sigma0 * sqrt(fabs(Q_ENU[8]));

                    printf("  => 测速结果: Vx=%.4f Vy=%.4f Vz=%.4f m/s\n",
                           fin_vr_x, fin_vr_y, fin_vr_z);
                    printf("  => ENU中误差: sigmaE=%.4f sigmaN=%.4f sigmaU=%.4f\n",
                           sigma_E, sigma_N, sigma_U);
                    break;
                }
            }

            if (vel_converged) {
                prev_vr_x = fin_vr_x;
                prev_vr_y = fin_vr_y;
                prev_vr_z = fin_vr_z;
                prev_cdt_r = fin_cdt_r;
            }
        }

        /* ============================================================ */
        /*  Part C: 输出结果                                              */
        /* ============================================================ */

        char time_str[22];
        sprintf(time_str, "%04d-%02d-%02d %02d:%02d:%05.2f",
                ep->year, ep->month, ep->day,
                ep->hour, ep->minute, ep->second);

        if (spp_converged && vel_converged) {
            fprintf(fp_out, "%-22s %-5s %-4d %-8.3f %-8.3f %-8.3f %-8.3f %-8.4f %-8.4f %-8.4f %-8.4f %-8.4f %-8.4f\n",
                    time_str, "GPS", n_eq_spp,
                    pdop, dE, dN, dU,
                    fin_vr_x, fin_vr_y, fin_vr_z,
                    sigma_E, sigma_N, sigma_U);
        }
        else if (spp_converged) {
            fprintf(fp_out, "%-22s %-5s %-4d %-8.3f %-8.3f %-8.3f %-8.3f %-8s %-8s %-8s %-8s %-8s %-8s\n",
                    time_str, "GPS", n_eq_spp,
                    pdop, dE, dN, dU,
                    "N/A", "N/A", "N/A", "N/A", "N/A", "N/A");
        }
        else {
            fprintf(fp_out, "%-22s %-5s %-4s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s\n",
                    time_str, "GPS", "0",
                    "N/A", "N/A", "N/A", "N/A",
                    "N/A", "N/A", "N/A", "N/A", "N/A", "N/A");
        }
    }

    /* ============================================================ */
    /*  Part D: 精度统计（RMS）                                       */
    /* ============================================================ */

    printf("\n===== 重新读取输出文件进行精度统计 =====\n");

    fclose(fp_out);

    /* 重新打开输出文件进行统计 */
    fp_out = fopen("output.txt", "r");
    if (!fp_out) {
        printf("无法打开 output.txt 进行统计\n");
        return 1;
    }

    /* 跳过表头（2行） */
    char line[MAX_LINE_LEN];
    fgets(line, MAX_LINE_LEN, fp_out);
    fgets(line, MAX_LINE_LEN, fp_out);

    double sum_dE = 0, sum_dN = 0, sum_dU = 0;
    double sum_dE2 = 0, sum_dN2 = 0, sum_dU2 = 0;
    double sum_vx = 0, sum_vy = 0, sum_vz = 0;
    double sum_vx2 = 0, sum_vy2 = 0, sum_vz2 = 0;
    int count_spp = 0, count_vel = 0;

    while (fgets(line, MAX_LINE_LEN, fp_out)) {
        /* 行格式: "YYYY-MM-DD HH:MM:SS.ss GPS N PDOP dE dN dU Vx Vy Vz sigmaE sigmaN sigmaU" */
        char date_buf[12], time_buf[14], sys_buf[8];
        int nsat;
        double pdop_val, de, dn, du, vx, vy, vz, se, sn, su;

        int n = sscanf(line, "%s %s %s %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                       date_buf, time_buf, sys_buf, &nsat,
                       &pdop_val, &de, &dn, &du,
                       &vx, &vy, &vz, &se, &sn, &su);

        if (n >= 7 && strcmp(sys_buf, "GPS") == 0 && nsat > 0) {
            /* SPP 统计 */
            sum_dE  += de;
            sum_dN  += dn;
            sum_dU  += du;
            sum_dE2 += de * de;
            sum_dN2 += dn * dn;
            sum_dU2 += du * du;
            count_spp++;

            if (n >= 10 && vx > -9999.0) {
                /* Velocity 统计 */
                sum_vx  += vx;
                sum_vy  += vy;
                sum_vz  += vz;
                sum_vx2 += vx * vx;
                sum_vy2 += vy * vy;
                sum_vz2 += vz * vz;
                count_vel++;
            }
        }
    }

    fclose(fp_out);

    /* 追加统计到文件 */
    fp_out = fopen("output.txt", "a");
    if (fp_out) {
        fprintf(fp_out, "\n");
        fprintf(fp_out, "========================================\n");
        fprintf(fp_out, "           精度统计结果\n");
        fprintf(fp_out, "========================================\n");

        if (count_spp > 0) {
            double rms_E = sqrt(sum_dE2 / count_spp);
            double rms_N = sqrt(sum_dN2 / count_spp);
            double rms_U = sqrt(sum_dU2 / count_spp);
            double mean_E = sum_dE / count_spp;
            double mean_N = sum_dN / count_spp;
            double mean_U = sum_dU / count_spp;

            fprintf(fp_out, "\n--- 伪距单点定位精度 (SPP) ---\n");
            fprintf(fp_out, "有效历元数: %d\n", count_spp);
            fprintf(fp_out, "  ENU方向    Mean(m)    RMS(m)\n");
            fprintf(fp_out, "  E          %8.3f  %8.3f\n", mean_E, rms_E);
            fprintf(fp_out, "  N          %8.3f  %8.3f\n", mean_N, rms_N);
            fprintf(fp_out, "  U          %8.3f  %8.3f\n", mean_U, rms_U);
        }

        if (count_vel > 0) {
            double rms_vx = sqrt(sum_vx2 / count_vel);
            double rms_vy = sqrt(sum_vy2 / count_vel);
            double rms_vz = sqrt(sum_vz2 / count_vel);
            double mean_vx = sum_vx / count_vel;
            double mean_vy = sum_vy / count_vel;
            double mean_vz = sum_vz / count_vel;

            fprintf(fp_out, "\n--- 多普勒测速统计 (Velocity) ---\n");
            fprintf(fp_out, "有效历元数: %d\n", count_vel);
            fprintf(fp_out, "  XYZ方向    Mean(m/s)  RMS(m/s)\n");
            fprintf(fp_out, "  Vx         %8.4f  %8.4f\n", mean_vx, rms_vx);
            fprintf(fp_out, "  Vy         %8.4f  %8.4f\n", mean_vy, rms_vy);
            fprintf(fp_out, "  Vz         %8.4f  %8.4f\n", mean_vz, rms_vz);
        }

        fprintf(fp_out, "========================================\n");
        fclose(fp_out);
    }

    printf("\n===== 处理完成！结果已保存到 output.txt =====\n");
    printf("SPP有效历元: %d, 测速有效历元: %d\n", count_spp, count_vel);

    return 0;
}