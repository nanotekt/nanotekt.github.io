#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <vector>
#include <algorithm>
#include <random>
#include <thread>
#include <atomic>
#include <string>
#include <map>
#include <sys/stat.h>

struct Vec {
    double x, y, z;
    Vec(double x=0, double y=0, double z=0) : x(x), y(y), z(z) {}
    Vec operator+(const Vec& b) const { return {x+b.x, y+b.y, z+b.z}; }
    Vec operator-(const Vec& b) const { return {x-b.x, y-b.y, z-b.z}; }
    Vec operator*(double s) const { return {x*s, y*s, z*s}; }
    double dot(const Vec& b) const { return x*b.x + y*b.y + z*b.z; }
    Vec cross(const Vec& b) const { return {y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x}; }
    double len() const { return std::sqrt(x*x + y*y + z*z); }
    Vec norm() const { double d = len(); return d > 0 ? *this * (1.0/d) : Vec(); }
    Vec cmul(const Vec& b) const { return {x*b.x, y*b.y, z*b.z}; }
};

struct Mat3 {
    double m[3][3];
    Vec mul(const Vec& v) const {
        return Vec(m[0][0]*v.x + m[0][1]*v.y + m[0][2]*v.z,
                   m[1][0]*v.x + m[1][1]*v.y + m[1][2]*v.z,
                   m[2][0]*v.x + m[2][1]*v.y + m[2][2]*v.z);
    }
    Mat3 transpose() const {
        Mat3 r = {};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                r.m[i][j] = m[j][i];
        return r;
    }
};

static Mat3 make_rotation(double ax, double ay, double az) {
    double cx = std::cos(ax), sx = std::sin(ax);
    double cy = std::cos(ay), sy = std::sin(ay);
    double cz = std::cos(az), sz = std::sin(az);
    Mat3 r;
    r.m[0][0] = cy*cz; r.m[0][1] = cz*sx*sy - cx*sz; r.m[0][2] = cx*cz*sy + sx*sz;
    r.m[1][0] = cy*sz; r.m[1][1] = cx*cz + sx*sy*sz; r.m[1][2] = cx*sy*sz - cz*sx;
    r.m[2][0] = -sy;   r.m[2][1] = cy*sx;             r.m[2][2] = cx*cy;
    return r;
}

// --- Minimal YAML parser (no external deps) ---

struct YamlConfig {
    std::map<std::string, std::string> data;

    static std::string trim(const std::string& s) {
        size_t a = s.find_first_not_of(" \t\r\n");
        if (a == std::string::npos) return "";
        size_t b = s.find_last_not_of(" \t\r\n");
        return s.substr(a, b - a + 1);
    }

    bool load(const char* filename) {
        FILE* f = fopen(filename, "r");
        if (!f) return false;
        char buf[512];
        std::string prefix[8];
        while (fgets(buf, sizeof(buf), f)) {
            std::string line(buf);
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r'))
                line.pop_back();
            // Strip comments
            size_t hash = line.find('#');
            if (hash != std::string::npos)
                line = line.substr(0, hash);
            std::string trimmed = trim(line);
            if (trimmed.empty()) continue;
            // Indent level
            int indent = 0;
            for (char c : line) { if (c == ' ') indent++; else break; }
            int level = indent / 2;
            size_t colon = trimmed.find(':');
            if (colon == std::string::npos) continue;
            std::string key = trim(trimmed.substr(0, colon));
            std::string value = trim(trimmed.substr(colon + 1));
            prefix[level] = key;
            if (value.empty()) continue;  // parent key
            std::string full_key;
            for (int i = 0; i < level; i++) full_key += prefix[i] + ".";
            full_key += key;
            data[full_key] = value;
        }
        fclose(f);
        return true;
    }

    int get_int(const char* key, int def) const {
        auto it = data.find(key);
        if (it == data.end()) return def;
        return (int)std::strtol(it->second.c_str(), nullptr, 10);
    }

    double get_double(const char* key, double def) const {
        auto it = data.find(key);
        if (it == data.end()) return def;
        return std::strtod(it->second.c_str(), nullptr);
    }

    Vec get_vec(const char* key, Vec def) const {
        auto it = data.find(key);
        if (it == data.end()) return def;
        double x, y, z;
        if (sscanf(it->second.c_str(), "[%lf, %lf, %lf]", &x, &y, &z) == 3)
            return Vec(x, y, z);
        if (sscanf(it->second.c_str(), "[%lf,%lf,%lf]", &x, &y, &z) == 3)
            return Vec(x, y, z);
        return def;
    }

    std::string get_string(const char* key, const std::string& def = "") const {
        auto it = data.find(key);
        if (it == data.end()) return def;
        return it->second;
    }

    void get_int_array(const char* key, int* out, int count) const {
        auto it = data.find(key);
        if (it == data.end()) return;
        const char* p = it->second.c_str();
        if (*p == '[') p++;
        for (int i = 0; i < count; i++) {
            out[i] = (int)strtol(p, (char**)&p, 10);
            while (*p == ',' || *p == ' ') p++;
        }
    }
};

// --- Scene configuration ---

struct SceneConfig {
    // Rendering
    int width = 1024, height = 1024, samples = 4;
    double fog_max_dist = 60.0;

    // Camera
    Vec cam{0, 3.5, 14}, target{0, 3, -12};

    // Text
    std::vector<std::string> bitmap_rows;
    double sph_r = 1.0;
    double text_x_offset = -8.5, text_y_base = 8.8, text_y_spacing = 1.0, text_z = -12.0;
    double slash_x_base = -6.5, slash_y_base = 6.3;
    int slash_count = 7;

    // Rounded cube
    Vec rcube_center{0.0, 8.3, -20.0};
    double rcube_half = 4.0, rcube_round = 0.75;
    Vec rcube_rotation{0.45, 0.65, 0.3};

    // Red disc
    Vec disc_center{0.0, 8.3, -28.0};
    Vec disc_normal{0.0, 0.0, 1.0};
    double disc_r = 12.0;
    Vec disc_emissive{4.0, 0.04, 0.04};

    // Floor
    int tile_lo_x = -15, tile_hi_x = 15, tile_lo_z = -25, tile_hi_z = 15;
    double tile_max_y = 1.25;
    Vec floor_white{0.95, 0.95, 0.95}, floor_red{0.85, 0.08, 0.08};
    double floor_diffuse = 0.5, floor_ambient = 0.12;
    double floor_refl_base = 0.7, floor_refl_mirror = 0.3;
    int floor_max_depth = 2;
    double floor_spec_power = 48, floor_spec_intensity = 0.4;

    // Glass
    double glass_edge = 1.2, glass_ior = 1.5, glass_r0 = 0.04;
    double glass_tint_scale = 0.15;
    Vec glass_tint_rgb{0.05, 0.02, 0.01};
    double glass_spec_power = 96, glass_spec_intensity = 0.3;
    int glass_max_depth = 4, glass_passthrough_depth = 6;

    // Lights
    Vec light1_pos{6, 12, -4}, light1_color{0.4, 0.6, 1.0};
    Vec light2_pos{-6, 12, 4}, light2_color{1.0, 0.9, 0.3};
    Vec light3_color{2.0, 0.04, 0.04};
    double shadow_ambient = 0.15;

    // Chrome material
    Vec chrome_color{0.85, 0.85, 0.9};
    double chrome_diffuse = 0.3, chrome_ambient = 0.08;
    double chrome_refl_base = 0.3, chrome_refl_mirror = 0.7;
    int chrome_max_depth = 3;
    double chrome_spec_power = 64, chrome_spec_intensity = 0.6;

    // Silver material
    Vec silver_color{0.9, 0.9, 0.95};
    Vec scratched_color{0.45, 0.45, 0.5};
    double scratch_discolor = 2.5;
    double silver_diffuse = 0.5, silver_ambient = 0.15;
    double silver_refl_base = 0.55, silver_refl_mirror = 0.45;
    int silver_max_depth = 3;
    double silver_spec_power = 32, silver_spec_intensity = 0.5;

    // Scratches
    double scratch_uv = 8.0;
    double s1_freq = 12.0, s1_power = 4.0, s1_intensity = 0.35;
    double s2_angle = 1.0, s2_freq = 20.0, s2_phase = 1.7, s2_power = 5.0, s2_intensity = 0.2;
    double s3_angle = 2.3, s3_freq = 30.0, s3_phase = 0.6, s3_power = 6.0, s3_intensity = 0.12;

    // Sky
    Vec sky_base{0, 0, 0};
    double sky_horizon = 0.005, ceiling_h = 35.0, grid_spacing = 10.0;
    double glow_width = 0.6, grid_fade = 0.015;
    Vec neon_color{1.2, 0.2, 0.6};
    double neon_intensity = 1.8;

    // Fog
    double fog_density = 0.015;
    int fog_steps = 32;
    double fog_nx = 0.03, fog_nxz = 0.015, fog_ny = 0.05, fog_nz = 0.04;

    // Bloom
    int bloom_r = 20;
    double bloom_threshold = 0.8, bloom_intensity = 0.6, gamma = 0.45;

    // Derived values
    double sph_r2, rcube_bound_r, rcube_bound_r2, disc_r2;
    double glass_max_y;
    int tile_nx, tile_nz;
    Mat3 rcube_rot, rcube_rot_inv;

    void compute_derived() {
        sph_r2 = sph_r * sph_r;
        rcube_bound_r = std::sqrt(3.0) * rcube_half + rcube_round;
        rcube_bound_r2 = rcube_bound_r * rcube_bound_r;
        disc_r2 = disc_r * disc_r;
        glass_max_y = tile_max_y + glass_edge;
        tile_nx = tile_hi_x - tile_lo_x + 1;
        tile_nz = tile_hi_z - tile_lo_z + 1;
        rcube_rot = make_rotation(rcube_rotation.x, rcube_rotation.y, rcube_rotation.z);
        rcube_rot_inv = rcube_rot.transpose();
    }

    void load_from_yaml(const YamlConfig& cfg) {
        width = cfg.get_int("rendering.width", width);
        height = cfg.get_int("rendering.height", height);
        samples = cfg.get_int("rendering.samples", samples);
        fog_max_dist = cfg.get_double("rendering.fog_max_dist", fog_max_dist);

        cam = cfg.get_vec("camera.position", cam);
        target = cfg.get_vec("camera.target", target);

        bitmap_rows.clear();
        for (int i = 0; ; i++) {
            char key[64];
            snprintf(key, sizeof(key), "text.bitmap.%d", i);
            std::string row = cfg.get_string(key);
            if (row.empty()) break;
            bitmap_rows.push_back(row);
        }
        sph_r = cfg.get_double("text.sphere_radius", sph_r);
        text_x_offset = cfg.get_double("text.x_offset", text_x_offset);
        text_y_base = cfg.get_double("text.y_base", text_y_base);
        text_y_spacing = cfg.get_double("text.y_spacing", text_y_spacing);
        text_z = cfg.get_double("text.z_depth", text_z);
        slash_x_base = cfg.get_double("text.slash_x_base", slash_x_base);
        slash_y_base = cfg.get_double("text.slash_y_base", slash_y_base);
        slash_count = cfg.get_int("text.slash_count", slash_count);

        rcube_center = cfg.get_vec("rounded_cube.center", rcube_center);
        rcube_half = cfg.get_double("rounded_cube.half_edge", rcube_half);
        rcube_round = cfg.get_double("rounded_cube.rounding", rcube_round);
        rcube_rotation = cfg.get_vec("rounded_cube.rotation", rcube_rotation);

        disc_center = cfg.get_vec("red_disc.center", disc_center);
        disc_normal = cfg.get_vec("red_disc.normal", disc_normal);
        disc_r = cfg.get_double("red_disc.radius", disc_r);
        disc_emissive = cfg.get_vec("red_disc.emissive_color", disc_emissive);

        Vec bx = cfg.get_vec("floor.bounds_x", Vec(tile_lo_x, tile_hi_x, 0));
        tile_lo_x = (int)bx.x; tile_hi_x = (int)bx.y;
        Vec bz = cfg.get_vec("floor.bounds_z", Vec(tile_lo_z, tile_hi_z, 0));
        tile_lo_z = (int)bz.x; tile_hi_z = (int)bz.y;
        tile_max_y = cfg.get_double("floor.max_height", tile_max_y);
        floor_white = cfg.get_vec("floor.white_color", floor_white);
        floor_red = cfg.get_vec("floor.red_color", floor_red);
        floor_diffuse = cfg.get_double("floor.diffuse", floor_diffuse);
        floor_ambient = cfg.get_double("floor.ambient", floor_ambient);
        floor_refl_base = cfg.get_double("floor.reflect_base", floor_refl_base);
        floor_refl_mirror = cfg.get_double("floor.reflect_mirror", floor_refl_mirror);
        floor_max_depth = cfg.get_int("floor.max_depth", floor_max_depth);
        floor_spec_power = cfg.get_double("floor.spec_power", floor_spec_power);
        floor_spec_intensity = cfg.get_double("floor.spec_intensity", floor_spec_intensity);

        glass_edge = cfg.get_double("glass.edge", glass_edge);
        glass_ior = cfg.get_double("glass.ior", glass_ior);
        glass_r0 = cfg.get_double("glass.fresnel_r0", glass_r0);
        glass_tint_scale = cfg.get_double("glass.tint_scale", glass_tint_scale);
        glass_tint_rgb = cfg.get_vec("glass.tint_rgb", glass_tint_rgb);
        glass_spec_power = cfg.get_double("glass.spec_power", glass_spec_power);
        glass_spec_intensity = cfg.get_double("glass.spec_intensity", glass_spec_intensity);
        glass_max_depth = cfg.get_int("glass.max_depth", glass_max_depth);
        glass_passthrough_depth = cfg.get_int("glass.passthrough_depth", glass_passthrough_depth);

        light1_pos = cfg.get_vec("lights.light1.position", light1_pos);
        light1_color = cfg.get_vec("lights.light1.color", light1_color);
        light2_pos = cfg.get_vec("lights.light2.position", light2_pos);
        light2_color = cfg.get_vec("lights.light2.color", light2_color);
        light3_color = cfg.get_vec("lights.light3_color", light3_color);
        shadow_ambient = cfg.get_double("lights.shadow_ambient", shadow_ambient);

        chrome_color = cfg.get_vec("materials.chrome.color", chrome_color);
        chrome_diffuse = cfg.get_double("materials.chrome.diffuse", chrome_diffuse);
        chrome_ambient = cfg.get_double("materials.chrome.ambient", chrome_ambient);
        chrome_refl_base = cfg.get_double("materials.chrome.reflect_base", chrome_refl_base);
        chrome_refl_mirror = cfg.get_double("materials.chrome.reflect_mirror", chrome_refl_mirror);
        chrome_max_depth = cfg.get_int("materials.chrome.max_depth", chrome_max_depth);
        chrome_spec_power = cfg.get_double("materials.chrome.spec_power", chrome_spec_power);
        chrome_spec_intensity = cfg.get_double("materials.chrome.spec_intensity", chrome_spec_intensity);

        silver_color = cfg.get_vec("materials.silver.color", silver_color);
        scratched_color = cfg.get_vec("materials.silver.scratched_color", scratched_color);
        scratch_discolor = cfg.get_double("materials.silver.discolor_scale", scratch_discolor);
        silver_diffuse = cfg.get_double("materials.silver.diffuse", silver_diffuse);
        silver_ambient = cfg.get_double("materials.silver.ambient", silver_ambient);
        silver_refl_base = cfg.get_double("materials.silver.reflect_base", silver_refl_base);
        silver_refl_mirror = cfg.get_double("materials.silver.reflect_mirror", silver_refl_mirror);
        silver_max_depth = cfg.get_int("materials.silver.max_depth", silver_max_depth);
        silver_spec_power = cfg.get_double("materials.silver.spec_power", silver_spec_power);
        silver_spec_intensity = cfg.get_double("materials.silver.spec_intensity", silver_spec_intensity);

        scratch_uv = cfg.get_double("scratches.uv_scale", scratch_uv);
        s1_freq = cfg.get_double("scratches.layer1.frequency", s1_freq);
        s1_power = cfg.get_double("scratches.layer1.power", s1_power);
        s1_intensity = cfg.get_double("scratches.layer1.intensity", s1_intensity);
        s2_angle = cfg.get_double("scratches.layer2.angle_offset", s2_angle);
        s2_freq = cfg.get_double("scratches.layer2.frequency", s2_freq);
        s2_phase = cfg.get_double("scratches.layer2.phase_mult", s2_phase);
        s2_power = cfg.get_double("scratches.layer2.power", s2_power);
        s2_intensity = cfg.get_double("scratches.layer2.intensity", s2_intensity);
        s3_angle = cfg.get_double("scratches.layer3.angle_offset", s3_angle);
        s3_freq = cfg.get_double("scratches.layer3.frequency", s3_freq);
        s3_phase = cfg.get_double("scratches.layer3.phase_mult", s3_phase);
        s3_power = cfg.get_double("scratches.layer3.power", s3_power);
        s3_intensity = cfg.get_double("scratches.layer3.intensity", s3_intensity);

        sky_base = cfg.get_vec("sky.base_color", sky_base);
        sky_horizon = cfg.get_double("sky.horizon_threshold", sky_horizon);
        ceiling_h = cfg.get_double("sky.ceiling_height", ceiling_h);
        grid_spacing = cfg.get_double("sky.grid_spacing", grid_spacing);
        glow_width = cfg.get_double("sky.glow_width", glow_width);
        grid_fade = cfg.get_double("sky.fade_rate", grid_fade);
        neon_color = cfg.get_vec("sky.neon_color", neon_color);
        neon_intensity = cfg.get_double("sky.neon_intensity", neon_intensity);

        fog_density = cfg.get_double("fog.density", fog_density);
        fog_steps = cfg.get_int("fog.steps", fog_steps);
        fog_nx = cfg.get_double("fog.noise_x", fog_nx);
        fog_nxz = cfg.get_double("fog.noise_xz", fog_nxz);
        fog_ny = cfg.get_double("fog.noise_y", fog_ny);
        fog_nz = cfg.get_double("fog.noise_z", fog_nz);

        bloom_r = cfg.get_int("bloom.radius", bloom_r);
        bloom_threshold = cfg.get_double("bloom.threshold", bloom_threshold);
        bloom_intensity = cfg.get_double("bloom.intensity", bloom_intensity);
        gamma = cfg.get_double("bloom.gamma", gamma);

        compute_derived();
    }
};

static SceneConfig SC;

// Sphere centers (text + slash)
struct SphCenter { double x, y, z; };
static std::vector<SphCenter> centers;

// Tile heights (dynamic — bounds come from YAML)
static std::vector<double> tile_heights;

static inline double get_tile_height(int ix, int iz) {
    int i = ix - SC.tile_lo_x;
    int j = iz - SC.tile_lo_z;
    if (i < 0 || i >= SC.tile_nx || j < 0 || j >= SC.tile_nz) return -1.0;
    return tile_heights[i * SC.tile_nz + j];
}

static void init_tiles() {
    tile_heights.resize(SC.tile_nx * SC.tile_nz);
    for (int ix = SC.tile_lo_x; ix <= SC.tile_hi_x; ix++) {
        for (int iz = SC.tile_lo_z; iz <= SC.tile_hi_z; iz++) {
            uint32_t h = ((uint32_t)((uint32_t)(ix * 374761393u) + (uint32_t)(iz * 668265263u)) ^ 0xDEADBEEFu);
            double height = (h & 0xFFFF) / 65535.0 * SC.tile_max_y;
            tile_heights[(ix - SC.tile_lo_x) * SC.tile_nz + (iz - SC.tile_lo_z)] = height;
        }
    }
}

static void init_centers() {
    for (int j = 0; j < (int)SC.bitmap_rows.size(); j++) {
        for (int i = 0; i < (int)SC.bitmap_rows[j].size(); i++) {
            if (SC.bitmap_rows[j][i] == '@') {
                centers.push_back({i + SC.text_x_offset, SC.text_y_base - j * SC.text_y_spacing, SC.text_z});
            }
        }
    }
    for (int k = 0; k < SC.slash_count; k++) {
        centers.push_back({SC.slash_x_base + (3 - k), SC.slash_y_base + (3 - k), SC.text_z});
    }
}

static inline double sdf_round_box(const Vec& p) {
    double qx = std::max(std::abs(p.x) - SC.rcube_half, 0.0);
    double qy = std::max(std::abs(p.y) - SC.rcube_half, 0.0);
    double qz = std::max(std::abs(p.z) - SC.rcube_half, 0.0);
    double outer = std::sqrt(qx*qx + qy*qy + qz*qz);
    double inner = std::min(std::max({std::abs(p.x) - SC.rcube_half, std::abs(p.y) - SC.rcube_half, std::abs(p.z) - SC.rcube_half}), 0.0);
    return outer + inner - SC.rcube_round;
}

static Vec sdf_normal(const Vec& p) {
    const double eps = 0.001;
    double dx = sdf_round_box(p + Vec(eps,0,0)) - sdf_round_box(p - Vec(eps,0,0));
    double dy = sdf_round_box(p + Vec(0,eps,0)) - sdf_round_box(p - Vec(0,eps,0));
    double dz = sdf_round_box(p + Vec(0,0,eps)) - sdf_round_box(p - Vec(0,0,eps));
    return Vec(dx, dy, dz).norm();
}

// Thread-local glass hit info
struct GlassInfo { int ix, iz; double th; };
static thread_local GlassInfo glass_hit = {0, 0, 0.0};

// Thread-local sphere index for scratch perturbation
static thread_local int hit_sphere_idx = 0;

struct HitResult {
    double t;
    int m;   // 0=sky, 1=white tile, 2=red tile, 3=text sphere, 4=glass, 5=rounded cube
    Vec n;
};

// AABB entry test (inlined). Returns entry t and sets normal. Returns false if no hit.
static inline bool aabb_test(double ox, double oy, double oz, double dx, double dy, double dz,
                              double x0, double x1, double y0, double y1, double z0, double z1,
                              bool has_dx, double inv_dx, double inv_dy,
                              bool has_dz, double inv_dz,
                              double max_t, double& out_t, Vec& out_n) {
    double tx0, tx1, ty0_, ty1_, tz0, tz1;
    if (has_dx) {
        tx0 = (x0 - ox) * inv_dx;
        tx1 = (x1 - ox) * inv_dx;
        if (tx0 > tx1) std::swap(tx0, tx1);
    } else {
        if (ox >= x0 && ox < x1) { tx0 = -1e18; tx1 = 1e18; }
        else return false;
    }
    ty0_ = (y0 - oy) * inv_dy;
    ty1_ = (y1 - oy) * inv_dy;
    if (ty0_ > ty1_) std::swap(ty0_, ty1_);
    if (has_dz) {
        tz0 = (z0 - oz) * inv_dz;
        tz1 = (z1 - oz) * inv_dz;
        if (tz0 > tz1) std::swap(tz0, tz1);
    } else {
        if (oz >= z0 && oz < z1) { tz0 = -1e18; tz1 = 1e18; }
        else return false;
    }
    double te = std::max({tx0, ty0_, tz0});
    double tx = std::min({tx1, ty1_, tz1});
    if (te >= tx || tx <= 0.01 || te >= max_t) return false;
    if (te < 0.01) te = 0.01;
    if (te >= max_t) return false;
    out_t = te;
    // Determine normal from which axis gave te
    if (te == tx0) out_n = Vec(dx > 0 ? -1 : 1, 0, 0);
    else if (te == ty0_) out_n = Vec(0, dy > 0 ? -1 : 1, 0);
    else out_n = Vec(0, 0, dz > 0 ? -1 : 1);
    return true;
}

static HitResult test(const Vec& o, const Vec& d, bool skip_tiles = false) {
    double t = 1e9;
    int m = 0;
    Vec n(0, 0, 1);
    double ox = o.x, oy = o.y, oz = o.z;
    double dx = d.x, dy = d.y, dz = d.z;

    // Floor cubes + glass cubes
    if (!skip_tiles && std::abs(dy) > 1e-6) {
        double inv_dy = 1.0 / dy;
        double t_y0 = -oy * inv_dy;
        double t_ymax = (SC.glass_max_y - oy) * inv_dy;
        double t_enter = std::min(t_y0, t_ymax);
        double t_exit_y = std::max(t_y0, t_ymax);
        if (t_exit_y > 0.01 && t_enter < t) {
            t_enter = std::max(t_enter, 0.01);
            double t_exit = std::min(t_exit_y, t);
            double x0 = ox + dx * t_enter, x1 = ox + dx * t_exit;
            double z0 = oz + dz * t_enter, z1 = oz + dz * t_exit;
            int ix_lo = std::max((int)std::floor(std::min(x0, x1)) - 1, SC.tile_lo_x);
            int ix_hi = std::min((int)std::floor(std::max(x0, x1)) + 1, SC.tile_hi_x);
            int iz_lo = std::max((int)std::floor(std::min(z0, z1)) - 1, SC.tile_lo_z);
            int iz_hi = std::min((int)std::floor(std::max(z0, z1)) + 1, SC.tile_hi_z);
            bool has_dx = std::abs(dx) > 1e-9;
            double inv_dx = has_dx ? 1.0 / dx : 0.0;
            bool has_dz = std::abs(dz) > 1e-9;
            double inv_dz = has_dz ? 1.0 / dz : 0.0;

            for (int ix = ix_lo; ix <= ix_hi; ix++) {
                for (int iz = iz_lo; iz <= iz_hi; iz++) {
                    double th = get_tile_height(ix, iz);
                    if (th < 0) continue;

                    // Floor cube AABB
                    double ft; Vec fn;
                    if (aabb_test(ox, oy, oz, dx, dy, dz,
                                  ix, ix + 1.0, 0, th, iz, iz + 1.0,
                                  has_dx, inv_dx, inv_dy, has_dz, inv_dz, t, ft, fn)) {
                        t = ft; n = fn;
                        m = ((ix + 100) + (iz + 100)) % 2 ? 1 : 2;
                    }

                    // Glass cube AABB
                    double gt; Vec gn;
                    if (aabb_test(ox, oy, oz, dx, dy, dz,
                                  ix - 0.1, ix + 1.1, th, th + SC.glass_edge, iz - 0.1, iz + 1.1,
                                  has_dx, inv_dx, inv_dy, has_dz, inv_dz, t, gt, gn)) {
                        t = gt; n = gn; m = 4;
                        glass_hit = {ix, iz, th};
                    }
                }
            }
        }
    }

    // Text spheres
    for (int ci = 0; ci < (int)centers.size(); ci++) {
        const auto& c = centers[ci];
        double px = ox - c.x, py = oy - c.y, pz = oz - c.z;
        double b = px*dx + py*dy + pz*dz;
        double disc = b*b - (px*px + py*py + pz*pz - SC.sph_r2);
        if (disc < 0) continue;
        double sq = std::sqrt(disc);
        double tt = -b - sq;
        if (tt < 0.01) tt = -b + sq;
        if (tt > 0.01 && tt < t) {
            t = tt;
            double inv_r = 1.0 / SC.sph_r;
            n = Vec((ox + dx*tt - c.x) * inv_r, (oy + dy*tt - c.y) * inv_r, (oz + dz*tt - c.z) * inv_r);
            m = 3;
            hit_sphere_idx = ci;
        }
    }

    // Rounded reflective cube (SDF ray march)
    {
        double px = ox - SC.rcube_center.x, py = oy - SC.rcube_center.y, pz = oz - SC.rcube_center.z;
        double b = px*dx + py*dy + pz*dz;
        double disc = b*b - (px*px + py*py + pz*pz - SC.rcube_bound_r2);
        if (disc >= 0) {
            double sq = std::sqrt(disc);
            double t_near = -b - sq;
            double t_far = -b + sq;
            if (t_far > 0.01 && t_near < t) {
                Vec local_o = SC.rcube_rot.mul(Vec(ox, oy, oz) - SC.rcube_center);
                Vec local_d = SC.rcube_rot.mul(Vec(dx, dy, dz));
                double tt = std::max(t_near, 0.01);
                for (int step = 0; step < 64; step++) {
                    Vec lp = local_o + local_d * tt;
                    double dist = sdf_round_box(lp);
                    if (dist < 0.001) {
                        if (tt < t) {
                            t = tt;
                            n = SC.rcube_rot_inv.mul(sdf_normal(lp));
                            m = 5;
                        }
                        break;
                    }
                    tt += dist;
                    if (tt > std::min(t_far, t)) break;
                }
            }
        }
    }

    // Red emissive disc
    if (std::abs(dz) > 1e-9) {
        double tt = (SC.disc_center.z - oz) / dz;
        if (tt > 0.01 && tt < t) {
            double hx = ox + dx * tt - SC.disc_center.x;
            double hy = oy + dy * tt - SC.disc_center.y;
            if (hx * hx + hy * hy <= SC.disc_r2) {
                t = tt; m = 6;
                n = SC.disc_normal;
            }
        }
    }

    return {t, m, n};
}

// Thread-local RNG
static thread_local std::mt19937 rng;

static inline double randf() {
    return std::uniform_real_distribution<double>(0.0, 1.0)(rng);
}

static Vec trace(const Vec& o, const Vec& d, int depth = 0);

static Vec sky_color(const Vec& d) {
    Vec sky = SC.sky_base;

    if (d.y > SC.sky_horizon) {
        double t = SC.ceiling_h / d.y;
        double gx = d.x * t;
        double gz = d.z * t;

        double spacing = SC.grid_spacing;
        double fx = gx / spacing;
        fx = std::min(fx - std::floor(fx), 1.0 - (fx - std::floor(fx))) * spacing;
        double fz = gz / spacing;
        fz = std::min(fz - std::floor(fz), 1.0 - (fz - std::floor(fz))) * spacing;

        double gw2 = SC.glow_width * SC.glow_width;
        double lx = std::exp(-fx * fx / gw2);
        double lz = std::exp(-fz * fz / gw2);
        double grid = std::max(lx, lz);

        double fade = 1.0 / (1.0 + t * SC.grid_fade);
        grid *= fade;

        sky = sky + SC.neon_color * (grid * SC.neon_intensity);
    }

    return sky;
}

static Vec trace(const Vec& o, const Vec& d, int depth) {
    auto [t, m, n] = test(o, d);
    if (m == 0) return sky_color(d);

    Vec hit = o + d * t;

    // Emissive disc
    if (m == 6) return SC.disc_emissive;

    Vec light1 = Vec(SC.light1_pos.x + randf(), SC.light1_pos.y + randf(), SC.light1_pos.z).norm();
    Vec light2 = Vec(SC.light2_pos.x + randf(), SC.light2_pos.y + randf(), SC.light2_pos.z).norm();
    double angle = randf() * 2.0 * M_PI;
    double rr = std::sqrt(randf()) * SC.disc_r;
    Vec disc_pt(SC.disc_center.x + rr * std::cos(angle), SC.disc_center.y + rr * std::sin(angle), SC.disc_center.z);
    Vec light3 = (disc_pt - hit).norm();
    Vec L1 = SC.light1_color, L2 = SC.light2_color, L3 = SC.light3_color;

    // Shadow
    auto [_, lm1, __] = test(hit + n * 0.01, light1, true);
    auto [_2, lm2, _3] = test(hit + n * 0.01, light2, true);
    auto [_4, lm3, _5] = test(hit + n * 0.01, light3, true);
    double s1 = lm1 ? SC.shadow_ambient : 1.0;
    double s2 = lm2 ? SC.shadow_ambient : 1.0;
    double s3 = (lm3 && lm3 != 6) ? SC.shadow_ambient : 1.0;
    double d1 = std::max(0.0, light1.dot(n)) * s1;
    double d2 = std::max(0.0, light2.dot(n)) * s2;
    double d3 = std::max(0.0, light3.dot(n)) * s3;

    if (m == 4 && depth < SC.glass_max_depth) {
        // Glass cube — refraction + reflection
        double cos_i = -(d.dot(n));
        Vec nn = n;
        if (cos_i < 0) { cos_i = -cos_i; nn = n * -1; }
        double fresnel = SC.glass_r0 + (1 - SC.glass_r0) * std::pow(1 - cos_i, 5);
        Vec r = d - nn * (2 * d.dot(nn));

        Vec transmitted(0, 0, 0);
        double eta = 1.0 / SC.glass_ior;
        double k = 1.0 - eta * eta * (1.0 - cos_i * cos_i);
        if (k > 0) {
            Vec refr = d * eta + nn * (eta * cos_i - std::sqrt(k));
            auto gi = glass_hit;
            double gx0 = gi.ix - 0.1, gx1 = gi.ix + 1.1;
            double gy0 = gi.th, gy1 = gi.th + SC.glass_edge;
            double gz0 = gi.iz - 0.1, gz1 = gi.iz + 1.1;
            Vec entry_pt = hit + refr * 0.02;
            double rx = refr.x, ry = refr.y, rz = refr.z;
            double exit_t = 1e9;
            Vec en(0, 1, 0);
            if (std::abs(rx) > 1e-9) {
                double te = ((rx > 0 ? gx1 : gx0) - entry_pt.x) / rx;
                if (te > 0 && te < exit_t) { exit_t = te; en = Vec(rx > 0 ? 1 : -1, 0, 0); }
            }
            if (std::abs(ry) > 1e-9) {
                double te = ((ry > 0 ? gy1 : gy0) - entry_pt.y) / ry;
                if (te > 0 && te < exit_t) { exit_t = te; en = Vec(0, ry > 0 ? 1 : -1, 0); }
            }
            if (std::abs(rz) > 1e-9) {
                double te = ((rz > 0 ? gz1 : gz0) - entry_pt.z) / rz;
                if (te > 0 && te < exit_t) { exit_t = te; en = Vec(0, 0, rz > 0 ? 1 : -1); }
            }
            Vec exit_pt = entry_pt + refr * exit_t;
            double cos_o = -(refr.dot(en));
            if (cos_o < 0) { cos_o = -cos_o; en = en * -1; }
            double eta2 = SC.glass_ior;
            double k2 = 1.0 - eta2 * eta2 * (1.0 - cos_o * cos_o);
            if (k2 > 0) {
                Vec exit_refr = refr * eta2 + en * (eta2 * cos_o - std::sqrt(k2));
                transmitted = trace(exit_pt + exit_refr * 0.02, exit_refr, depth + 1);
            } else {
                Vec refl_int = refr - en * (2 * refr.dot(en));
                transmitted = trace(exit_pt + refl_int * 0.02, refl_int, depth + 1);
            }
            double pl = exit_t * SC.glass_tint_scale;
            transmitted = Vec(transmitted.x * (1 - pl*SC.glass_tint_rgb.x),
                              transmitted.y * (1 - pl*SC.glass_tint_rgb.y),
                              transmitted.z * (1 - pl*SC.glass_tint_rgb.z));
        } else {
            fresnel = 1.0;
        }
        Vec reflected = fresnel > 0.02 ? trace(hit + nn * 0.02, r, depth + 1) : Vec();
        Vec base = reflected * fresnel + transmitted * (1.0 - fresnel);
        double sp1 = std::pow(std::max(0.0, light1.dot(r)), SC.glass_spec_power) * s1;
        double sp2 = std::pow(std::max(0.0, light2.dot(r)), SC.glass_spec_power) * s2;
        return base + L1 * (sp1 * SC.glass_spec_intensity) + L2 * (sp2 * SC.glass_spec_intensity);
    }
    if (m == 4) {
        return depth < SC.glass_passthrough_depth ? trace(hit + d * 0.02, d, depth + 1) : Vec();
    }
    if (m == 5) {
        // Rounded reflective cube — chrome
        Vec r = d - n * (2 * d.dot(n));
        Vec color = SC.chrome_color;
        Vec lit = L1 * d1 + L2 * d2 + L3 * d3;
        Vec base(color.x * (lit.x * SC.chrome_diffuse + SC.chrome_ambient),
                 color.y * (lit.y * SC.chrome_diffuse + SC.chrome_ambient),
                 color.z * (lit.z * SC.chrome_diffuse + SC.chrome_ambient));
        if (depth < SC.chrome_max_depth) {
            Vec ref = trace(hit + n * 0.01, r, depth + 1);
            base = base * SC.chrome_refl_base + ref * SC.chrome_refl_mirror;
        }
        double sp1 = std::pow(std::max(0.0, light1.dot(r)), SC.chrome_spec_power) * s1;
        double sp2 = std::pow(std::max(0.0, light2.dot(r)), SC.chrome_spec_power) * s2;
        double sp3 = std::pow(std::max(0.0, light3.dot(r)), SC.chrome_spec_power) * s3;
        return base + L1 * (sp1 * SC.chrome_spec_intensity) + L2 * (sp2 * SC.chrome_spec_intensity) + L3 * (sp3 * SC.chrome_spec_intensity);
    }
    if (m == 3) {
        // Text sphere — silver reflective with scratches
        double scratch_amount = 0.0;
        {
            Vec local = hit - Vec(centers[hit_sphere_idx].x, centers[hit_sphere_idx].y, centers[hit_sphere_idx].z);
            uint32_t seed = (uint32_t)(hit_sphere_idx * 2654435761u);
            double sa = (seed & 0xFFFF) / 65535.0 * M_PI;
            double sp = ((seed >> 16) & 0xFFFF) / 65535.0 * 10.0;
            double theta = std::atan2(local.z, local.x);
            double phi = std::asin(std::clamp(local.y / SC.sph_r, -1.0, 1.0));
            double u = theta * SC.scratch_uv, v = phi * SC.scratch_uv;
            // Primary scratches
            double rot_u = u * std::cos(sa) - v * std::sin(sa);
            double scratch = std::sin(rot_u * SC.s1_freq + sp) * 0.5 + 0.5;
            scratch = std::pow(scratch, SC.s1_power) * SC.s1_intensity;
            // Secondary scratches
            double sa2 = sa + SC.s2_angle;
            double rot_u2 = u * std::cos(sa2) - v * std::sin(sa2);
            double scratch2 = std::sin(rot_u2 * SC.s2_freq + sp * SC.s2_phase) * 0.5 + 0.5;
            scratch2 = std::pow(scratch2, SC.s2_power) * SC.s2_intensity;
            // Tertiary scratches
            double sa3 = sa + SC.s3_angle;
            double rot_u3 = u * std::cos(sa3) - v * std::sin(sa3);
            double scratch3 = std::sin(rot_u3 * SC.s3_freq + sp * SC.s3_phase) * 0.5 + 0.5;
            scratch3 = std::pow(scratch3, SC.s3_power) * SC.s3_intensity;
            double total = scratch + scratch2 + scratch3;
            scratch_amount = total;
            Vec tang1 = n.cross(Vec(0, 1, 0)).norm();
            if (tang1.len() < 0.1) tang1 = n.cross(Vec(1, 0, 0)).norm();
            Vec tang2 = n.cross(tang1);
            double pa = std::cos(sa) * total, pb = std::sin(sa) * total;
            n = (n + tang1 * pa + tang2 * pb).norm();
        }
        Vec r = d - n * (2 * d.dot(n));
        double dis = std::min(scratch_amount * SC.scratch_discolor, 1.0);
        Vec color = SC.silver_color * (1.0 - dis) + SC.scratched_color * dis;
        Vec lit = L1 * d1 + L2 * d2 + L3 * d3;
        Vec base(color.x * (lit.x * SC.silver_diffuse + SC.silver_ambient),
                 color.y * (lit.y * SC.silver_diffuse + SC.silver_ambient),
                 color.z * (lit.z * SC.silver_diffuse + SC.silver_ambient));
        if (depth < SC.silver_max_depth) {
            Vec ref = trace(hit + n * 0.01, r, depth + 1);
            base = base * SC.silver_refl_base + ref * SC.silver_refl_mirror;
        }
        double sp1 = std::pow(std::max(0.0, light1.dot(r)), SC.silver_spec_power) * s1;
        double sp2 = std::pow(std::max(0.0, light2.dot(r)), SC.silver_spec_power) * s2;
        double sp3 = std::pow(std::max(0.0, light3.dot(r)), SC.silver_spec_power) * s3;
        return base + L1 * (sp1 * SC.silver_spec_intensity) + L2 * (sp2 * SC.silver_spec_intensity) + L3 * (sp3 * SC.silver_spec_intensity);
    }
    // Floor tiles
    Vec color = (m == 1) ? SC.floor_white : SC.floor_red;
    Vec r = d - n * (2 * d.dot(n));
    Vec lit = L1 * d1 + L2 * d2 + L3 * d3;
    Vec base(color.x * (lit.x * SC.floor_diffuse + SC.floor_ambient),
             color.y * (lit.y * SC.floor_diffuse + SC.floor_ambient),
             color.z * (lit.z * SC.floor_diffuse + SC.floor_ambient));
    if (depth < SC.floor_max_depth) {
        Vec ref = trace(hit + n * 0.01, r, depth + 1);
        base = base * SC.floor_refl_base + ref * SC.floor_refl_mirror;
    }
    double sp1 = std::pow(std::max(0.0, light1.dot(r)), SC.floor_spec_power) * s1;
    double sp2 = std::pow(std::max(0.0, light2.dot(r)), SC.floor_spec_power) * s2;
    double sp3 = std::pow(std::max(0.0, light3.dot(r)), SC.floor_spec_power) * s3;
    return base + L1 * (sp1 * SC.floor_spec_intensity) + L2 * (sp2 * SC.floor_spec_intensity) + L3 * (sp3 * SC.floor_spec_intensity);
}

// 3D noise for volumetric fog
static inline double hash3d(int x, int y, int z) {
    uint32_t h = (uint32_t)(x * 374761393u + y * 668265263u + z * 1274126177u);
    h = (h ^ (h >> 13)) * 1274126177u;
    return (h & 0xFFFF) / 65535.0;
}

static double noise3d(double x, double y, double z) {
    int ix = (int)std::floor(x), iy = (int)std::floor(y), iz = (int)std::floor(z);
    double fx = x - ix, fy = y - iy, fz = z - iz;
    double sx = fx * fx * (3 - 2 * fx);
    double sy = fy * fy * (3 - 2 * fy);
    double sz = fz * fz * (3 - 2 * fz);
    double c000 = hash3d(ix, iy, iz),     c100 = hash3d(ix+1, iy, iz);
    double c010 = hash3d(ix, iy+1, iz),   c110 = hash3d(ix+1, iy+1, iz);
    double c001 = hash3d(ix, iy, iz+1),   c101 = hash3d(ix+1, iy, iz+1);
    double c011 = hash3d(ix, iy+1, iz+1), c111 = hash3d(ix+1, iy+1, iz+1);
    double c00 = c000 + (c100 - c000) * sx, c01 = c001 + (c101 - c001) * sx;
    double c10 = c010 + (c110 - c010) * sx, c11 = c011 + (c111 - c011) * sx;
    double c0 = c00 + (c10 - c00) * sy,     c1 = c01 + (c11 - c01) * sy;
    return c0 + (c1 - c0) * sz;
}

static double fbm3d(double x, double y, double z) {
    double val = 0, amp = 0.5;
    for (int i = 0; i < 4; i++) {
        val += noise3d(x, y, z) * amp;
        x *= 2; y *= 2; z *= 2;
        amp *= 0.5;
    }
    return val;
}

// Volumetric fog
static Vec apply_fog(const Vec& o, const Vec& d, double t_max, const Vec& surface_color) {
    if (t_max <= 0.01) return surface_color;
    double step_size = t_max / SC.fog_steps;
    double transmittance = 1.0;
    Vec inscatter(0, 0, 0);

    for (int i = 0; i < SC.fog_steps; i++) {
        double t = (i + randf()) * step_size;
        Vec pos = o + d * t;

        Vec l1 = Vec(SC.light1_pos.x + randf(), SC.light1_pos.y + randf(), SC.light1_pos.z).norm();
        Vec l2 = Vec(SC.light2_pos.x + randf(), SC.light2_pos.y + randf(), SC.light2_pos.z).norm();
        double ang = randf() * 2.0 * M_PI;
        double rr = std::sqrt(randf()) * SC.disc_r;
        Vec disc_pt(SC.disc_center.x + rr * std::cos(ang),
                    SC.disc_center.y + rr * std::sin(ang),
                    SC.disc_center.z);
        Vec l3 = (disc_pt - pos).norm();

        Vec L1 = SC.light1_color, L2 = SC.light2_color, L3 = SC.light3_color;

        auto [t1, m1, n1] = test(pos, l1, false);
        auto [t2, m2, n2] = test(pos, l2, false);
        auto [t3, m3, n3] = test(pos, l3, false);

        double s1 = m1 ? 0.0 : 1.0;
        double s2 = m2 ? 0.0 : 1.0;
        double s3 = (m3 && m3 != 6) ? 0.0 : 1.0;

        Vec light_contrib = L1 * s1 + L2 * s2 + L3 * s3;

        double nx = pos.x * SC.fog_nx + pos.z * SC.fog_nxz;
        double ny = pos.y * SC.fog_ny;
        double nz = pos.z * SC.fog_nz;
        double local_density = SC.fog_density * fbm3d(nx, ny, nz);

        inscatter = inscatter + light_contrib * (local_density * transmittance * step_size);
        transmittance *= std::exp(-local_density * step_size);
    }

    return surface_color * transmittance + inscatter;
}

static inline int clamp(double v) {
    int i = (int)(v * 255.0);
    return i < 0 ? 0 : (i > 255 ? 255 : i);
}

int main(int argc, char* argv[]) {
    // Load scene from YAML
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <scene.yaml>\n", argv[0]);
        return 1;
    }
    const char* scene_file = argv[1];
    YamlConfig yaml;
    if (!yaml.load(scene_file)) {
        fprintf(stderr, "Error: could not load %s\n", scene_file);
        return 1;
    }
    std::string version = yaml.get_string("version");
    if (version != "nanore@1") {
        fprintf(stderr, "Error: unsupported scene version '%s' (expected nanore@1)\n", version.c_str());
        return 1;
    }
    SC.load_from_yaml(yaml);
    printf("Loaded scene: %s\n", scene_file);

    init_tiles();
    init_centers();

    int W = SC.width, H = SC.height;
    printf("nanore ray tracer\n");
    printf("  Resolution: %dx%d, %d samples/pixel\n", W, H, SC.samples);
    printf("  %zu spheres (r=%.1f), %d displaced cubes + glass\n", centers.size(), SC.sph_r, SC.tile_nx * SC.tile_nz);
    printf("  Rounded reflective cube (edge=%.0f, round=%.1f) behind text\n", SC.rcube_half * 2, SC.rcube_round);
    printf("  Blue + yellow dual lighting, shiny floor\n");
    fflush(stdout);

    Vec fwd = (SC.target - SC.cam).norm();
    Vec right = fwd.cross(Vec(0, 1, 0)).norm();
    Vec up = right.cross(fwd);

    auto* pixels = new uint8_t[W * H * 3];
    auto* hdr = new float[W * H * 3];

    int num_threads = std::max(1, (int)std::thread::hardware_concurrency());
    printf("  Using %d threads\n", num_threads);
    fflush(stdout);

    std::atomic<int> next_row(0);
    std::atomic<int> rows_done(0);
    auto start = clock();

    auto worker = [&]() {
        rng.seed(std::random_device{}());
        while (true) {
            int y = next_row.fetch_add(1);
            if (y >= H) break;
            for (int x = 0; x < W; x++) {
                Vec color;
                for (int s = 0; s < SC.samples; s++) {
                    double dx = (x + randf() - W / 2.0) / H;
                    double dy = (H / 2.0 - y + randf()) / H;
                    Vec dir = (right * dx + up * dy + fwd).norm();
                    auto [ht, hm, hn] = test(SC.cam, dir);
                    Vec c = trace(SC.cam, dir);
                    double fog_dist = (hm == 0) ? SC.fog_max_dist : ht;
                    c = apply_fog(SC.cam, dir, fog_dist, c);
                    color = color + c;
                }
                color = color * (1.0 / SC.samples);
                int idx = (y * W + x) * 3;
                hdr[idx]   = (float)color.x;
                hdr[idx+1] = (float)color.y;
                hdr[idx+2] = (float)color.z;
            }
            int done = rows_done.fetch_add(1) + 1;
            if (done % 32 == 0) {
                double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
                double eta_s = elapsed / done * (H - done);
                printf("  Rendering... %d%% (ETA: %.0fs)\n", done * 100 / H, eta_s);
                fflush(stdout);
            }
        }
    };

    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; i++) threads.emplace_back(worker);
    for (auto& t : threads) t.join();

    double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("  Render time: %.1fs CPU (%.1fs wall)\n", elapsed, elapsed / num_threads);
    fflush(stdout);

    // Bloom post-processing
    printf("  Applying glow filter (r=%d)...\n", SC.bloom_r); fflush(stdout);
    double bloom_sigma = SC.bloom_r / 3.0;

    std::vector<double> weights(SC.bloom_r + 1);
    double wsum = 0;
    for (int i = 0; i <= SC.bloom_r; i++) {
        weights[i] = std::exp(-0.5 * (i * i) / (bloom_sigma * bloom_sigma));
        wsum += (i == 0) ? weights[i] : 2 * weights[i];
    }
    for (int i = 0; i <= SC.bloom_r; i++) weights[i] /= wsum;

    auto* bloom = new float[W * H * 3]();
    for (int i = 0; i < W * H * 3; i++) {
        float v = hdr[i];
        bloom[i] = (v > SC.bloom_threshold) ? (v - (float)SC.bloom_threshold) : 0.0f;
    }

    // Horizontal blur
    auto* temp = new float[W * H * 3]();
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            double r = 0, g = 0, b = 0;
            for (int k = -SC.bloom_r; k <= SC.bloom_r; k++) {
                int sx = std::clamp(x + k, 0, W - 1);
                double w = weights[std::abs(k)];
                int si = (y * W + sx) * 3;
                r += bloom[si] * w;
                g += bloom[si + 1] * w;
                b += bloom[si + 2] * w;
            }
            int di = (y * W + x) * 3;
            temp[di] = (float)r; temp[di+1] = (float)g; temp[di+2] = (float)b;
        }
    }

    // Vertical blur
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            double r = 0, g = 0, b = 0;
            for (int k = -SC.bloom_r; k <= SC.bloom_r; k++) {
                int sy = std::clamp(y + k, 0, H - 1);
                double w = weights[std::abs(k)];
                int si = (sy * W + x) * 3;
                r += temp[si] * w;
                g += temp[si + 1] * w;
                b += temp[si + 2] * w;
            }
            int di = (y * W + x) * 3;
            bloom[di] = (float)r; bloom[di+1] = (float)g; bloom[di+2] = (float)b;
        }
    }

    // Combine + gamma
    for (int i = 0; i < W * H * 3; i++) {
        double v = hdr[i] + bloom[i] * SC.bloom_intensity;
        pixels[i] = clamp(std::pow(v, SC.gamma));
    }
    delete[] bloom;
    delete[] temp;
    printf("  Glow applied.\n"); fflush(stdout);

    // Build output path: renders/render-<timestamp>.bmp
    mkdir("renders", 0755);
    time_t now = time(nullptr);
    struct tm* t = localtime(&now);
    char outpath[256];
    strftime(outpath, sizeof(outpath), "renders/render-%Y%m%d-%H%M%S.bmp", t);

    // Save BMP
    FILE* f = fopen(outpath, "wb");
    int row_size = W * 3;
    int pad = (4 - row_size % 4) % 4;
    int img_size = (row_size + pad) * H;
    uint8_t bmp_header[54] = {};
    bmp_header[0] = 'B'; bmp_header[1] = 'M';
    *(uint32_t*)(bmp_header + 2) = 54 + img_size;
    *(uint32_t*)(bmp_header + 10) = 54;
    *(uint32_t*)(bmp_header + 14) = 40;
    *(int32_t*)(bmp_header + 18) = W;
    *(int32_t*)(bmp_header + 22) = H;
    *(uint16_t*)(bmp_header + 26) = 1;
    *(uint16_t*)(bmp_header + 28) = 24;
    *(uint32_t*)(bmp_header + 34) = img_size;
    *(int32_t*)(bmp_header + 38) = 2835;
    *(int32_t*)(bmp_header + 42) = 2835;
    fwrite(bmp_header, 1, 54, f);
    uint8_t padding[4] = {};
    for (int y = H - 1; y >= 0; y--) {
        for (int x = 0; x < W; x++) {
            int idx = (y * W + x) * 3;
            uint8_t bgr[3] = {pixels[idx+2], pixels[idx+1], pixels[idx]};
            fwrite(bgr, 1, 3, f);
        }
        if (pad) fwrite(padding, 1, pad, f);
    }
    fclose(f);
    printf("  Saved to %s\n", outpath);

    delete[] hdr;
    delete[] pixels;
    return 0;
}
