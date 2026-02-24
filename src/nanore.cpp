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

    int get_int(const std::string& key, int def) const {
        auto it = data.find(key);
        if (it == data.end()) return def;
        return (int)std::strtol(it->second.c_str(), nullptr, 10);
    }

    double get_double(const std::string& key, double def) const {
        auto it = data.find(key);
        if (it == data.end()) return def;
        return std::strtod(it->second.c_str(), nullptr);
    }

    Vec get_vec(const std::string& key, Vec def) const {
        auto it = data.find(key);
        if (it == data.end()) return def;
        double x, y, z;
        if (sscanf(it->second.c_str(), "[%lf, %lf, %lf]", &x, &y, &z) == 3)
            return Vec(x, y, z);
        if (sscanf(it->second.c_str(), "[%lf,%lf,%lf]", &x, &y, &z) == 3)
            return Vec(x, y, z);
        return def;
    }

    std::string get_string(const std::string& key, const std::string& def = "") const {
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

    std::vector<std::string> get_keys_under(const std::string& prefix) const {
        std::vector<std::string> result;
        std::string pfx = prefix + ".";
        std::string last;
        for (const auto& kv : data) {
            if (kv.first.compare(0, pfx.size(), pfx) != 0) continue;
            size_t dot = kv.first.find('.', pfx.size());
            std::string child = (dot != std::string::npos)
                ? kv.first.substr(pfx.size(), dot - pfx.size())
                : kv.first.substr(pfx.size());
            if (child != last) {
                result.push_back(child);
                last = child;
            }
        }
        return result;
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
    int rcube_count = 1;
    double rcube_half_min = 0.8, rcube_half_max = 3.2;
    double rcube_spread = 8.0;
    double rcube_round_ratio = 0.2;
    int rcube_seed = 42;

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

    // Reflective inner spheres (inside glass spheres)
    double rsph_radius_mult = 0.625; // half the glass sphere diameter = gsph * 0.5

    // Glass spheres (around text spheres)
    double gsph_radius_mult = 1.25;
    double gsph_ior = 1.5, gsph_r0 = 0.04;
    double gsph_tint_scale = 0.1;
    Vec gsph_tint_rgb{0.02, 0.01, 0.005};
    double gsph_spec_power = 128, gsph_spec_intensity = 0.3;
    int gsph_max_depth = 4, gsph_passthrough_depth = 6;

    // Emissive tori (around text spheres)
    double torus_major_mult = 1.0;
    double torus_minor = 0.07;
    double torus_emissive_intensity = 2.0;

    // Glass
    double glass_edge = 1.2, glass_ior = 1.5, glass_r0 = 0.04;
    double glass_tint_scale = 0.15;
    Vec glass_tint_rgb{0.05, 0.02, 0.01};
    double glass_spec_power = 96, glass_spec_intensity = 0.3;
    int glass_max_depth = 4, glass_passthrough_depth = 6;
    int glass_cluster_min = 0, glass_cluster_max = 0;
    double glass_size_min = 0.2, glass_size_max = 0.6;
    double glass_spread = 0.5;

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

    // Font-rendered text blocks
    struct TextBlock {
        std::string content;
        double x, y, z;
        double col_spacing = 1.0;
        double row_spacing = 1.0;
    };
    std::vector<TextBlock> text_blocks;

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
        glass_max_y = tile_max_y + glass_edge; // overridden by init_glass if clustering
        tile_nx = tile_hi_x - tile_lo_x + 1;
        tile_nz = tile_hi_z - tile_lo_z + 1;
        rcube_rot = make_rotation(rcube_rotation.x, rcube_rotation.y, rcube_rotation.z);
        rcube_rot_inv = rcube_rot.transpose();
    }

    void load_from_yaml(const YamlConfig& cfg) {
        // Rendering
        width = cfg.get_int("rendering.width", width);
        height = cfg.get_int("rendering.height", height);
        samples = cfg.get_int("rendering.samples", samples);
        fog_max_dist = cfg.get_double("rendering.fog_max_dist", fog_max_dist);

        // Camera
        cam = cfg.get_vec("camera.position", cam);
        target = cfg.get_vec("camera.target", target);

        // Objects
        int light_idx = 0;
        text_blocks.clear();
        for (const auto& name : cfg.get_keys_under("objects")) {
            std::string p = "objects." + name + ".";
            std::string type = cfg.get_string(p + "type");

            if (type == "bitmap_text") {
                bitmap_rows.clear();
                for (int i = 0; ; i++) {
                    std::string row = cfg.get_string(p + "bitmap." + std::to_string(i));
                    if (row.empty()) break;
                    bitmap_rows.push_back(row);
                }
                sph_r = cfg.get_double(p + "sphere_radius", sph_r);
                text_x_offset = cfg.get_double(p + "x_offset", text_x_offset);
                text_y_base = cfg.get_double(p + "y_base", text_y_base);
                text_y_spacing = cfg.get_double(p + "y_spacing", text_y_spacing);
                text_z = cfg.get_double(p + "z_depth", text_z);
                slash_x_base = cfg.get_double(p + "slash_x_base", slash_x_base);
                slash_y_base = cfg.get_double(p + "slash_y_base", slash_y_base);
                slash_count = cfg.get_int(p + "slash_count", slash_count);
            } else if (type == "rounded_cube") {
                rcube_center = cfg.get_vec(p + "center", rcube_center);
                rcube_half = cfg.get_double(p + "half_edge", rcube_half);
                rcube_round = cfg.get_double(p + "rounding", rcube_round);
                rcube_rotation = cfg.get_vec(p + "rotation", rcube_rotation);
                rcube_count = cfg.get_int(p + "count", rcube_count);
                rcube_half_min = cfg.get_double(p + "half_min", rcube_half_min);
                rcube_half_max = cfg.get_double(p + "half_max", rcube_half_max);
                rcube_spread = cfg.get_double(p + "spread", rcube_spread);
                rcube_round_ratio = cfg.get_double(p + "round_ratio", rcube_round_ratio);
                rcube_seed = cfg.get_int(p + "seed", rcube_seed);
            } else if (type == "disc") {
                disc_center = cfg.get_vec(p + "center", disc_center);
                disc_normal = cfg.get_vec(p + "normal", disc_normal);
                disc_r = cfg.get_double(p + "radius", disc_r);
                disc_emissive = cfg.get_vec(p + "emissive_color", disc_emissive);
                light3_color = cfg.get_vec(p + "light_color", light3_color);
            } else if (type == "floor") {
                Vec bx = cfg.get_vec(p + "bounds_x", Vec(tile_lo_x, tile_hi_x, 0));
                tile_lo_x = (int)bx.x; tile_hi_x = (int)bx.y;
                Vec bz = cfg.get_vec(p + "bounds_z", Vec(tile_lo_z, tile_hi_z, 0));
                tile_lo_z = (int)bz.x; tile_hi_z = (int)bz.y;
                tile_max_y = cfg.get_double(p + "max_height", tile_max_y);
                floor_white = cfg.get_vec(p + "white_color", floor_white);
                floor_red = cfg.get_vec(p + "red_color", floor_red);
                floor_diffuse = cfg.get_double(p + "diffuse", floor_diffuse);
                floor_ambient = cfg.get_double(p + "ambient", floor_ambient);
                floor_refl_base = cfg.get_double(p + "reflect_base", floor_refl_base);
                floor_refl_mirror = cfg.get_double(p + "reflect_mirror", floor_refl_mirror);
                floor_max_depth = cfg.get_int(p + "max_depth", floor_max_depth);
                floor_spec_power = cfg.get_double(p + "spec_power", floor_spec_power);
                floor_spec_intensity = cfg.get_double(p + "spec_intensity", floor_spec_intensity);
            } else if (type == "glass") {
                glass_edge = cfg.get_double(p + "edge", glass_edge);
                glass_ior = cfg.get_double(p + "ior", glass_ior);
                glass_r0 = cfg.get_double(p + "fresnel_r0", glass_r0);
                glass_tint_scale = cfg.get_double(p + "tint_scale", glass_tint_scale);
                glass_tint_rgb = cfg.get_vec(p + "tint_rgb", glass_tint_rgb);
                glass_spec_power = cfg.get_double(p + "spec_power", glass_spec_power);
                glass_spec_intensity = cfg.get_double(p + "spec_intensity", glass_spec_intensity);
                glass_max_depth = cfg.get_int(p + "max_depth", glass_max_depth);
                glass_passthrough_depth = cfg.get_int(p + "passthrough_depth", glass_passthrough_depth);
                glass_cluster_min = cfg.get_int(p + "cluster_min", glass_cluster_min);
                glass_cluster_max = cfg.get_int(p + "cluster_max", glass_cluster_max);
                glass_size_min = cfg.get_double(p + "size_min", glass_size_min);
                glass_size_max = cfg.get_double(p + "size_max", glass_size_max);
                glass_spread = cfg.get_double(p + "spread", glass_spread);
            } else if (type == "light") {
                Vec pos = cfg.get_vec(p + "position", Vec());
                Vec col = cfg.get_vec(p + "color", Vec());
                if (light_idx == 0) { light1_pos = pos; light1_color = col; }
                else if (light_idx == 1) { light2_pos = pos; light2_color = col; }
                light_idx++;
            } else if (type == "text") {
                TextBlock tb;
                tb.content = cfg.get_string(p + "string");
                Vec pos = cfg.get_vec(p + "position", Vec(0, 0, 0));
                tb.x = pos.x; tb.y = pos.y; tb.z = pos.z;
                tb.col_spacing = cfg.get_double(p + "col_spacing", 1.0);
                tb.row_spacing = cfg.get_double(p + "row_spacing", 1.0);
                text_blocks.push_back(tb);
            } else if (type == "reflective_sphere") {
                rsph_radius_mult = cfg.get_double(p + "radius_mult", rsph_radius_mult);
            } else if (type == "glass_sphere") {
                gsph_radius_mult = cfg.get_double(p + "radius_mult", gsph_radius_mult);
                gsph_ior = cfg.get_double(p + "ior", gsph_ior);
                gsph_r0 = cfg.get_double(p + "fresnel_r0", gsph_r0);
                gsph_tint_scale = cfg.get_double(p + "tint_scale", gsph_tint_scale);
                gsph_tint_rgb = cfg.get_vec(p + "tint_rgb", gsph_tint_rgb);
                gsph_spec_power = cfg.get_double(p + "spec_power", gsph_spec_power);
                gsph_spec_intensity = cfg.get_double(p + "spec_intensity", gsph_spec_intensity);
                gsph_max_depth = cfg.get_int(p + "max_depth", gsph_max_depth);
                gsph_passthrough_depth = cfg.get_int(p + "passthrough_depth", gsph_passthrough_depth);
            } else if (type == "emissive_torus") {
                torus_major_mult = cfg.get_double(p + "major_mult", torus_major_mult);
                torus_minor = cfg.get_double(p + "minor", torus_minor);
                torus_emissive_intensity = cfg.get_double(p + "intensity", torus_emissive_intensity);
            }
        }

        // Materials
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

        scratch_uv = cfg.get_double("materials.scratches.uv_scale", scratch_uv);
        s1_freq = cfg.get_double("materials.scratches.layer1.frequency", s1_freq);
        s1_power = cfg.get_double("materials.scratches.layer1.power", s1_power);
        s1_intensity = cfg.get_double("materials.scratches.layer1.intensity", s1_intensity);
        s2_angle = cfg.get_double("materials.scratches.layer2.angle_offset", s2_angle);
        s2_freq = cfg.get_double("materials.scratches.layer2.frequency", s2_freq);
        s2_phase = cfg.get_double("materials.scratches.layer2.phase_mult", s2_phase);
        s2_power = cfg.get_double("materials.scratches.layer2.power", s2_power);
        s2_intensity = cfg.get_double("materials.scratches.layer2.intensity", s2_intensity);
        s3_angle = cfg.get_double("materials.scratches.layer3.angle_offset", s3_angle);
        s3_freq = cfg.get_double("materials.scratches.layer3.frequency", s3_freq);
        s3_phase = cfg.get_double("materials.scratches.layer3.phase_mult", s3_phase);
        s3_power = cfg.get_double("materials.scratches.layer3.power", s3_power);
        s3_intensity = cfg.get_double("materials.scratches.layer3.intensity", s3_intensity);

        // Environment
        shadow_ambient = cfg.get_double("environment.shadow_ambient", shadow_ambient);
        sky_base = cfg.get_vec("environment.sky.base_color", sky_base);
        sky_horizon = cfg.get_double("environment.sky.horizon_threshold", sky_horizon);
        ceiling_h = cfg.get_double("environment.sky.ceiling_height", ceiling_h);
        grid_spacing = cfg.get_double("environment.sky.grid_spacing", grid_spacing);
        glow_width = cfg.get_double("environment.sky.glow_width", glow_width);
        grid_fade = cfg.get_double("environment.sky.fade_rate", grid_fade);
        neon_color = cfg.get_vec("environment.sky.neon_color", neon_color);
        neon_intensity = cfg.get_double("environment.sky.neon_intensity", neon_intensity);
        fog_density = cfg.get_double("environment.fog.density", fog_density);
        fog_steps = cfg.get_int("environment.fog.steps", fog_steps);
        fog_nx = cfg.get_double("environment.fog.noise_x", fog_nx);
        fog_nxz = cfg.get_double("environment.fog.noise_xz", fog_nxz);
        fog_ny = cfg.get_double("environment.fog.noise_y", fog_ny);
        fog_nz = cfg.get_double("environment.fog.noise_z", fog_nz);

        // Post-processing
        bloom_r = cfg.get_int("post_processing.bloom.radius", bloom_r);
        bloom_threshold = cfg.get_double("post_processing.bloom.threshold", bloom_threshold);
        bloom_intensity = cfg.get_double("post_processing.bloom.intensity", bloom_intensity);
        gamma = cfg.get_double("post_processing.bloom.gamma", gamma);

        compute_derived();
    }
};

static SceneConfig SC;

// Sphere centers (text + slash)
struct SphCenter { double x, y, z; };
static std::vector<SphCenter> centers;
static std::vector<Vec> torus_colors;

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

// Built-in 5x6 pixel font — each glyph is 30 chars (6 rows x 5 cols), '@'=filled
struct Glyph { char ch; char p[31]; };
static const Glyph FONT[] = {
    // Uppercase
    {'A', ".@@@." "@...@" "@@@@@" "@...@" "@...@" "@...@"},
    {'B', "@@@@." "@...@" "@@@@." "@...@" "@...@" "@@@@."},
    {'C', ".@@@@" "@...." "@...." "@...." "@...." ".@@@@"},
    {'D', "@@@@." "@...@" "@...@" "@...@" "@...@" "@@@@."},
    {'E', "@@@@@" "@...." "@@@@." "@...." "@...." "@@@@@"},
    {'F', "@@@@@" "@...." "@@@@." "@...." "@...." "@...."},
    {'G', ".@@@@" "@...." "@..@@" "@...@" "@...@" ".@@@@"},
    {'H', "@...@" "@...@" "@@@@@" "@...@" "@...@" "@...@"},
    {'I', "@@@@@" "..@.." "..@.." "..@.." "..@.." "@@@@@"},
    {'J', ".@@@@" "...@." "...@." "...@." "@..@." ".@@.."},
    {'K', "@...@" "@..@." "@@@.." "@..@." "@...@" "@...@"},
    {'L', "@...." "@...." "@...." "@...." "@...." "@@@@@"},
    {'M', "@...@" "@@.@@" "@.@.@" "@...@" "@...@" "@...@"},
    {'N', "@...@" "@@..@" "@.@.@" "@..@@" "@...@" "@...@"},
    {'O', ".@@@." "@...@" "@...@" "@...@" "@...@" ".@@@."},
    {'P', "@@@@." "@...@" "@@@@." "@...." "@...." "@...."},
    {'Q', ".@@@." "@...@" "@...@" "@.@.@" "@..@@" ".@@.@"},
    {'R', "@@@@." "@...@" "@@@@." "@..@." "@...@" "@...@"},
    {'S', ".@@@@" "@...." ".@@@." "....@" "....@" "@@@@."},
    {'T', "@@@@@" "..@.." "..@.." "..@.." "..@.." "..@.."},
    {'U', "@...@" "@...@" "@...@" "@...@" "@...@" ".@@@."},
    {'V', "@...@" "@...@" "@...@" ".@.@." ".@.@." "..@.."},
    {'W', "@...@" "@...@" "@.@.@" "@.@.@" "@@.@@" "@...@"},
    {'X', "@...@" ".@.@." "..@.." ".@.@." "@...@" "@...@"},
    {'Y', "@...@" ".@.@." "..@.." "..@.." "..@.." "..@.."},
    {'Z', "@@@@@" "...@." "..@.." ".@..." "@...." "@@@@@"},
    // Lowercase
    {'a', "....." "....." ".@@@." "...@@" ".@@@@" ".@@.@"},
    {'b', "@...." "@...." "@@@@." "@...@" "@...@" "@@@@."},
    {'c', "....." "....." ".@@@@" "@...." "@...." ".@@@@"},
    {'d', "....@" "....@" ".@@@@" "@...@" "@...@" ".@@@@"},
    {'e', "....." "....." ".@@@." "@@@@@" "@...." ".@@@@"},
    {'f', "..@@." ".@..." "@@@@." ".@..." ".@..." ".@..."},
    {'g', "....." ".@@@@" "@...@" ".@@@@" "....@" ".@@@."},
    {'h', "@...." "@...." "@@@@." "@...@" "@...@" "@...@"},
    {'i', "..@.." "....." "..@.." "..@.." "..@.." "..@.."},
    {'j', "...@." "....." "...@." "...@." "@..@." ".@@.."},
    {'k', "@...." "@...." "@..@." "@@@.." "@..@." "@...@"},
    {'l', ".@@.." "..@.." "..@.." "..@.." "..@.." "..@@."},
    {'m', "....." "....." "@@.@@" "@.@.@" "@.@.@" "@...@"},
    {'n', "....." "....." "@@@@." "@...@" "@...@" "@...@"},
    {'o', "....." "....." ".@@@." "@...@" "@...@" ".@@@."},
    {'p', "....." "@@@@." "@...@" "@@@@." "@...." "@...."},
    {'q', "....." ".@@@@" "@...@" ".@@@@" "....@" "....@"},
    {'r', "....." "....." ".@@@@" "@...." "@...." "@...."},
    {'s', "....." ".@@@@" "@...." ".@@@." "....@" "@@@@."},
    {'t', ".@..." ".@..." "@@@@." ".@..." ".@..." "..@@."},
    {'u', "....." "....." "@...@" "@...@" "@...@" ".@@@@"},
    {'v', "....." "....." "@...@" "@...@" ".@.@." "..@.."},
    {'w', "....." "....." "@...@" "@.@.@" "@.@.@" ".@.@."},
    {'x', "....." "....." "@...@" ".@.@." ".@.@." "@...@"},
    {'y', "....." "@...@" "@...@" ".@@@@" "....@" ".@@@."},
    {'z', "....." "....." "@@@@@" "..@.." ".@..." "@@@@@"},
    // Digits
    {'0', ".@@@." "@..@@" "@.@.@" "@@..@" "@...@" ".@@@."},
    {'1', "..@.." ".@@.." "..@.." "..@.." "..@.." "@@@@@"},
    {'2', ".@@@." "@...@" "...@." "..@.." ".@..." "@@@@@"},
    {'3', ".@@@." "@...@" "..@@." "....@" "@...@" ".@@@."},
    {'4', "@...@" "@...@" "@@@@@" "....@" "....@" "....@"},
    {'5', "@@@@@" "@...." "@@@@." "....@" "....@" "@@@@."},
    {'6', ".@@@." "@...." "@@@@." "@...@" "@...@" ".@@@."},
    {'7', "@@@@@" "....@" "...@." "..@.." ".@..." ".@..."},
    {'8', ".@@@." "@...@" ".@@@." "@...@" "@...@" ".@@@."},
    {'9', ".@@@." "@...@" ".@@@@" "....@" "....@" ".@@@."},
    // Symbols
    {' ', "....." "....." "....." "....." "....." "....."},
    {',', "....." "....." "....." "....." "..@.." ".@..."},
    {'.', "....." "....." "....." "....." "....." "..@.."},
    {'!', "..@.." "..@.." "..@.." "..@.." "....." "..@.."},
    {'?', ".@@@." "@...@" "...@." "..@.." "....." "..@.."},
    {'<', "...@." "..@.." ".@..." "..@.." "...@." "....."},
    {'>', ".@..." "..@.." "...@." "..@.." ".@..." "....."},
    {'-', "....." "....." "....." "@@@@@" "....." "....."},
    {'_', "....." "....." "....." "....." "....." "@@@@@"},
    {'(', "..@.." ".@..." "@...." "@...." ".@..." "..@.."},
    {')', "..@.." "...@." "....@" "....@" "...@." "..@.."},
    {'/', "....@" "...@." "..@.." ".@..." "@...." "....."},
    {':', "....." "..@.." "....." "....." "..@.." "....."},
    {';', "....." "..@.." "....." "....." "..@.." ".@..."},
    {'\'', "..@.." "..@.." "....." "....." "....." "....."},
    {'"', ".@.@." ".@.@." "....." "....." "....." "....."},
    {'@', ".@@@." "@...@" "@.@@@" "@.@@@" "@...." ".@@@@"},
    {'*', "....." "@.@.@" ".@@@." "@@@@@" ".@@@." "@.@.@"},
    {'+', "....." "..@.." "..@.." "@@@@@" "..@.." "..@.."},
    {'=', "....." "....." "@@@@@" "....." "@@@@@" "....."},
    {'[', ".@@@." ".@..." ".@..." ".@..." ".@..." ".@@@."},
    {']', ".@@@." "...@." "...@." "...@." "...@." ".@@@."},
    {'#', ".@.@." "@@@@@" ".@.@." ".@.@." "@@@@@" ".@.@."},
};
static const int FONT_COUNT = sizeof(FONT) / sizeof(FONT[0]);

static const char* font_glyph(char ch) {
    for (int i = 0; i < FONT_COUNT; i++)
        if (FONT[i].ch == ch) return FONT[i].p;
    return nullptr;
}

static void init_centers() {
    // Bitmap text (raw pixel-art from YAML)
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

    // Font-rendered text blocks
    for (auto& tb : SC.text_blocks) {
        int cursor = 0;
        for (char ch : tb.content) {
            const char* g = font_glyph(ch);
            if (!g) { cursor += 6; continue; }
            for (int row = 0; row < 6; row++) {
                for (int col = 0; col < 5; col++) {
                    if (g[row * 5 + col] == '@') {
                        centers.push_back({
                            tb.x + (cursor + col) * tb.col_spacing,
                            tb.y - row * tb.row_spacing,
                            tb.z
                        });
                    }
                }
            }
            cursor += 6;
        }
    }
}

static inline double sdf_round_box(const Vec& p, double half_edge, double rounding) {
    double qx = std::max(std::abs(p.x) - half_edge, 0.0);
    double qy = std::max(std::abs(p.y) - half_edge, 0.0);
    double qz = std::max(std::abs(p.z) - half_edge, 0.0);
    double outer = std::sqrt(qx*qx + qy*qy + qz*qz);
    double inner = std::min(std::max({std::abs(p.x) - half_edge, std::abs(p.y) - half_edge, std::abs(p.z) - half_edge}), 0.0);
    return outer + inner - rounding;
}

static Vec sdf_normal(const Vec& p, double half_edge, double rounding) {
    const double eps = 0.001;
    double dx = sdf_round_box(p + Vec(eps,0,0), half_edge, rounding) - sdf_round_box(p - Vec(eps,0,0), half_edge, rounding);
    double dy = sdf_round_box(p + Vec(0,eps,0), half_edge, rounding) - sdf_round_box(p - Vec(0,eps,0), half_edge, rounding);
    double dz = sdf_round_box(p + Vec(0,0,eps), half_edge, rounding) - sdf_round_box(p - Vec(0,0,eps), half_edge, rounding);
    return Vec(dx, dy, dz).norm();
}

static inline double sdf_torus(const Vec& p, double major_r, double minor_r) {
    double q = std::sqrt(p.x*p.x + p.z*p.z) - major_r;
    return std::sqrt(q*q + p.y*p.y) - minor_r;
}

static Vec sdf_torus_normal(const Vec& p, double major_r, double minor_r) {
    const double eps = 0.001;
    double dx = sdf_torus(p + Vec(eps,0,0), major_r, minor_r) - sdf_torus(p - Vec(eps,0,0), major_r, minor_r);
    double dy = sdf_torus(p + Vec(0,eps,0), major_r, minor_r) - sdf_torus(p - Vec(0,eps,0), major_r, minor_r);
    double dz = sdf_torus(p + Vec(0,0,eps), major_r, minor_r) - sdf_torus(p - Vec(0,0,eps), major_r, minor_r);
    return Vec(dx, dy, dz).norm();
}

// Rounded cube sub-cubes (clustered)
struct RCubeSub {
    Vec center;
    double half_edge, rounding, bound_r, bound_r2;
    Mat3 rot, rot_inv;
};
static std::vector<RCubeSub> rcube_subs;
static Vec rcube_cluster_center;
static double rcube_cluster_bound_r2;

// Glass sub-cubes (clustered)
struct GlassSub { double x0, y0, z0, x1, y1, z1; };
static std::vector<GlassSub> glass_subs;
struct TileGlass { int start, count; };
static std::vector<TileGlass> tile_glass;

static inline uint32_t glass_rand(uint32_t& s) {
    s = s * 1664525u + 1013904223u;
    return s;
}
static inline double glass_randf(uint32_t& s) {
    return (glass_rand(s) & 0xFFFF) / 65535.0;
}

static void init_glass() {
    tile_glass.resize(SC.tile_nx * SC.tile_nz);
    bool cluster = SC.glass_cluster_min > 0 && SC.glass_cluster_max > 0;
    double max_y = 0;
    for (int ix = SC.tile_lo_x; ix <= SC.tile_hi_x; ix++) {
        for (int iz = SC.tile_lo_z; iz <= SC.tile_hi_z; iz++) {
            double th = get_tile_height(ix, iz);
            if (th < 0) {
                tile_glass[(ix - SC.tile_lo_x) * SC.tile_nz + (iz - SC.tile_lo_z)] = {0, 0};
                continue;
            }
            int start = (int)glass_subs.size();
            if (cluster) {
                uint32_t seed = (uint32_t)(ix * 374761393u + iz * 668265263u + 12345u);
                int range = SC.glass_cluster_max - SC.glass_cluster_min + 1;
                int count = SC.glass_cluster_min + (int)(glass_rand(seed) % range);
                for (int k = 0; k < count; k++) {
                    double size = SC.glass_size_min + glass_randf(seed) * (SC.glass_size_max - SC.glass_size_min);
                    double hs = size * 0.5;
                    double cx = ix + 0.5 + (glass_randf(seed) * 2.0 - 1.0) * SC.glass_spread;
                    double cz = iz + 0.5 + (glass_randf(seed) * 2.0 - 1.0) * SC.glass_spread;
                    double y1 = th + size;
                    if (y1 > max_y) max_y = y1;
                    glass_subs.push_back({cx - hs, th, cz - hs, cx + hs, y1, cz + hs});
                }
                tile_glass[(ix - SC.tile_lo_x) * SC.tile_nz + (iz - SC.tile_lo_z)] = {start, count};
            } else {
                double y1 = th + SC.glass_edge;
                if (y1 > max_y) max_y = y1;
                glass_subs.push_back({ix - 0.1, th, iz - 0.1, ix + 1.1, y1, iz + 1.1});
                tile_glass[(ix - SC.tile_lo_x) * SC.tile_nz + (iz - SC.tile_lo_z)] = {start, 1};
            }
        }
    }
    SC.glass_max_y = max_y;
}

static inline uint32_t rcube_rand(uint32_t& s) {
    s = s * 1664525u + 1013904223u;
    return s;
}
static inline double rcube_randf(uint32_t& s) {
    return (rcube_rand(s) & 0xFFFF) / 65535.0;
}

static void init_rcubes() {
    rcube_subs.clear();
    if (SC.rcube_count <= 1) {
        RCubeSub sub;
        sub.center = SC.rcube_center;
        sub.half_edge = SC.rcube_half;
        sub.rounding = SC.rcube_round;
        sub.bound_r = std::sqrt(3.0) * sub.half_edge + sub.rounding;
        sub.bound_r2 = sub.bound_r * sub.bound_r;
        sub.rot = SC.rcube_rot;
        sub.rot_inv = SC.rcube_rot_inv;
        rcube_subs.push_back(sub);
        rcube_cluster_center = SC.rcube_center;
        rcube_cluster_bound_r2 = sub.bound_r2;
        return;
    }
    uint32_t seed = (uint32_t)SC.rcube_seed;
    double max_reach = 0;
    for (int i = 0; i < SC.rcube_count; i++) {
        RCubeSub sub;
        sub.half_edge = SC.rcube_half_min + rcube_randf(seed) * (SC.rcube_half_max - SC.rcube_half_min);
        sub.rounding = sub.half_edge * SC.rcube_round_ratio;
        sub.bound_r = std::sqrt(3.0) * sub.half_edge + sub.rounding;
        sub.bound_r2 = sub.bound_r * sub.bound_r;
        // Random direction (uniform on sphere) and random distance
        double theta = rcube_randf(seed) * 2.0 * M_PI;
        double cos_phi = rcube_randf(seed) * 2.0 - 1.0;
        double sin_phi = std::sqrt(1.0 - cos_phi * cos_phi);
        double dist = rcube_randf(seed) * SC.rcube_spread;
        Vec offset(sin_phi * std::cos(theta) * dist,
                   sin_phi * std::sin(theta) * dist,
                   cos_phi * dist);
        sub.center = SC.rcube_center + offset;
        double reach = offset.len() + sub.bound_r;
        if (reach > max_reach) max_reach = reach;
        // Use original cube orientation for all sub-cubes
        sub.rot = SC.rcube_rot;
        sub.rot_inv = SC.rcube_rot_inv;
        rcube_subs.push_back(sub);
    }
    rcube_cluster_center = SC.rcube_center;
    rcube_cluster_bound_r2 = max_reach * max_reach;
}

// Thread-local glass hit info — stores actual AABB bounds
struct GlassInfo { double x0, y0, z0, x1, y1, z1; };
static thread_local GlassInfo glass_hit = {0, 0, 0, 0, 0, 0};

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

                    // Glass sub-cubes
                    int ti = (ix - SC.tile_lo_x) * SC.tile_nz + (iz - SC.tile_lo_z);
                    auto& tg = tile_glass[ti];
                    for (int gi = tg.start; gi < tg.start + tg.count; gi++) {
                        auto& gs = glass_subs[gi];
                        double gt; Vec gn;
                        if (aabb_test(ox, oy, oz, dx, dy, dz,
                                      gs.x0, gs.x1, gs.y0, gs.y1, gs.z0, gs.z1,
                                      has_dx, inv_dx, inv_dy, has_dz, inv_dz, t, gt, gn)) {
                            t = gt; n = gn; m = 4;
                            glass_hit = {gs.x0, gs.y0, gs.z0, gs.x1, gs.y1, gs.z1};
                        }
                    }
                }
            }
        }
    }

    // Text spheres + glass spheres + reflective inner spheres + emissive tori
    {
        double gsph_r = SC.sph_r * SC.gsph_radius_mult;
        double gsph_r2 = gsph_r * gsph_r;
        double rsph_r = SC.sph_r * SC.rsph_radius_mult;
        double rsph_r2 = rsph_r * rsph_r;
        double torus_major = SC.sph_r * SC.torus_major_mult;
        double torus_bound_r = torus_major + SC.torus_minor;
        double torus_bound_r2 = torus_bound_r * torus_bound_r;
        for (int ci = 0; ci < (int)centers.size(); ci++) {
            const auto& c = centers[ci];
            double px = ox - c.x, py = oy - c.y, pz = oz - c.z;
            double b = px*dx + py*dy + pz*dz;
            double d2 = px*px + py*py + pz*pz;

            // Text sphere (m=3)
            double disc = b*b - (d2 - SC.sph_r2);
            if (disc >= 0) {
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

            // Glass sphere (m=7)
            double gdisc = b*b - (d2 - gsph_r2);
            if (gdisc >= 0) {
                double gsq = std::sqrt(gdisc);
                double gtt = -b - gsq;
                if (gtt < 0.01) gtt = -b + gsq;
                if (gtt > 0.01 && gtt < t) {
                    t = gtt;
                    double inv_r = 1.0 / gsph_r;
                    n = Vec((ox + dx*gtt - c.x) * inv_r, (oy + dy*gtt - c.y) * inv_r, (oz + dz*gtt - c.z) * inv_r);
                    m = 7;
                    hit_sphere_idx = ci;
                }
            }

            // Reflective inner sphere (m=9)
            double rdisc = b*b - (d2 - rsph_r2);
            if (rdisc >= 0) {
                double rsq = std::sqrt(rdisc);
                double rtt = -b - rsq;
                if (rtt < 0.01) rtt = -b + rsq;
                if (rtt > 0.01 && rtt < t) {
                    t = rtt;
                    double inv_r = 1.0 / rsph_r;
                    n = Vec((ox + dx*rtt - c.x) * inv_r, (oy + dy*rtt - c.y) * inv_r, (oz + dz*rtt - c.z) * inv_r);
                    m = 9;
                    hit_sphere_idx = ci;
                }
            }

            // Emissive torus (m=8) — SDF ray march
            double tdisc = b*b - (d2 - torus_bound_r2);
            if (tdisc >= 0) {
                double tsq = std::sqrt(tdisc);
                double t_near = -b - tsq;
                double t_far = -b + tsq;
                if (t_far > 0.01 && t_near < t) {
                    double tt = std::max(t_near, 0.01);
                    for (int step = 0; step < 48; step++) {
                        Vec lp(ox + dx*tt - c.x, oy + dy*tt - c.y, oz + dz*tt - c.z);
                        double dist = sdf_torus(lp, torus_major, SC.torus_minor);
                        if (dist < 0.001) {
                            if (tt < t) {
                                t = tt;
                                n = sdf_torus_normal(lp, torus_major, SC.torus_minor);
                                m = 8;
                                hit_sphere_idx = ci;
                            }
                            break;
                        }
                        tt += dist;
                        if (tt > std::min(t_far, t)) break;
                    }
                }
            }
        }
    }

    // Rounded reflective cubes (SDF ray march, clustered)
    {
        double cpx = ox - rcube_cluster_center.x, cpy = oy - rcube_cluster_center.y, cpz = oz - rcube_cluster_center.z;
        double cb = cpx*dx + cpy*dy + cpz*dz;
        double cdisc = cb*cb - (cpx*cpx + cpy*cpy + cpz*cpz - rcube_cluster_bound_r2);
        if (cdisc >= 0) {
            double csq = std::sqrt(cdisc);
            double ct_far = -cb + csq;
            if (ct_far > 0.01) {
                for (const auto& sub : rcube_subs) {
                    double spx = ox - sub.center.x, spy = oy - sub.center.y, spz = oz - sub.center.z;
                    double sb = spx*dx + spy*dy + spz*dz;
                    double sdisc = sb*sb - (spx*spx + spy*spy + spz*spz - sub.bound_r2);
                    if (sdisc < 0) continue;
                    double ssq = std::sqrt(sdisc);
                    double t_near = -sb - ssq;
                    double t_far = -sb + ssq;
                    if (t_far <= 0.01 || t_near >= t) continue;
                    Vec local_o = sub.rot.mul(Vec(ox, oy, oz) - sub.center);
                    Vec local_d = sub.rot.mul(Vec(dx, dy, dz));
                    double tt = std::max(t_near, 0.01);
                    for (int step = 0; step < 64; step++) {
                        Vec lp = local_o + local_d * tt;
                        double dist = sdf_round_box(lp, sub.half_edge, sub.rounding);
                        if (dist < 0.001) {
                            if (tt < t) {
                                t = tt;
                                n = sub.rot_inv.mul(sdf_normal(lp, sub.half_edge, sub.rounding));
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

    // Emissive torus
    if (m == 8) return torus_colors[hit_sphere_idx] * SC.torus_emissive_intensity;

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

    if (m == 7 && depth < SC.gsph_max_depth) {
        // Glass sphere — refraction + reflection
        double cos_i = -(d.dot(n));
        Vec nn = n;
        if (cos_i < 0) { cos_i = -cos_i; nn = n * -1; }
        double fresnel = SC.gsph_r0 + (1 - SC.gsph_r0) * std::pow(1 - cos_i, 5);
        Vec r = d - nn * (2 * d.dot(nn));

        Vec transmitted(0, 0, 0);
        double eta = 1.0 / SC.gsph_ior;
        double k = 1.0 - eta * eta * (1.0 - cos_i * cos_i);
        if (k > 0) {
            Vec refr = d * eta + nn * (eta * cos_i - std::sqrt(k));
            // Find exit point: ray from inside sphere hits back surface
            const auto& c = centers[hit_sphere_idx];
            double gsph_r = SC.sph_r * SC.gsph_radius_mult;
            Vec entry = hit + refr * 0.02;
            double epx = entry.x - c.x, epy = entry.y - c.y, epz = entry.z - c.z;
            double eb = epx*refr.x + epy*refr.y + epz*refr.z;
            double edisc = eb*eb - (epx*epx + epy*epy + epz*epz - gsph_r*gsph_r);
            if (edisc > 0) {
                double exit_t = -eb + std::sqrt(edisc);
                if (exit_t > 0) {
                    Vec exit_pt = entry + refr * exit_t;
                    Vec en = Vec(exit_pt.x - c.x, exit_pt.y - c.y, exit_pt.z - c.z) * (1.0 / gsph_r);
                    double cos_o = -(refr.dot(en));
                    if (cos_o < 0) { cos_o = -cos_o; en = en * -1; }
                    double eta2 = SC.gsph_ior;
                    double k2 = 1.0 - eta2 * eta2 * (1.0 - cos_o * cos_o);
                    if (k2 > 0) {
                        Vec exit_refr = refr * eta2 + en * (eta2 * cos_o - std::sqrt(k2));
                        transmitted = trace(exit_pt + exit_refr * 0.02, exit_refr, depth + 1);
                    } else {
                        Vec refl_int = refr - en * (2 * refr.dot(en));
                        transmitted = trace(exit_pt + refl_int * 0.02, refl_int, depth + 1);
                    }
                    double pl = exit_t * SC.gsph_tint_scale;
                    transmitted = Vec(transmitted.x * (1 - pl*SC.gsph_tint_rgb.x),
                                      transmitted.y * (1 - pl*SC.gsph_tint_rgb.y),
                                      transmitted.z * (1 - pl*SC.gsph_tint_rgb.z));
                }
            }
        } else {
            fresnel = 1.0;
        }
        Vec reflected = fresnel > 0.02 ? trace(hit + nn * 0.02, r, depth + 1) : Vec();
        Vec base = reflected * fresnel + transmitted * (1.0 - fresnel);
        double sp1 = std::pow(std::max(0.0, light1.dot(r)), SC.gsph_spec_power) * s1;
        double sp2 = std::pow(std::max(0.0, light2.dot(r)), SC.gsph_spec_power) * s2;
        return base + L1 * (sp1 * SC.gsph_spec_intensity) + L2 * (sp2 * SC.gsph_spec_intensity);
    }
    if (m == 7) {
        return depth < SC.gsph_passthrough_depth ? trace(hit + d * 0.02, d, depth + 1) : Vec();
    }
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
            double gx0 = gi.x0, gx1 = gi.x1;
            double gy0 = gi.y0, gy1 = gi.y1;
            double gz0 = gi.z0, gz1 = gi.z1;
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
    if (m == 5 || m == 9) {
        // Chrome reflective surface (rounded cubes + inner spheres)
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
    if (version != "nanore@2") {
        fprintf(stderr, "Error: unsupported scene version '%s' (expected nanore@2)\n", version.c_str());
        return 1;
    }
    SC.load_from_yaml(yaml);
    printf("Loaded scene: %s\n", scene_file);

    init_tiles();
    init_glass();
    init_rcubes();
    init_centers();

    // Generate random emissive colors for tori
    torus_colors.clear();
    {
        uint32_t tseed = 777u;
        for (size_t i = 0; i < centers.size(); i++) {
            torus_colors.push_back(Vec(
                rcube_randf(tseed) * 2.0 + 0.5,
                rcube_randf(tseed) * 2.0 + 0.5,
                rcube_randf(tseed) * 2.0 + 0.5));
        }
    }

    int W = SC.width, H = SC.height;
    printf("nanore ray tracer\n");
    printf("  Resolution: %dx%d, %d samples/pixel\n", W, H, SC.samples);
    printf("  %zu spheres (r=%.1f), %d displaced cubes + glass\n", centers.size(), SC.sph_r, SC.tile_nx * SC.tile_nz);
    if (SC.rcube_count > 1)
        printf("  %d rounded reflective sub-cubes (edge %.1f-%.1f, spread=%.1f)\n",
               SC.rcube_count, SC.rcube_half_min * 2, SC.rcube_half_max * 2, SC.rcube_spread);
    else
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

    // Also save as renders/latest.bmp
    FILE* fl = fopen("renders/latest.bmp", "wb");
    if (fl) {
        fwrite(bmp_header, 1, 54, fl);
        for (int y = H - 1; y >= 0; y--) {
            for (int x = 0; x < W; x++) {
                int idx = (y * W + x) * 3;
                uint8_t bgr[3] = {pixels[idx+2], pixels[idx+1], pixels[idx]};
                fwrite(bgr, 1, 3, fl);
            }
            if (pad) fwrite(padding, 1, pad, fl);
        }
        fclose(fl);
        printf("  Saved to renders/latest.bmp\n");
    }

    delete[] hdr;
    delete[] pixels;
    return 0;
}
