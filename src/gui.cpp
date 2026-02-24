#include "gui.h"
#include <SDL2/SDL.h>
#include <GL/glew.h>
#include "imgui.h"
#include "imgui_internal.h"
#include "imgui_impl_sdl2.h"
#include "imgui_impl_opengl3.h"
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <sys/stat.h>

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

struct GuiState {
    SDL_Window* window = nullptr;
    SDL_GLContext gl_ctx = nullptr;

    // Scene file
    std::string yaml_path;
    std::string yaml_content;
    time_t yaml_mtime = 0;
    bool file_changed = false;
    bool path_dirty = false;  // set when yaml_path changes, to sync path_buf in panel

    // Editor buffer (64 KB max)
    static constexpr int BUF_SZ = 65536;
    char yaml_buf[BUF_SZ] = {};

    // Scene / render state
    bool scene_loaded = false;  // true once init_scene_from_yaml_str succeeds (wireframes available)
    bool rendering = false;
    bool render_valid = false;
    int render_w = 0, render_h = 0;

    // GL textures
    GLuint tex_render = 0;
    GLuint tex_depth = 0;

    // Render/depth viewport pan/zoom (shared between color and depth)
    float rv_zoom = 1.0f, rv_pan_x = 0.0f, rv_pan_y = 0.0f;

    // Wireframe view pan/zoom
    float top_zoom = 8.0f, top_pan_x = 0.0f, top_pan_y = 0.0f;
    float front_zoom = 8.0f, front_pan_x = 0.0f, front_pan_y = 0.0f;

    // File watcher timing
    Uint64 last_file_check = 0;

    // First frame docking layout
    bool first_frame = true;
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static time_t file_mtime(const std::string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return 0;
    return st.st_mtime;
}

static std::string read_file(const std::string& path) {
    std::ifstream f(path);
    if (!f) return "";
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

static GLuint make_texture(int w, int h) {
    GLuint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
    return tex;
}

static void upload_rgb(GLuint tex, int w, int h, const unsigned char* px) {
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, px);
}

static void upload_depth(GLuint tex, int w, int h, const float* zbuf) {
    float zmin = zbuf[0], zmax = zbuf[0];
    for (int i = 1; i < w * h; i++) {
        if (zbuf[i] < zmin) zmin = zbuf[i];
        if (zbuf[i] > zmax) zmax = zbuf[i];
    }
    float zrange = (zmax - zmin > 1e-6f) ? (zmax - zmin) : 1.0f;
    std::vector<unsigned char> gray(w * h * 3);
    for (int i = 0; i < w * h; i++) {
        unsigned char v = (unsigned char)(std::clamp(1.0f - (zbuf[i] - zmin) / zrange, 0.0f, 1.0f) * 255.0f);
        gray[i * 3] = gray[i * 3 + 1] = gray[i * 3 + 2] = v;
    }
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, gray.data());
}

// ---------------------------------------------------------------------------
// File watching
// ---------------------------------------------------------------------------

static void check_file(GuiState& s) {
    if (s.yaml_path.empty()) return;
    Uint64 now = SDL_GetTicks64();
    if (now - s.last_file_check < 500) return;
    s.last_file_check = now;
    time_t mt = file_mtime(s.yaml_path);
    if (mt != s.yaml_mtime && mt != 0) {
        s.yaml_mtime = mt;
        s.yaml_content = read_file(s.yaml_path);
        std::strncpy(s.yaml_buf, s.yaml_content.c_str(), GuiState::BUF_SZ - 1);
        s.yaml_buf[GuiState::BUF_SZ - 1] = '\0';
        s.file_changed = true;
    }
}

// ---------------------------------------------------------------------------
// Render trigger
// ---------------------------------------------------------------------------

static void init_scene_for_wireframe(GuiState& s) {
    if (s.yaml_buf[0] == '\0') return;
    if (init_scene_from_yaml_str(s.yaml_buf) == 0)
        s.scene_loaded = true;
}

static void trigger_render(GuiState& s) {
    if (s.yaml_buf[0] == '\0') return;
    if (s.rendering) {
        cancel_render();
        while (!is_render_done()) SDL_Delay(5);
    }
    // init_scene runs again inside render thread, but we already called it
    // for wireframe. The render thread will re-init (safe, same data).
    start_render(s.yaml_buf);
    s.rendering = true;
    s.file_changed = false;
}

static void save_file(GuiState& s) {
    if (s.yaml_path.empty()) return;
    std::ofstream f(s.yaml_path);
    if (f) {
        f << s.yaml_buf;
        f.close();
        s.yaml_mtime = file_mtime(s.yaml_path);
    }
}

static void load_file(GuiState& s, const char* path) {
    s.yaml_path = path;
    s.yaml_content = read_file(s.yaml_path);
    s.yaml_mtime = file_mtime(s.yaml_path);
    std::strncpy(s.yaml_buf, s.yaml_content.c_str(), GuiState::BUF_SZ - 1);
    s.yaml_buf[GuiState::BUF_SZ - 1] = '\0';
    s.path_dirty = true;
    init_scene_for_wireframe(s);
    trigger_render(s);
}

// ---------------------------------------------------------------------------
// Wireframe drawing
// ---------------------------------------------------------------------------

static ImVec2 proj_top(double x, double z, float zoom, float px, float py,
                       ImVec2 org, ImVec2 sz) {
    return ImVec2(org.x + sz.x * 0.5f + (float)x * zoom + px,
                  org.y + sz.y * 0.5f + (float)z * zoom + py);
}

static ImVec2 proj_front(double x, double y, float zoom, float px, float py,
                         ImVec2 org, ImVec2 sz) {
    return ImVec2(org.x + sz.x * 0.5f + (float)x * zoom + px,
                  org.y + sz.y * 0.5f - (float)y * zoom + py);
}

// Draw sphere as wireframe circle
static void draw_circle(ImDrawList* dl, ImVec2 c, float r, ImU32 col, int segs = 16) {
    if (r < 0.5f) r = 0.5f;
    dl->AddCircle(c, r, col, segs);
}

// Draw axis-aligned box as wireframe rect
static void draw_rect(ImDrawList* dl, ImVec2 a, ImVec2 b, ImU32 col) {
    dl->AddRect(a, b, col);
}

typedef ImVec2 (*ProjFn)(double, double, float, float, float, ImVec2, ImVec2);

static void draw_wireframe(ImDrawList* dl, ImVec2 org, ImVec2 sz,
                            float zoom, float pan_x, float pan_y,
                            ProjFn proj, bool is_top) {
    ImU32 col_sph   = IM_COL32(0, 200, 255, 180);
    ImU32 col_cube  = IM_COL32(255, 200, 0, 180);
    ImU32 col_glass = IM_COL32(100, 200, 255, 80);
    ImU32 col_disc  = IM_COL32(255, 50, 50, 200);
    ImU32 col_floor = IM_COL32(80, 80, 80, 100);
    ImU32 col_light = IM_COL32(255, 255, 100, 255);
    ImU32 col_cam   = IM_COL32(0, 255, 0, 255);
    ImU32 col_ss    = IM_COL32(0, 180, 220, 150);

    dl->PushClipRect(org, ImVec2(org.x + sz.x, org.y + sz.y), true);

    // Grid lines (subtle)
    {
        ImU32 col_grid = IM_COL32(40, 40, 50, 60);
        float step = 5.0f * zoom;
        if (step > 4.0f) {
            float cx = org.x + sz.x * 0.5f + pan_x;
            float cy = org.y + sz.y * 0.5f + pan_y;
            float left = org.x, right = org.x + sz.x;
            float top = org.y, bottom = org.y + sz.y;
            for (float gx = cx; gx <= right; gx += step)
                dl->AddLine(ImVec2(gx, top), ImVec2(gx, bottom), col_grid);
            for (float gx = cx - step; gx >= left; gx -= step)
                dl->AddLine(ImVec2(gx, top), ImVec2(gx, bottom), col_grid);
            for (float gy = cy; gy <= bottom; gy += step)
                dl->AddLine(ImVec2(left, gy), ImVec2(right, gy), col_grid);
            for (float gy = cy - step; gy >= top; gy -= step)
                dl->AddLine(ImVec2(left, gy), ImVec2(right, gy), col_grid);
        }
        // Origin crosshair
        ImU32 col_axis = IM_COL32(60, 60, 70, 100);
        float cx = org.x + sz.x * 0.5f + pan_x;
        float cy = org.y + sz.y * 0.5f + pan_y;
        dl->AddLine(ImVec2(cx, org.y), ImVec2(cx, org.y + sz.y), col_axis);
        dl->AddLine(ImVec2(org.x, cy), ImVec2(org.x + sz.x, cy), col_axis);
    }

    // Floor bounds
    {
        int lox, hix, loz, hiz;
        get_floor_bounds(&lox, &hix, &loz, &hiz);
        double a1 = lox, a2 = hix + 1.0;
        double b1 = is_top ? (double)loz : 0.0;
        double b2 = is_top ? (double)(hiz + 1.0) : 0.0;
        if (is_top) {
            draw_rect(dl, proj(a1, b1, zoom, pan_x, pan_y, org, sz),
                          proj(a2, b2, zoom, pan_x, pan_y, org, sz), col_floor);
        } else {
            // Front: floor is a horizontal line at y≈0
            dl->AddLine(proj(a1, 0, zoom, pan_x, pan_y, org, sz),
                        proj(a2, 0, zoom, pan_x, pan_y, org, sz), col_floor, 2.0f);
        }
    }

    // Text spheres
    {
        double sph_r = get_sphere_radius();
        int n = get_sphere_count();
        for (int i = 0; i < n; i++) {
            double x, y, z;
            get_sphere_center(i, &x, &y, &z);
            ImVec2 c = is_top ? proj(x, z, zoom, pan_x, pan_y, org, sz)
                              : proj(x, y, zoom, pan_x, pan_y, org, sz);
            draw_circle(dl, c, (float)(sph_r * zoom), col_sph, 8);
        }
    }

    // Rounded cubes
    {
        int n = get_rcube_count();
        for (int i = 0; i < n; i++) {
            double cx, cy, cz, half, round_;
            get_rcube_info(i, &cx, &cy, &cz, &half, &round_);
            double ext = half + round_;
            if (is_top) {
                draw_rect(dl, proj(cx - ext, cz - ext, zoom, pan_x, pan_y, org, sz),
                              proj(cx + ext, cz + ext, zoom, pan_x, pan_y, org, sz), col_cube);
            } else {
                draw_rect(dl, proj(cx - ext, cy + ext, zoom, pan_x, pan_y, org, sz),
                              proj(cx + ext, cy - ext, zoom, pan_x, pan_y, org, sz), col_cube);
            }
        }
    }

    // Glass sub-cubes
    {
        int n = get_glass_sub_count();
        for (int i = 0; i < n; i++) {
            double x0, y0, z0, x1, y1, z1;
            get_glass_sub(i, &x0, &y0, &z0, &x1, &y1, &z1);
            if (is_top)
                draw_rect(dl, proj(x0, z0, zoom, pan_x, pan_y, org, sz),
                              proj(x1, z1, zoom, pan_x, pan_y, org, sz), col_glass);
            else
                draw_rect(dl, proj(x0, y1, zoom, pan_x, pan_y, org, sz),
                              proj(x1, y0, zoom, pan_x, pan_y, org, sz), col_glass);
        }
    }

    // Disc
    {
        double cx, cy, cz, r;
        get_disc_info(&cx, &cy, &cz, &r);
        if (is_top) {
            draw_circle(dl, proj(cx, cz, zoom, pan_x, pan_y, org, sz),
                        (float)(r * zoom), col_disc, 32);
        } else {
            // Front: disc appears as a horizontal line at its y position
            dl->AddLine(proj(cx - r, cy, zoom, pan_x, pan_y, org, sz),
                        proj(cx + r, cy, zoom, pan_x, pan_y, org, sz), col_disc, 2.0f);
        }
    }

    // Lights
    {
        int n = get_light_count();
        for (int i = 0; i < n; i++) {
            double lx, ly, lz;
            get_light_info(i, &lx, &ly, &lz);
            ImVec2 lp = is_top ? proj(lx, lz, zoom, pan_x, pan_y, org, sz)
                               : proj(lx, ly, zoom, pan_x, pan_y, org, sz);
            dl->AddCircleFilled(lp, 4.0f, col_light);
            // Cross marker
            dl->AddLine(ImVec2(lp.x - 6, lp.y), ImVec2(lp.x + 6, lp.y), col_light, 1.5f);
            dl->AddLine(ImVec2(lp.x, lp.y - 6), ImVec2(lp.x, lp.y + 6), col_light, 1.5f);
        }
    }

    // Camera
    {
        double cpx, cpy, cpz, ctx_, cty, ctz;
        get_camera_info(&cpx, &cpy, &cpz, &ctx_, &cty, &ctz);
        ImVec2 cam = is_top ? proj(cpx, cpz, zoom, pan_x, pan_y, org, sz)
                            : proj(cpx, cpy, zoom, pan_x, pan_y, org, sz);
        ImVec2 tgt = is_top ? proj(ctx_, ctz, zoom, pan_x, pan_y, org, sz)
                            : proj(ctx_, cty, zoom, pan_x, pan_y, org, sz);
        dl->AddCircleFilled(cam, 5.0f, col_cam);
        dl->AddLine(cam, tgt, col_cam, 2.0f);
        // FOV triangle hint
        double dx = (is_top ? ctx_ - cpx : ctx_ - cpx);
        double dy = (is_top ? ctz - cpz : cty - cpy);
        double len = std::sqrt(dx * dx + dy * dy);
        if (len > 0.01) {
            double nx = -dy / len, ny = dx / len;
            double fov_spread = len * 0.3;
            ImVec2 fl = is_top ? proj(ctx_ + nx * fov_spread, ctz + ny * fov_spread, zoom, pan_x, pan_y, org, sz)
                               : proj(ctx_ + nx * fov_spread, cty + ny * fov_spread, zoom, pan_x, pan_y, org, sz);
            ImVec2 fr = is_top ? proj(ctx_ - nx * fov_spread, ctz - ny * fov_spread, zoom, pan_x, pan_y, org, sz)
                               : proj(ctx_ - nx * fov_spread, cty - ny * fov_spread, zoom, pan_x, pan_y, org, sz);
            dl->AddLine(cam, fl, IM_COL32(0, 255, 0, 80));
            dl->AddLine(cam, fr, IM_COL32(0, 255, 0, 80));
        }
    }

    // Standalone spheres
    {
        int n = get_standalone_sphere_count();
        for (int i = 0; i < n; i++) {
            double cx, cy, cz, gr, cr;
            get_standalone_sphere(i, &cx, &cy, &cz, &gr, &cr);
            ImVec2 c = is_top ? proj(cx, cz, zoom, pan_x, pan_y, org, sz)
                              : proj(cx, cy, zoom, pan_x, pan_y, org, sz);
            draw_circle(dl, c, (float)(gr * zoom), col_ss, 16);
            draw_circle(dl, c, (float)(cr * zoom), col_sph, 12);
        }
    }

    dl->PopClipRect();
}

// ---------------------------------------------------------------------------
// Native file dialog (Linux: zenity, kdialog fallback)
// ---------------------------------------------------------------------------

static std::string open_file_dialog() {
    // Try zenity first, then kdialog
    const char* cmds[] = {
        "zenity --file-selection --title='Open Scene YAML' --file-filter='YAML files|*.yaml *.yml' 2>/dev/null",
        "kdialog --getopenfilename . 'YAML files (*.yaml *.yml)' 2>/dev/null",
    };
    for (auto cmd : cmds) {
        FILE* fp = popen(cmd, "r");
        if (!fp) continue;
        char buf[1024] = {};
        if (fgets(buf, sizeof(buf), fp)) {
            int status = pclose(fp);
            if (status == 0) {
                // Strip trailing newline
                size_t len = strlen(buf);
                while (len > 0 && (buf[len-1] == '\n' || buf[len-1] == '\r')) buf[--len] = '\0';
                if (len > 0) return std::string(buf);
            }
        } else {
            pclose(fp);
        }
    }
    return "";
}

// ---------------------------------------------------------------------------
// Panel: YAML editor
// ---------------------------------------------------------------------------

static void panel_yaml(GuiState& s) {
    ImGui::Begin("Scene YAML");

    // File path + Open/Load buttons
    static char path_buf[1024] = {};
    // Sync path_buf when file path changes (load, drag-drop, file watcher reload)
    if (s.path_dirty) {
        std::strncpy(path_buf, s.yaml_path.c_str(), sizeof(path_buf) - 1);
        path_buf[sizeof(path_buf) - 1] = '\0';
        s.path_dirty = false;
    }

    if (ImGui::Button("Open...")) {
        std::string chosen = open_file_dialog();
        if (!chosen.empty()) {
            std::strncpy(path_buf, chosen.c_str(), sizeof(path_buf) - 1);
            path_buf[sizeof(path_buf) - 1] = '\0';
            load_file(s, path_buf);
        }
    }
    ImGui::SameLine();
    ImGui::SetNextItemWidth(-130);
    if (ImGui::InputText("##path", path_buf, sizeof(path_buf), ImGuiInputTextFlags_EnterReturnsTrue)) {
        load_file(s, path_buf);
    }
    ImGui::SameLine();
    if (ImGui::Button("Load", ImVec2(55, 0))) {
        load_file(s, path_buf);
    }
    ImGui::SameLine();
    if (ImGui::Button("Save", ImVec2(-1, 0))) {
        save_file(s);
    }

    // Render / Cancel + progress bar (single row)
    bool active = s.rendering && !is_render_done();
    if (active) {
        if (ImGui::Button("Cancel", ImVec2(80, 0)))
            cancel_render();
    } else {
        if (ImGui::Button("Render", ImVec2(80, 0)))
            trigger_render(s);
    }
    ImGui::SameLine();
    float pct = active ? get_progress() / 100.0f : (s.render_valid ? 1.0f : 0.0f);
    char overlay[32];
    if (active) snprintf(overlay, sizeof(overlay), "%d%%", (int)(pct * 100));
    else if (s.render_valid) snprintf(overlay, sizeof(overlay), "%dx%d", s.render_w, s.render_h);
    else snprintf(overlay, sizeof(overlay), "--");
    ImGui::ProgressBar(pct, ImVec2(-1, 0), overlay);

    ImGui::Separator();

    // YAML text editor
    ImVec2 avail = ImGui::GetContentRegionAvail();
    ImGui::InputTextMultiline("##yaml", s.yaml_buf, GuiState::BUF_SZ,
                              avail, ImGuiInputTextFlags_AllowTabInput);

    ImGui::End();
}

// ---------------------------------------------------------------------------
// Panel: image viewport (render or depth)
// ---------------------------------------------------------------------------

static void panel_image(const char* name, GLuint tex, int w, int h,
                        float& zoom, float& pan_x, float& pan_y) {
    ImGui::Begin(name);
    if (tex && w > 0 && h > 0) {
        ImVec2 avail = ImGui::GetContentRegionAvail();
        ImVec2 cursor = ImGui::GetCursorScreenPos();
        float aspect = (float)w / (float)h;

        // Base size: fit image in available space
        float base_w = avail.x, base_h = avail.x / aspect;
        if (base_h > avail.y) { base_h = avail.y; base_w = avail.y * aspect; }

        // Apply zoom
        float dw = base_w * zoom, dh = base_h * zoom;

        // Compute top-left with pan (centered when zoom=1, pan=0)
        float cx = cursor.x + avail.x * 0.5f;
        float cy = cursor.y + avail.y * 0.5f;
        float x0 = cx - dw * 0.5f + pan_x;
        float y0 = cy - dh * 0.5f + pan_y;

        // UV clipping to visible region
        float clip_x0 = cursor.x, clip_y0 = cursor.y;
        float clip_x1 = cursor.x + avail.x, clip_y1 = cursor.y + avail.y;

        // UV coordinates for the visible portion
        float uv_x0 = (clip_x0 - x0) / dw, uv_y0 = (clip_y0 - y0) / dh;
        float uv_x1 = (clip_x1 - x0) / dw, uv_y1 = (clip_y1 - y0) / dh;
        uv_x0 = std::max(0.0f, uv_x0); uv_y0 = std::max(0.0f, uv_y0);
        uv_x1 = std::min(1.0f, uv_x1); uv_y1 = std::min(1.0f, uv_y1);

        // Screen coordinates for the visible portion
        float sx0 = std::max(clip_x0, x0), sy0 = std::max(clip_y0, y0);
        float sx1 = std::min(clip_x1, x0 + dw), sy1 = std::min(clip_y1, y0 + dh);

        if (sx1 > sx0 && sy1 > sy0) {
            ImDrawList* dl = ImGui::GetWindowDrawList();
            dl->AddImage((ImTextureID)(intptr_t)tex,
                         ImVec2(sx0, sy0), ImVec2(sx1, sy1),
                         ImVec2(uv_x0, uv_y0), ImVec2(uv_x1, uv_y1));
        }

        // Interaction area
        ImGui::InvisibleButton("##img_interact", avail);
        ImGuiIO& io = ImGui::GetIO();
        if (ImGui::IsItemHovered()) {
            if (io.MouseWheel != 0.0f) {
                float old_zoom = zoom;
                zoom *= 1.0f + io.MouseWheel * 0.15f;
                zoom = std::clamp(zoom, 0.1f, 50.0f);
                // Zoom towards mouse pointer
                float mx = io.MousePos.x - cx - pan_x;
                float my = io.MousePos.y - cy - pan_y;
                float scale = zoom / old_zoom;
                pan_x -= mx * (scale - 1.0f);
                pan_y -= my * (scale - 1.0f);
            }
            if (ImGui::IsMouseDragging(ImGuiMouseButton_Left) ||
                ImGui::IsMouseDragging(ImGuiMouseButton_Middle)) {
                pan_x += io.MouseDelta.x;
                pan_y += io.MouseDelta.y;
            }
            if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                zoom = 1.0f; pan_x = 0; pan_y = 0;
            }
        }
    }
    ImGui::End();
}

// ---------------------------------------------------------------------------
// Panel: wireframe viewport
// ---------------------------------------------------------------------------

static void panel_wireframe(const char* name, float& zoom, float& pan_x, float& pan_y,
                             bool render_valid, bool is_top) {
    ImGui::Begin(name);
    ImVec2 pos = ImGui::GetCursorScreenPos();
    ImVec2 avail = ImGui::GetContentRegionAvail();
    ImDrawList* dl = ImGui::GetWindowDrawList();

    // Background
    dl->AddRectFilled(pos, ImVec2(pos.x + avail.x, pos.y + avail.y), IM_COL32(15, 15, 25, 255));

    if (render_valid) {
        ProjFn fn = is_top ? proj_top : proj_front;
        draw_wireframe(dl, pos, avail, zoom, pan_x, pan_y, fn, is_top);
    }

    // Legend
    {
        float lx = pos.x + 8, ly = pos.y + avail.y - 16;
        ImGui::PushStyleVar(ImGuiStyleVar_Alpha, 0.5f);
        dl->AddText(ImVec2(lx, ly), IM_COL32(180, 180, 180, 128),
                    is_top ? "Top (X-Z)" : "Front (X-Y)");
        ImGui::PopStyleVar();
    }

    // Interaction
    ImGui::InvisibleButton(is_top ? "##top_interact" : "##front_interact", avail);
    ImGuiIO& io = ImGui::GetIO();
    if (ImGui::IsItemHovered()) {
        if (io.MouseWheel != 0.0f) {
            float old_zoom = zoom;
            zoom *= 1.0f + io.MouseWheel * 0.15f;
            if (zoom < 0.5f) zoom = 0.5f;
            if (zoom > 200.0f) zoom = 200.0f;
            // Zoom towards mouse pointer
            float mx = io.MousePos.x - (pos.x + avail.x * 0.5f) - pan_x;
            float my = io.MousePos.y - (pos.y + avail.y * 0.5f) - pan_y;
            float scale = zoom / old_zoom;
            pan_x -= mx * (scale - 1.0f);
            pan_y -= my * (scale - 1.0f);
        }
        if (ImGui::IsMouseDragging(ImGuiMouseButton_Left) ||
            ImGui::IsMouseDragging(ImGuiMouseButton_Middle)) {
            pan_x += io.MouseDelta.x;
            pan_y += io.MouseDelta.y;
        }
        if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
            zoom = 8.0f; pan_x = 0; pan_y = 0;
        }
    }

    ImGui::End();
}

// ---------------------------------------------------------------------------
// Default docking layout (first frame only)
// ---------------------------------------------------------------------------

static void setup_default_layout(ImGuiID dockspace_id) {
    ImGui::DockBuilderRemoveNode(dockspace_id);
    ImGui::DockBuilderAddNode(dockspace_id, ImGuiDockNodeFlags_DockSpace);
    ImGui::DockBuilderSetNodeSize(dockspace_id, ImGui::GetMainViewport()->Size);

    // Split: left (yaml) | right (viewports)
    ImGuiID left, right;
    ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Left, 0.30f, &left, &right);

    // Split right into top-row and bottom-row
    ImGuiID right_top, right_bot;
    ImGui::DockBuilderSplitNode(right, ImGuiDir_Up, 0.5f, &right_top, &right_bot);

    // Split top-row: render | depth
    ImGuiID rt_left, rt_right;
    ImGui::DockBuilderSplitNode(right_top, ImGuiDir_Left, 0.5f, &rt_left, &rt_right);

    // Split bottom-row: top-view | front-view
    ImGuiID rb_left, rb_right;
    ImGui::DockBuilderSplitNode(right_bot, ImGuiDir_Left, 0.5f, &rb_left, &rb_right);

    ImGui::DockBuilderDockWindow("Scene YAML", left);
    ImGui::DockBuilderDockWindow("Render", rt_left);
    ImGui::DockBuilderDockWindow("Depth", rt_right);
    ImGui::DockBuilderDockWindow("Top View", rb_left);
    ImGui::DockBuilderDockWindow("Front View", rb_right);

    ImGui::DockBuilderFinish(dockspace_id);
}

// ---------------------------------------------------------------------------
// Main GUI loop
// ---------------------------------------------------------------------------

int gui_main(int argc, char* argv[]) {
    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) != 0) {
        fprintf(stderr, "SDL_Init error: %s\n", SDL_GetError());
        return 1;
    }

    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);

    SDL_Window* window = SDL_CreateWindow("nanore",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        1600, 1000,
        SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
    if (!window) {
        fprintf(stderr, "SDL_CreateWindow error: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_GLContext gl_ctx = SDL_GL_CreateContext(window);
    SDL_GL_MakeCurrent(window, gl_ctx);
    SDL_GL_SetSwapInterval(1);

    glewExperimental = GL_TRUE;
    GLenum glewErr = glewInit();
    if (glewErr != GLEW_OK) {
        fprintf(stderr, "glewInit error: %s\n", glewGetErrorString(glewErr));
        SDL_GL_DeleteContext(gl_ctx);
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    // Enable drag-and-drop
    SDL_EventState(SDL_DROPFILE, SDL_ENABLE);

    // ImGui setup
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
    io.IniFilename = nullptr; // Don't save layout to file

    // Dark theme with neon accent (matching web UI)
    ImGui::StyleColorsDark();
    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowRounding = 2.0f;
    style.FrameRounding = 2.0f;
    style.Colors[ImGuiCol_WindowBg] = ImVec4(0.05f, 0.05f, 0.10f, 1.0f);
    style.Colors[ImGuiCol_TitleBg] = ImVec4(0.08f, 0.06f, 0.14f, 1.0f);
    style.Colors[ImGuiCol_TitleBgActive] = ImVec4(0.16f, 0.10f, 0.24f, 1.0f);
    style.Colors[ImGuiCol_Tab] = ImVec4(0.10f, 0.10f, 0.18f, 1.0f);
    style.Colors[ImGuiCol_TabSelected] = ImVec4(0.16f, 0.10f, 0.24f, 1.0f);
    style.Colors[ImGuiCol_FrameBg] = ImVec4(0.04f, 0.04f, 0.08f, 1.0f);
    style.Colors[ImGuiCol_Button] = ImVec4(0.88f, 0.25f, 0.50f, 1.0f);
    style.Colors[ImGuiCol_ButtonHovered] = ImVec4(1.0f, 0.31f, 0.56f, 1.0f);
    style.Colors[ImGuiCol_PlotHistogram] = ImVec4(0.88f, 0.25f, 0.50f, 1.0f);

    ImGui_ImplSDL2_InitForOpenGL(window, gl_ctx);
    ImGui_ImplOpenGL3_Init("#version 330");

    GuiState state;
    state.window = window;
    state.gl_ctx = gl_ctx;

    // Load default scene
    load_file(state, "scenes/nanotekt.yaml");

    bool running = true;
    while (running) {
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) {
            ImGui_ImplSDL2_ProcessEvent(&ev);
            if (ev.type == SDL_QUIT)
                running = false;
            if (ev.type == SDL_DROPFILE) {
                load_file(state, ev.drop.file);
                std::strncpy(state.yaml_buf, state.yaml_content.c_str(), GuiState::BUF_SZ - 1);
                // Also update path_buf for display (via load_file which sets yaml_path)
                SDL_free(ev.drop.file);
            }
        }

        // File watching
        check_file(state);
        if (state.file_changed) {
            init_scene_for_wireframe(state);
            trigger_render(state);
        }

        // Check render completion
        if (state.rendering && is_render_done()) {
            state.rendering = false;
            if (get_render_result() == 0) {
                state.render_valid = true;
                state.render_w = get_width();
                state.render_h = get_height();
                // (Re)create textures
                if (state.tex_render) glDeleteTextures(1, &state.tex_render);
                if (state.tex_depth) glDeleteTextures(1, &state.tex_depth);
                state.tex_render = make_texture(state.render_w, state.render_h);
                state.tex_depth = make_texture(state.render_w, state.render_h);
                upload_rgb(state.tex_render, state.render_w, state.render_h, get_pixels_ptr());
                upload_depth(state.tex_depth, state.render_w, state.render_h, get_zbuf_ptr());
            }
        }

        // Begin ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplSDL2_NewFrame();
        ImGui::NewFrame();

        // Full-window dockspace
        ImGuiID dockspace_id = ImGui::DockSpaceOverViewport(0, ImGui::GetMainViewport());

        // Setup default layout on first frame
        if (state.first_frame) {
            setup_default_layout(dockspace_id);
            state.first_frame = false;
        }

        // Panels
        panel_yaml(state);
        panel_image("Render", state.tex_render, state.render_w, state.render_h,
                    state.rv_zoom, state.rv_pan_x, state.rv_pan_y);
        panel_image("Depth", state.tex_depth, state.render_w, state.render_h,
                    state.rv_zoom, state.rv_pan_x, state.rv_pan_y);
        panel_wireframe("Top View", state.top_zoom, state.top_pan_x, state.top_pan_y,
                        state.scene_loaded, true);
        panel_wireframe("Front View", state.front_zoom, state.front_pan_x, state.front_pan_y,
                        state.scene_loaded, false);

        // Render
        int disp_w, disp_h;
        SDL_GL_GetDrawableSize(window, &disp_w, &disp_h);
        glViewport(0, 0, disp_w, disp_h);
        glClearColor(0.06f, 0.06f, 0.10f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        SDL_GL_SwapWindow(window);
    }

    // Cleanup
    if (state.rendering) {
        cancel_render();
        while (!is_render_done()) SDL_Delay(5);
    }
    if (state.tex_render) glDeleteTextures(1, &state.tex_render);
    if (state.tex_depth) glDeleteTextures(1, &state.tex_depth);
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplSDL2_Shutdown();
    ImGui::DestroyContext();
    SDL_GL_DeleteContext(gl_ctx);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
