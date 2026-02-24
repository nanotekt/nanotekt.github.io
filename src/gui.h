#pragma once

// GUI entry point — called from main() when no scene file is specified
int gui_main(int argc, char* argv[]);

// Async render API (defined in nanore.cpp)
extern "C" {
int start_render(const char* yaml_str);
void cancel_render();
int is_render_done();
int get_render_result();
int get_progress();
unsigned char* get_pixels_ptr();
float* get_zbuf_ptr();
int get_width();
int get_height();
}

// Geometry accessors for wireframe views (defined in nanore.cpp)
extern "C" {
int get_sphere_count();
void get_sphere_center(int i, double* x, double* y, double* z);
double get_sphere_radius();
int get_rcube_count();
void get_rcube_info(int i, double* cx, double* cy, double* cz, double* half, double* round_);
void get_disc_info(double* cx, double* cy, double* cz, double* r);
void get_floor_bounds(int* lox, int* hix, int* loz, int* hiz);
void get_camera_info(double* px, double* py, double* pz, double* tx, double* ty, double* tz);
int get_light_count();
void get_light_info(int i, double* x, double* y, double* z);
int get_standalone_sphere_count();
void get_standalone_sphere(int i, double* cx, double* cy, double* cz, double* gr, double* cr);
int get_glass_sub_count();
void get_glass_sub(int i, double* x0, double* y0, double* z0, double* x1, double* y1, double* z1);
// Initialize scene from YAML (for wireframe preview without rendering)
int init_scene_from_yaml_str(const char* yaml_str);
}
