# nanore — AI Workflow

## Project Overview

CPU-based C++17 raytracer. Single source file `src/nanore.cpp`, YAML-driven scene configuration.

## Build & Render Workflow

Render a scene (builds automatically):
```bash
./render.sh scenes/<scene>.yaml
```
Output is saved to `renders/render-<timestamp>.bmp`.

To build without rendering:
```bash
bash build.sh
```

## Creating New Scenes

- Create YAML files in `scenes/`. Use `scenes/qbd.yaml` as the reference for all available parameters.
- Scene files control everything without recompilation: resolution, camera, geometry, materials, lighting, fog, bloom.
- After creating or editing a scene file, run `./render.sh scenes/<scene>.yaml` to render.

### Quick iteration tips

- Use low `samples: 1` and small resolution (`width: 256`, `height: 256`) for fast preview renders.
- Increase to `samples: 4`+ and full resolution for final output.

## Key Parameters (all optional, have defaults)

| Section | What it controls |
|---------|-----------------|
| `rendering` | Resolution, sample count, fog distance |
| `camera` | Position and look-at target |
| `text` | Bitmap-encoded text spheres (`@` = filled, `.` = empty) |
| `rounded_cube` | Reflective rounded cube (center, size, rotation) |
| `red_disc` | Emissive disc light |
| `floor` | Tiled floor grid (bounds, colors, reflectivity) |
| `glass` | Glass cube properties (IOR, Fresnel, tint) |
| `lights` | Point light positions and colors |
| `materials` | Chrome and silver material properties |
| `scratches` | Procedural scratch layers on silver |
| `sky` | Neon grid ceiling |
| `fog` | Volumetric Perlin noise fog |
| `bloom` | Glow post-processing and gamma |

## File Layout

- `src/nanore.cpp` — entire renderer (single file)
- `scenes/*.yaml` — scene definitions (must start with `version: nanore@1`)
- `renders/*.bmp` — rendered output (timestamped, not tracked by git)
- `render.sh` — build + render script
- `build.sh` — build-only script
