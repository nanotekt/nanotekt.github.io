# nanore

A CPU-based software raytracer written in C++17 with no external dependencies. Renders scenes with reflective metals, refractive glass, volumetric fog, bloom post-processing, and procedural geometry — all driven by YAML scene configuration.

## Prerequisites

- **g++** with C++17 support
- **POSIX threads** (standard on Linux/macOS)
- No external libraries required

## Rendering

```bash
./render.sh scenes/qbd.yaml
```

This builds the project (if needed) and renders the scene. Output is saved to `renders/render-<timestamp>.bmp`.

## Scene Configuration

Scenes are defined in YAML files (see [scenes/qbd.yaml](scenes/qbd.yaml) for an example). Configurable parameters include:

- **Resolution and sampling** — image dimensions, samples per pixel
- **Camera** — position and target
- **Objects** — spheres, rounded cubes, discs, floor tiles, glass cubes
- **Materials** — chrome, silver, glass (with Fresnel/refraction), diffuse, emissive
- **Lighting** — multiple point lights with color
- **Post-processing** — bloom radius, gamma correction
- **Volumetric fog** — Perlin noise-based density and distance
