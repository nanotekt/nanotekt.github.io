#!/bin/bash
set -e

mkdir -p dist

echo "Building nanore for WebAssembly..."

em++ src/nanore.cpp \
  -O2 \
  -std=c++17 \
  -pthread \
  -sALLOW_MEMORY_GROWTH=1 \
  -sINITIAL_MEMORY=134217728 \
  -sPTHREAD_POOL_SIZE='navigator.hardwareConcurrency+1' \
  -sEXPORTED_FUNCTIONS='["_main","_start_render","_cancel_render","_is_render_done","_get_render_result","_get_progress","_get_pixels_ptr","_get_zbuf_ptr","_get_width","_get_height","_malloc","_free"]' \
  -sEXPORTED_RUNTIME_METHODS='["ccall","cwrap","stringToUTF8","lengthBytesUTF8","wasmMemory"]' \
  -sENVIRONMENT=web,worker \
  -o dist/nanore.js

cp web/coi-serviceworker.js dist/

echo "Embedding default scene YAML into index.html..."

python3 -c "
with open('scenes/nanotekt.yaml', 'r') as f:
    yaml = f.read()
with open('web/index.html', 'r') as f:
    html = f.read()
import re
yaml = re.sub(r'^  width: \d+', '  width: 256', yaml, count=1, flags=re.MULTILINE)
yaml = re.sub(r'^  height: \d+', '  height: 256', yaml, count=1, flags=re.MULTILINE)
yaml = re.sub(r'^  samples: \d+', '  samples: 1', yaml, count=1, flags=re.MULTILINE)
html = html.replace('__DEFAULT_YAML__', yaml.rstrip())
with open('dist/index.html', 'w') as f:
    f.write(html)
"

# Create a minimal dev server with required COOP/COEP headers
cat > dist/serve.py << 'SERVE'
#!/usr/bin/env python3
import http.server
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

class Handler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Cross-Origin-Opener-Policy', 'same-origin')
        self.send_header('Cross-Origin-Embedder-Policy', 'require-corp')
        super().end_headers()

print('Serving nanore at http://localhost:8080')
http.server.HTTPServer(('', 8080), Handler).serve_forever()
SERVE

echo ""
echo "Build complete. Output in dist/"
echo "  dist/nanore.js"
echo "  dist/nanore.wasm"
echo "  dist/nanore.worker.js"
echo "  dist/index.html"
echo "  dist/serve.py"
echo ""
echo "To run:  python3 dist/serve.py"
