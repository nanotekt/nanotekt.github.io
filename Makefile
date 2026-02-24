CXX := g++
CXXFLAGS := -O3 -std=c++17 -pthread
LDFLAGS := -lm

nanore: src/nanore.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

# --- GUI build (SDL2 + GLEW + ImGui) ---
GUI_CXXFLAGS := $(CXXFLAGS) -DNANORE_GUI $(shell sdl2-config --cflags) -Ilib/imgui -Ilib/imgui/backends
GUI_LDFLAGS  := $(LDFLAGS) $(shell sdl2-config --libs) -lGLEW -lGL

IMGUI_SRC := lib/imgui/imgui.cpp lib/imgui/imgui_draw.cpp lib/imgui/imgui_tables.cpp \
             lib/imgui/imgui_widgets.cpp lib/imgui/backends/imgui_impl_sdl2.cpp \
             lib/imgui/backends/imgui_impl_opengl3.cpp

nanore-gui: src/nanore.cpp src/gui.cpp $(IMGUI_SRC)
	$(CXX) $(GUI_CXXFLAGS) -o $@ $^ $(GUI_LDFLAGS)

clean:
	rm -f nanore nanore-gui

.PHONY: clean
