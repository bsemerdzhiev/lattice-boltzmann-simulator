#include "render_simulation.hpp"
#include "lbm.hpp"
#include "lbm_c.hpp"

#include "glad/glad.h"
#include <GLFW/glfw3.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

// Called whenever the window is resized; keeps the OpenGL viewport
// matching the actual framebuffer size.
void framebuffer_size_callback(GLFWwindow *window, int width, int height) {
  glViewport(0, 0, width, height);
}

const char *vertex_shader_src = R"(
#version 330 core
layout (location = 0) in vec2 aPos;
layout (location = 1) in vec2 aTexCoord;

out vec2 TexCoord;

void main() {
    gl_Position = vec4(aPos, 0.0, 1.0);
    TexCoord = aTexCoord;
}
)";

const char *fragment_shader_src = R"(
#version 330 core
in vec2 TexCoord;
out vec4 FragColor;

uniform sampler2D gridTexture;

void main() {
    FragColor = texture(gridTexture, vec2(TexCoord.x, 1 - TexCoord.y));
}
)";

GLuint compile_shader(GLenum type, const char *source) {
  GLuint shader = glCreateShader(type);
  glShaderSource(shader, 1, &source, nullptr);
  glCompileShader(shader);

  int success;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
  if (!success) {
    char log[512];
    glGetShaderInfoLog(shader, 512, nullptr, log);
    std::fprintf(stderr, "Shader compile error: %s\n", log);
  }
  return shader;
}

GLuint build_shader_program() {
  GLuint vertex = compile_shader(GL_VERTEX_SHADER, vertex_shader_src);
  GLuint fragment = compile_shader(GL_FRAGMENT_SHADER, fragment_shader_src);

  GLuint program = glCreateProgram();
  glAttachShader(program, vertex);
  glAttachShader(program, fragment);
  glLinkProgram(program);

  int success;
  glGetProgramiv(program, GL_LINK_STATUS, &success);
  if (!success) {
    char log[512];
    glGetProgramInfoLog(program, 512, nullptr, log);
    std::fprintf(stderr, "Shader link error: %s\n", log);
  }

  glDeleteShader(vertex);
  glDeleteShader(fragment);
  return program;
}

std::vector<unsigned char> compute_grid_colors() {
  std::vector<unsigned char> pixels(LBM_CONSTANTS::WIDTH *
                                    LBM_CONSTANTS::HEIGHT * 3);

  float max_speed = 1e-6f; // avoid divide-by-zero on the first frame
  for (std::size_t i = 0; i < LBM_CONSTANTS::HEIGHT; i++) {
    for (std::size_t j = 0; j < LBM_CONSTANTS::WIDTH; j++) {
      if (j == 0 && (i == 0 || i == LBM_CONSTANTS::HEIGHT - 1)) {
        continue;
      }
      const auto &v = Cell::velocity[i][j];
      float speed = std::sqrt(v.dot_product(v));
      max_speed = std::max(max_speed, speed);
    }
  }

  for (std::size_t i = 0; i < LBM_CONSTANTS::HEIGHT; i++) {
    for (std::size_t j = 0; j < LBM_CONSTANTS::WIDTH; j++) {
      const auto &v = Cell::velocity[i][j];
      float speed = std::sqrt(v.dot_product(v));
      float t = speed / max_speed; // normalized 0..1
      //
      // std::cout << t << " ";

      std::size_t idx = (i * LBM_CONSTANTS::WIDTH + j) * 3;
      if (Cell::blockade[i][j]) {
        pixels[idx + 0] = 254; // static_cast<unsigned char>(255.f * t); // R
        pixels[idx + 1] = 254; // G
        pixels[idx + 2] =
            0; // static_cast<unsigned char>(255.f * (1.f - t)); // B
      } else {
        pixels[idx + 0] = static_cast<unsigned char>(255.f * t);         // R
        pixels[idx + 1] = 0;                                             // G
        pixels[idx + 2] = static_cast<unsigned char>(255.f * (1.f - t)); // B
      }
    }
    // std::cout << "\n";
  }

  return pixels;
}

void RenderSimulation::render() {
  if (!glfwInit()) {
    return;
  }
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  const int WINDOW_WIDTH = 1000;
  const int WINDOW_HEIGHT = 450;
  GLFWwindow *window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT,
                                        "LBM Visualizer", nullptr, nullptr);
  if (window == nullptr) {
    glfwTerminate();
    return;
  }
  glfwMakeContextCurrent(window);
  glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

  // --- Initialize GLAD (loads OpenGL function pointers from the driver) ---
  if (!gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress))) {
    glfwTerminate();
    return;
  }
  glViewport(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT);

  GLuint shader_program = build_shader_program();

  // Full-window quad: 4 vertices (position.xy, texcoord.xy), two
  // triangles via an element buffer. You shouldn't need to touch this.
  float vertices[] = {
      -1.0f, -1.0f, 0.0f, 0.0f, // bottom-left
      1.0f,  -1.0f, 1.0f, 0.0f, // bottom-right
      1.0f,  1.0f,  1.0f, 1.0f, // top-right
      -1.0f, 1.0f,  0.0f, 1.0f, // top-left
  };
  unsigned int indices[] = {0, 1, 2, 2, 3, 0};

  GLuint vao, vbo, ebo;
  glGenVertexArrays(1, &vao);
  glGenBuffers(1, &vbo);
  glGenBuffers(1, &ebo);

  glBindVertexArray(vao);

  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices,
               GL_STATIC_DRAW);

  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void *)0);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float),
                        (void *)(2 * sizeof(float)));
  glEnableVertexAttribArray(1);

  GLuint texture;

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glGenTextures(1, &texture);
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, LBM_CONSTANTS::WIDTH,
               LBM_CONSTANTS::HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);

  glUseProgram(shader_program);
  glUniform1i(glGetUniformLocation(shader_program, "gridTexture"), 0);

  LBM::initialize();

  bool k = 0;
  int32_t cnt = 0;
  // --- Render loop ---
  while (!glfwWindowShouldClose(window)) {
    if (glfwGetKey(window, GLFW_KEY_CAPS_LOCK) == GLFW_PRESS) {
      glfwSetWindowShouldClose(window, true);
    }

    LBM::update(k);
    if (cnt % 100 == 0) {

      std::cout << cnt << "\n";
    }
    cnt++;
    k ^= 1;

    std::vector<unsigned char> pixels = compute_grid_colors();
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, LBM_CONSTANTS::WIDTH,
                    LBM_CONSTANTS::HEIGHT, GL_RGB, GL_UNSIGNED_BYTE,
                    pixels.data());

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glActiveTexture(GL_TEXTURE0);
    glUseProgram(shader_program);
    glBindVertexArray(vao);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    glfwSwapBuffers(window);
    glfwPollEvents();
  }

  glfwTerminate();
  return;
}
