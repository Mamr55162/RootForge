#include <iostream>
#include <GLFW/glfw3.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

// Optional error callback for GLFW
static void glfw_error_callback(int error, const char* description) {
    std::cerr << "GLFW Error " << error << ": " << description << std::endl;
}

int main() {
    // 1. Initialize GLFW and create a window
    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit()) return 1;

    const char* glsl_version = "#version 130";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

    GLFWwindow* window = glfwCreateWindow(800, 600, "Numerical Analysis Engine", nullptr, nullptr);
    if (window == nullptr) return 1;
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync

    // 2. Initialize Dear ImGui Context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    
    ImGui::StyleColorsDark(); // Apply dark theme

    // 3. Initialize Platform/Renderer Backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

    // Application State Variables
    float x_input = 0.0f;
    float polynomial_result = 0.0f;

    // 4. The Main Application Loop
    while (!glfwWindowShouldClose(window)) {
        // Poll hardware events (keys, mouse)
        glfwPollEvents();

        // Start the ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // --- YOUR UI LOGIC STARTS HERE ---
        ImGui::Begin("Polynomial Evaluator"); // Create a window

        ImGui::Text("Function: f(x) = 3x^2 + 2x - 5");
        ImGui::Spacing();

        // Input widget tied directly to x_input
        ImGui::InputFloat("Enter x value", &x_input);

        // Button widget executing logic
        if (ImGui::Button("Calculate f(x)")) {
            polynomial_result = (3.0f * x_input * x_input) + (2.0f * x_input) - 5.0f;
        }

        ImGui::Separator();
        
        // Display the result
        ImGui::TextColored(ImVec4(0.0f, 1.0f, 0.0f, 1.0f), "Result: %.4f", polynomial_result);

        ImGui::End();
        // --- YOUR UI LOGIC ENDS HERE ---

        // 5. Render to the Screen
        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    // 6. Cleanup and Shutdown
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}