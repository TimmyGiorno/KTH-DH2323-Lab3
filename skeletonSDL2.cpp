// DH2323 skeleton code, Lab3 (SDL2 version)
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <glm/glm.hpp>
#include <glm/gtx/constants.hpp>
#include "SDL2auxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::vec2;
using glm::mat3;
using glm::ivec2;

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES
const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL2Aux *sdlAux;
int t;
vector<Triangle> triangles;

vec3 cameraPos(0, 0, -3.001); // Fixed camera position (from original)
float focalLength = SCREEN_HEIGHT;
float yaw = 0.0f; // Horizontal angle (around Y-axis)
float pitch = 0.0f; // Vertical angle (around X-axis)
float roll = 0.0f;
float rotationSpeed = 0.002f; // Radians per ms
float moveSpeed = 0.005f;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
vec3 lightPos(0.0f, -0.5f, -0.7f);
vec3 lightPower = 1.5f * vec3(1.0f, 1.0f, 1.0f);
vec3 indirectLightPowerPerArea = 0.3f * vec3(1.0f, 1.0f, 1.0f);
vec3 currentNormal;
vec3 currentReflectance;

struct Vertex
{
    vec3 position;
    vec3 normal;
    vec3 reflectance;
};

struct Pixel
{
    int x;
    int y;
    float zinv;
    vec3 illumination;
};

// ----------------------------------------------------------------------------
// FUNCTIONS Prototypes
void Update();
void Draw();
void PixelShader(const Pixel& p); // Task 7.1
void VertexShader(const Vertex& v, Pixel& p); // Task 7.2 / 7.5
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result);
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result); // Kept for DrawLineSDL
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
void DrawPolygon(const vector<Vertex>& vertices); // Task 7.3
void DrawLineSDL(ivec2 a, ivec2 b, vec3 color);
void DrawPolygonEdges(const vector<Vertex>& vertices);


// ----------------------------------------------------------------------------
// FUNCTIONS
int SDL_main(int argc, char* argv[])
{
    LoadTestModel(triangles); // Load model
    sdlAux = new SDL2Aux(SCREEN_WIDTH, SCREEN_HEIGHT);
    t = SDL_GetTicks(); // Set start value for timer.

    while (!sdlAux->quitEvent())
    {
        Update();
        Draw();
    }

    sdlAux->saveBMP("screenshot_illumination_dim.bmp"); // Changed filename
    delete sdlAux;
    return 0;
}


void Update()
{
    int t2 = SDL_GetTicks();
    float dt = float(t2 - t);
    t = t2;
    cout << "Render time: " << dt << " ms." << endl;

    const Uint8* keystate = SDL_GetKeyboardState(NULL);
    float move = moveSpeed * dt;
    float rot = rotationSpeed * dt;

    vec3 forward(sin(yaw), 0, cos(yaw));
    vec3 right(cos(yaw), 0, -sin(yaw));

    if (keystate[SDL_SCANCODE_UP])    cameraPos += move * forward;
    if (keystate[SDL_SCANCODE_DOWN])  cameraPos -= move * forward;
    if (keystate[SDL_SCANCODE_LEFT])  cameraPos -= move * right;
    if (keystate[SDL_SCANCODE_RIGHT]) cameraPos += move * right;

    if (keystate[SDL_SCANCODE_W]) pitch -= rot;
    if (keystate[SDL_SCANCODE_S]) pitch += rot;
    if (keystate[SDL_SCANCODE_A]) roll  += rot;
    if (keystate[SDL_SCANCODE_D]) roll  -= rot;
    if (keystate[SDL_SCANCODE_Q]) yaw   -= rot;
    if (keystate[SDL_SCANCODE_E]) yaw   += rot;

    // Clamp pitch (Matches original skeleton style)
    float limit = glm::pi<float>() / 2.0f - 0.01f; // Approx 89.9 degrees
    pitch = std::max(-limit, std::min(limit, pitch));
    // roll = 0; // Optional: uncomment to disable roll
}


void Draw()
{
    sdlAux->clearPixels();

    for(int y = 0; y < SCREEN_HEIGHT; ++y) {
        for(int x = 0; x < SCREEN_WIDTH; ++x) {
            depthBuffer[y][x] = 0;
        }
    }

    for (size_t i = 0; i < triangles.size(); ++i)
    {
        vector<Vertex> vertices(3); // Changed from vec3

        vertices[0].position = triangles[i].v0;
        vertices[0].normal = triangles[i].normal;
        vertices[0].reflectance = triangles[i].color;

        vertices[1].position = triangles[i].v1;
        vertices[1].normal = triangles[i].normal;
        vertices[1].reflectance = triangles[i].color;

        vertices[2].position = triangles[i].v2;
        vertices[2].normal = triangles[i].normal;
        vertices[2].reflectance = triangles[i].color;

        DrawPolygon(vertices);
        // DrawPolygonEdges(vertices); // Optional edge drawing
    }

    sdlAux->render();
}


void VertexShader(const Vertex& v, Pixel& p)
{
    // Illumination
    vec3 worldPos = v.position;
    vec3 worldNormal = v.normal;
    vec3 reflectance = v.reflectance;
    vec3 r_vec = lightPos - worldPos;
    float r_dist_sq = glm::dot(r_vec, r_vec);
    float r_dist = glm::sqrt(r_dist_sq);

    if (r_dist > 1e-5) {
        vec3 r_hat = r_vec / r_dist;
        vec3 n_hat = glm::normalize(worldNormal);
        float diffuseFactor = std::max(0.0f, glm::dot(r_hat, n_hat));
        // Using 1/r^2 falloff
        float denom = r_dist_sq;
        vec3 D = lightPower * diffuseFactor / std::max(denom, 1e-5f);
        p.illumination = reflectance * (D + indirectLightPowerPerArea);
    } else {
         p.illumination = reflectance * indirectLightPowerPerArea;
    }
    p.illumination = glm::clamp(p.illumination, 0.0f, 1.0f); // Clamp color

    vec3 relative = worldPos - cameraPos;

    mat3 Rx = mat3(1, 0, 0,
                   0, cos(pitch), -sin(pitch),
                   0, sin(pitch), cos(pitch));
    mat3 Ry = mat3(cos(yaw), 0, sin(yaw),
                   0, 1, 0,
                   -sin(yaw), 0, cos(yaw));
    mat3 Rz = mat3(cos(roll), -sin(roll), 0,
                   sin(roll), cos(roll), 0,
                   0, 0, 1);

    mat3 R = Ry * Rx * Rz; // Yaw -> Pitch -> Roll relative to camera

    vec3 rotated = R * relative;

    if (rotated.z <= 0) {
        p.x = -1; p.y = -1; p.zinv = 0; return; // Mark as invalid
    }

    p.zinv = 1.0f / rotated.z;
    p.x = static_cast<int>(focalLength * rotated.x * p.zinv + SCREEN_WIDTH / 2.0f);
    p.y = static_cast<int>(focalLength * rotated.y * p.zinv + SCREEN_HEIGHT / 2.0f);
}


void PixelShader( const Pixel& p )
{
    int x = p.x;
    int y = p.y;

    // Bounds check
    if (x < 0 || x >= SCREEN_WIDTH || y < 0 || y >= SCREEN_HEIGHT) {
        return;
    }

    // Depth Test
    if( p.zinv > depthBuffer[y][x] )
    {
        depthBuffer[y][x] = p.zinv; // Update depth buffer
        sdlAux->putPixel( x, y, p.illumination ); // Use interpolated illumination
    }
}


void Interpolate(Pixel a, Pixel b, vector<Pixel>& result)
{
    size_t N = result.size();
    if (N == 0) return;
    if (N == 1) { result[0] = a; return; }

    float steps = static_cast<float>(N - 1);
    float stepX = static_cast<float>(b.x - a.x) / steps;
    float stepY = static_cast<float>(b.y - a.y) / steps;
    float stepZinv = (b.zinv - a.zinv) / steps;
    vec3 stepIllum = (b.illumination - a.illumination) / steps;

    float currentX = static_cast<float>(a.x);
    float currentY = static_cast<float>(a.y);
    float currentZinv = a.zinv;
    vec3 currentIllum = a.illumination;

    for (size_t i = 0; i < N; ++i)
    {
        result[i].x = static_cast<int>(round(currentX));
        result[i].y = static_cast<int>(round(currentY));
        result[i].zinv = currentZinv;
        result[i].illumination = currentIllum;

        currentX += stepX;
        currentY += stepY;
        currentZinv += stepZinv;
        currentIllum += stepIllum;
    }
}

// Original Interpolate function kept for backward compatibility with DrawLineSDL
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result)
{
    size_t N = result.size();
    if (N == 0) return;
    if (N == 1) { result[0] = a; return; }

    float steps = static_cast<float>(N - 1);
    vec2 step = vec2(b - a) / steps;
    vec2 current(a);

    for (size_t i = 0; i < N; ++i)
    {
        result[i] = ivec2(round(current.x), round(current.y));
        current += step;
    }
}


void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels)
{
    // 1. Find min and max y-value of the polygon
    int minY = SCREEN_HEIGHT, maxY = -1;
    bool validVertexFound = false;
    for (const auto& vp : vertexPixels) {
        if (vp.x != -1) { // Check if vertex is valid
             minY = min(minY, vp.y);
             maxY = max(maxY, vp.y);
             validVertexFound = true;
        }
    }
    if (!validVertexFound || maxY < minY) { leftPixels.clear(); rightPixels.clear(); return; }

    // Clamp Y to screen bounds
    minY = max(0, minY);
    maxY = min(SCREEN_HEIGHT - 1, maxY);

    // 2. Resize leftPixels and rightPixels
    int numRows = maxY - minY + 1;
    if (numRows <= 0) { leftPixels.clear(); rightPixels.clear(); return; }
    leftPixels.resize(numRows);
    rightPixels.resize(numRows);

    // 3. Initialize the x-coordinates and other attributes
    for (int i = 0; i < numRows; ++i)
    {
        leftPixels[i].x = numeric_limits<int>::max();
        leftPixels[i].y = minY + i;
        leftPixels[i].zinv = 0.0f; // Initialize defaults

        rightPixels[i].x = numeric_limits<int>::min();
        rightPixels[i].y = minY + i;
        rightPixels[i].zinv = 0.0f; // Initialize defaults
    }

    // 4. Loop through all edges and interpolate x-coordinates and other attributes
    int V = vertexPixels.size();
    for (int i = 0; i < V; ++i)
    {
        Pixel a = vertexPixels[i];
        Pixel b = vertexPixels[(i + 1) % V];

        if(a.x == -1 || b.x == -1) continue; // Skip edges with invalid vertices

        int dy = abs(a.y - b.y);
        int dx = abs(a.x - b.x);
        int steps = max(dx, dy) + 1;
        steps = max(steps, 1); // Ensure at least one step

        vector<Pixel> line(steps);
        Interpolate(a, b, line); // Interpolate all Pixel attributes

        for (const auto& p : line)
        {
            if (p.y >= minY && p.y <= maxY)
            {
                int row = p.y - minY;
                if (p.x < leftPixels[row].x)
                {
                    leftPixels[row].x = p.x;
                    leftPixels[row].zinv = p.zinv;          // Store interpolated value
                    leftPixels[row].illumination = p.illumination; // Store interpolated value
                }
                if (p.x > rightPixels[row].x)
                {
                    rightPixels[row].x = p.x;
                    rightPixels[row].zinv = p.zinv;          // Store interpolated value
                    rightPixels[row].illumination = p.illumination; // Store interpolated value
                }
            }
        }
    }
}


void DrawRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels)
{
    if (leftPixels.size() != rightPixels.size() || leftPixels.empty()) return;

    for (size_t i = 0; i < leftPixels.size(); ++i)
    {
        Pixel left = leftPixels[i];
        Pixel right = rightPixels[i];

        // Skip invalid rows
        if (left.x > right.x) continue;

        // Original width for interpolation base
        int originalWidth = right.x - left.x + 1;
        if (originalWidth <= 0) continue;

        vector<Pixel> rowPixels(originalWidth);
        Interpolate(left, right, rowPixels); // Interpolate across the original width

        for (const auto& p : rowPixels)
        {
             PixelShader(p);
        }
    }
}

// Task 7.3: Modified DrawPolygon
void DrawPolygon(const vector<Vertex>& vertices)
{
    int V = vertices.size();
    if (V < 3) return;

    vector<Pixel> vertexPixels(V);
    bool anyVertexValid = false;
    for (int i = 0; i < V; ++i)
    {
        VertexShader(vertices[i], vertexPixels[i]);
        if(vertexPixels[i].x != -1) anyVertexValid = true;
    }

    // Simple full polygon clip if all vertices are invalid
    if (!anyVertexValid) return;

    vector<Pixel> leftPixels;
    vector<Pixel> rightPixels;
    ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
    DrawRows(leftPixels, rightPixels);
}


// ----------------------------------------------------------------------------
// UTILITY FUNCTIONS

void DrawLineSDL(ivec2 a, ivec2 b, vec3 color)
{
    ivec2 delta = glm::abs(a - b);
    int pixels = glm::max(delta.x, delta.y) + 1;

    vector<ivec2> line(pixels);
    Interpolate(a, b, line); // Use ivec2 version

    for (const auto& point : line)
    {
        if (point.x >= 0 && point.x < SCREEN_WIDTH && point.y >= 0 && point.y < SCREEN_HEIGHT)
        {
            sdlAux->putPixel(point.x, point.y, color);
        }
    }
}

void DrawPolygonEdges(const vector<Vertex>& vertices)
{
    int V = vertices.size();
    if (V < 2) return;

    vector<Pixel> projectedVertices(V);
    for (int i = 0; i < V; ++i)
    {
        VertexShader(vertices[i], projectedVertices[i]);
    }

    vec3 white(1.0f, 1.0f, 1.0f);
    for (int i = 0; i < V; ++i)
    {
        Pixel p1 = projectedVertices[i];
        Pixel p2 = projectedVertices[(i + 1) % V];
        // Only draw edges where both endpoints are valid
        if (p1.x != -1 && p2.x != -1) {
            DrawLineSDL(ivec2(p1.x, p1.y), ivec2(p2.x, p2.y), white);
        }
    }
}