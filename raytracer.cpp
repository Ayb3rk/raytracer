#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <cfloat>

typedef unsigned char RGB[3];

struct Ray {
    parser::Vec3f origin;
    parser::Vec3f direction;
};

parser::Vec3f multS(parser::Vec3f a,float s)
{
    parser::Vec3f result{};
    result.x = a.x*s;
    result.y = a.y*s;
    result.z = a.z*s;
    return result;
}

float dot(parser::Vec3f a,parser::Vec3f b)
{
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

parser::Vec3f normalize(parser::Vec3f a)
{
    return multS(a,1.0/sqrt(dot(a,a)));
}

parser::Vec3f hadamard(parser::Vec3f a,parser::Vec3f b)
{
    parser::Vec3f result{};
    result.x = a.x*b.x;
    result.y = a.y*b.y;
    result.z = a.z*b.z;
    return result;
}

parser::Vec3f add(parser::Vec3f a, parser::Vec3f b)
{
    parser::Vec3f result{};
    result.x = a.x+b.x;
    result.y = a.y+b.y;
    result.z = a.z+b.z;
    return result;
}

parser::Vec3f sub(parser::Vec3f a, parser::Vec3f b)
{
    parser::Vec3f result{};
    result.x = a.x-b.x;
    result.y = a.y-b.y;
    result.z = a.z-b.z;
    return result;
}
parser::Vec3f CrossProduct(parser::Vec3f a, parser::Vec3f b) {
    parser::Vec3f result{};
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

parser::Vec3f NormalOfFace(parser::Face face, std::vector<parser::Vec3f> vertices) {
    parser::Vec3f v0 = vertices[face.v0_id];
    parser::Vec3f v1 = vertices[face.v1_id];
    parser::Vec3f v2 = vertices[face.v2_id];
    parser::Vec3f v01 = sub(v1, v0);
    parser::Vec3f v02 = sub(v2, v0);
    return normalize(CrossProduct(v01, v02));
}

Ray GenerateRay(int i, int j, const parser::Camera& camera, int width, int height)
{
    Ray result{};
    float su,sv;
    parser::Vec3f m{},q{},s{};

    auto left = camera.near_plane.x;
    auto right = camera.near_plane.y;
    auto bottom = camera.near_plane.z;
    auto top = camera.near_plane.w;

    su = (float) (i+0.5)*(right - left)/(float) width;
    sv = (float) (j+0.5)*(top - bottom)/(float) height;
    parser::Vec3f v = camera.up;
    parser::Vec3f w = multS(camera.gaze,-1);
    parser::Vec3f u = CrossProduct(v,w);

    m = add(camera.position,multS(camera.gaze,camera.near_distance));
    q = add(m, add(multS(u,left),multS(v,top)));
    s = add(q,add(multS(u,su),multS(v,-sv)));
    result.origin = camera.position;
    result.direction = add(s,multS(camera.position,-1));
    return result;
}



float intersectSphere (parser::Vec3f sphereCenter, float sphereRadius, Ray ray) {
    float A,B,C; //constants for the quadratic equation
    float delta;
    float t,t1,t2;

    parser::Vec3f c = sphereCenter;

    C = (ray.origin.x-c.x)*(ray.origin.x-c.x)+(ray.origin.y-c.y)*(ray.origin.y-c.y)+(ray.origin.z-c.z)*(ray.origin.z-c.z)-sphereRadius*sphereRadius;

    B = 2*ray.direction.x*(ray.origin.x-c.x)+2*ray.direction.y*(ray.origin.y-c.y)+2*ray.direction.z*(ray.origin.z-c.z);

    A = ray.direction.x*ray.direction.x+ray.direction.y*ray.direction.y+ray.direction.z*ray.direction.z;

    delta = B*B-4*A*C;

    if (delta<0) return -1;
    else if (delta==0)
    {
        t = -B / (2*A);
    }
    else
    {
        delta = sqrt(delta);
        A = 2*A;
        t1 = (-B + delta) / A;
        t2 = (-B - delta) / A;

        t1<t2 ? t=t1 : t=t2;
    }

    return t;
}

float determinant(float matrix[3][3]) {
    float result = 0;
    result += matrix[0][0] * matrix[1][1] * matrix[2][2]; // a*e*i
    result += matrix[1][0] * matrix[2][1] * matrix[0][2]; // b*f*g
    result += matrix[2][0] * matrix[0][1] * matrix[1][2]; // c*d*h
    result -= matrix[2][0] * matrix[1][1] * matrix[0][2]; // c*e*g
    result -= matrix[1][0] * matrix[0][1] * matrix[2][2]; // b*d*i
    result -= matrix[0][0] * matrix[2][1] * matrix[1][2]; // a*f*h
    return result;
}

//intersect mesh function by using barycentric coordinates
float intersectFace (parser::Face face, Ray ray, parser::Scene scene) {
    parser::Vec3f v0 = scene.vertex_data[face.v0_id - 1];
    parser::Vec3f v1 = scene.vertex_data[face.v1_id - 1];
    parser::Vec3f v2 = scene.vertex_data[face.v2_id - 1];
    float beta, gamma, t;

    float AMatrix[3][3] = {
        {v0.x - v1.x, v0.x - v2.x, ray.direction.x}, //[a_x - b_x, a_x - c_x, d_x]
        {v0.y - v1.y, v0.y - v2.y, ray.direction.y}, //[a_y - b_y, a_y - c_y, d_y]
        {v0.z - v1.z, v0.z - v2.z, ray.direction.z} //[a_z - b_z, a_z - c_z, d_z]
    };

    float betaMatrix[3][3] = {
        {v0.x - ray.origin.x, v0.x - v2.x, ray.direction.x}, //[a_x - o_x, a_x - c_x, d_x]
        {v0.y - ray.origin.y, v0.y - v2.y, ray.direction.y}, //[a_y - o_y, a_y - c_y, d_y]
        {v0.z - ray.origin.z, v0.z - v2.z, ray.direction.z} //[a_z - o_z, a_z - c_z, d_z]
    };

    float gammaMatrix[3][3] = {
        {v0.x - v1.x, v0.x - ray.origin.x, ray.direction.x}, //[a_x - b_x, a_x - o_x, d_x]
        {v0.y - v1.y, v0.y - ray.origin.y, ray.direction.y}, //[a_y - b_y, a_y - o_y, d_y]
        {v0.z - v1.z, v0.z - ray.origin.z, ray.direction.z} //[a_z - b_z, a_z - o_z, d_z]
    };

    float tMatrix[3][3] = {
        {v0.x - v1.x, v0.x - v2.x, v0.x - ray.origin.x}, //[a_x - b_x, a_x - c_x, a_x - o_x]
        {v0.y - v1.y, v0.y - v2.y, v0.y - ray.origin.y}, //[a_y - b_y, a_y - c_y, a_y - o_y]
        {v0.z - v1.z, v0.z - v2.z, v0.z - ray.origin.z} //[a_z - b_z, a_z - c_z, a_z - o_z]
    };

    auto determinantA = determinant(AMatrix);
    if (determinantA == 0) return -1;

    beta = determinant(betaMatrix) / determinantA;
    gamma = determinant(gammaMatrix) / determinantA;
    t = determinant(tMatrix) / determinantA;

    if (beta >= 0 && gamma >= 0 && beta + gamma <= 1) return t;
    else return -1;

}

double CalculateLightDistance(parser::Vec3f lightPosition, parser::Vec3f intersectionPoint) {
    return sqrt(pow(lightPosition.x - intersectionPoint.x, 2) + pow(lightPosition.y - intersectionPoint.y, 2) + pow(lightPosition.z - intersectionPoint.z, 2));
}

parser::Vec3f ComputeColor(Ray ray, parser::Scene scene, int depth)
{

    if(depth == 0) {
        return {0, 0, 0};
    }
    float t = FLT_MAX;
    parser::Vec3f color{};
    parser::Sphere intersectedSphere{};
    parser::Mesh intersectedMesh{};
    for (const auto& sphere : scene.spheres) {
        auto sphereCenter = scene.vertex_data[sphere.center_vertex_id-1];
        float temp = intersectSphere(sphereCenter, sphere.radius, ray);
        if (temp > 0 && temp < t) {
            t = temp;
            intersectedSphere = sphere;
        }
    }
    for (const auto& mesh : scene.meshes) {
        for (auto face : mesh.faces) {
            auto temp = intersectFace(face, ray, scene);
            if (temp > 0 && temp < t) {
                t = temp;
                intersectedMesh = mesh;
                intersectedSphere = {};
            }
        }
    }
    if (t == FLT_MAX) {
        return color;
    }
    if (intersectedSphere.radius != 0) {
        auto P = add(add(ray.origin,multS(ray.direction,t)), {scene.shadow_ray_epsilon, scene.shadow_ray_epsilon, scene.shadow_ray_epsilon});
        auto N = normalize(add(P,multS(scene.vertex_data[intersectedSphere.center_vertex_id-1],-1)));
        auto W = normalize(add(ray.origin,multS(P,-1)));
        for (const auto& light : scene.point_lights) {
            auto L = normalize(add(light.position,multS(P,-1)));
            auto H = normalize(add(L,W));
            auto lightDistance = CalculateLightDistance(light.position, P);
            auto receivedIrradiance = multS(light.intensity,1/(lightDistance*lightDistance));
            auto cosAlphaSpecular = fmax(dot(N,H),0.0);
            auto cosAlphaDiffuse = fmax(dot(N,L),0.0);
            auto specular = multS(scene.materials[intersectedSphere.material_id-1].specular,pow(cosAlphaSpecular, scene.materials[intersectedSphere.material_id-1].phong_exponent));
            auto diffuse = multS(scene.materials[intersectedSphere.material_id-1].diffuse,cosAlphaDiffuse);
            color = add(color, add(hadamard(diffuse, receivedIrradiance), hadamard(specular, receivedIrradiance)));

            Ray reflectedRay{};
            reflectedRay.origin = P;
            reflectedRay.direction = normalize(add(ray.direction,multS(N,-2*dot(N,ray.direction))));
            auto reflectedColor = ComputeColor(reflectedRay, scene, depth-1);
            color = add(color, hadamard(reflectedColor, scene.materials[intersectedSphere.material_id-1].mirror));
        }
        color = add(color, hadamard(scene.materials[intersectedSphere.material_id-1].ambient, scene.ambient_light));
    }
    if (!intersectedMesh.faces.empty()) {
        auto P = add(add(ray.origin,multS(ray.direction,t)), {scene.shadow_ray_epsilon, scene.shadow_ray_epsilon, scene.shadow_ray_epsilon});
        auto N = normalize(NormalOfFace(intersectedMesh.faces[0], scene.vertex_data));
        auto W = normalize(add(ray.origin,multS(P,-1)));
        for (auto lights : scene.point_lights) {
            auto L = normalize(add(lights.position,multS(P,-1)));
            auto H = normalize(add(L,W));
            auto lightDistance = sqrt((lights.position.x-P.x)*(lights.position.x-P.x)+(lights.position.y-P.y)*(lights.position.y-P.y)+(lights.position.z-P.z)*(lights.position.z-P.z));
            auto receivedIrradiance = multS(lights.intensity,1/(lightDistance*lightDistance));
            auto cosAlphaSpecular = fmax(dot(N,H),0.0);
            auto cosAlphaDiffuse = fmax(dot(N,L),0.0);
            auto specular = multS(scene.materials[intersectedMesh.material_id-1].specular,pow(cosAlphaSpecular, scene.materials[intersectedMesh.material_id-1].phong_exponent));
            auto diffuse = multS(scene.materials[intersectedMesh.material_id-1].diffuse,(float) cosAlphaDiffuse);
            color = add(color, add(hadamard(diffuse, receivedIrradiance), hadamard(specular, receivedIrradiance)));

            Ray reflectedRay{};
            reflectedRay.origin = P;
            reflectedRay.direction = normalize(add(ray.direction,multS(N,-2*dot(N,ray.direction))));
            auto reflectedColor = ComputeColor(reflectedRay, scene, depth-1);
            color = add(color, hadamard(reflectedColor, scene.materials[intersectedMesh.material_id-1].mirror));
        }
        color = add(color, hadamard(scene.materials[intersectedMesh.material_id-1].ambient, scene.ambient_light));
    }
    return color;
}

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);
    auto spheres = scene.spheres;

    int width = scene.cameras[0].image_width;
    int height = scene.cameras[0].image_height;

    auto backroundColor = scene.background_color;

    auto* image = new unsigned char [width * height * 3]; // 3 channels

    for (int i = 0; i < width * height; ++i) //initialize image to backround color
    {
        image[i*3] = backroundColor.x;
        image[i*3+1] = backroundColor.y;
        image[i*3+2] = backroundColor.z;
    }


    for (int j = 0; j < height; j++) { //main loop
        for (int i = 0; i < width; i++) {
            Ray ray{};
            ray = GenerateRay(i, j, scene.cameras[0], width, height);
            auto color = ComputeColor(ray, scene, scene.max_recursion_depth);
            image[(j * width + i) * 3] = color.x > 255 ? 255 : round(color.x);
            image[(j * width + i) * 3 + 1] = color.y > 255 ? 255 : round(color.y);
            image[(j * width + i) * 3 + 2] = color.z > 255 ? 255 : round(color.z);
        }
    }

    write_ppm("test.ppm", image, width, height);

    delete[] image;

    return 0;
}
